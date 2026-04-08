#!/usr/bin/env python3
import os
import json
import time
import pickle
import argparse
import numpy as np
import laspy
import pdal
from concurrent.futures import ProcessPoolExecutor, as_completed
from pyproj import Transformer
from shapely.geometry import box
from shapely.ops import transform as shapely_transform
from shapely.strtree import STRtree
from osgeo import ogr
from shapely.wkb import loads as load_wkb

# --- Configuration & Constants ---
CHUNK_SIZE = 5_000_000  # points per chunk

# Mapping of laspy dimension names to PDAL-compatible names and types
PDAL_DIMENSION_MAP = {
    'intensity': ('Intensity', np.uint16),
    'return_number': ('ReturnNumber', np.uint8),
    'number_of_returns': ('NumberOfReturns', np.uint8),
    'classification': ('Classification', np.uint8),
    'scan_angle_rank': ('ScanAngleRank', np.float32),
    'user_data': ('UserData', np.uint8),
    'point_source_id': ('PointSourceId', np.uint16),
    'gps_time': ('GpsTime', np.float64),
    'red': ('Red', np.uint16),
    'green': ('Green', np.uint16),
    'blue': ('Blue', np.uint16),
}

def process_single_subtile(tile_data, input_epsg, output_epsg, out_dir, force_cube_id):
    """
    Worker function to process a single 1km x 1km subtile.
    """
    poly = tile_data['polygon']
    relevant_files = tile_data['files']

    minx, miny, maxx, maxy = poly.bounds
    out_path = os.path.join(out_dir, f"las_{round(minx)}_{round(miny)}_{force_cube_id.lower()}.copc.laz")

    if not relevant_files:
        return f"Skipped {tile_data['idx']}: No source files."

    transformer = Transformer.from_crs(f"EPSG:{input_epsg}",
                                       f"EPSG:{output_epsg}",
                                       always_xy=True)
    all_chunks =[]
    
    # Lock the schema on the first valid chunk to prevent np.concatenate crashes
    dtype = None
    active_dims =[]

    try:
        for f in relevant_files:
            with laspy.open(f) as lazf:
                for chunk in lazf.chunk_iterator(CHUNK_SIZE):
                    # 1) Reproject XY
                    tx, ty = transformer.transform(chunk.x, chunk.y)

                    # 2) Fast bbox mask in target CRS (EPSG:3035)
                    mask = (tx >= minx) & (tx <= maxx) & (ty >= miny) & (ty <= maxy)
                    if not np.any(mask):
                        continue

                    # Determine dtype dynamically on the FIRST chunk containing points
                    if dtype is None:
                        dtype =[('X', 'f8'), ('Y', 'f8'), ('Z', 'f8')]
                        for d in PDAL_DIMENSION_MAP.keys():
                            if hasattr(chunk, d):
                                dim_name, dim_type = PDAL_DIMENSION_MAP[d]
                                dtype.append((dim_name, dim_type))
                                active_dims.append((d, dim_name))

                    # 3) Build structured array
                    arr = np.empty(np.sum(mask), dtype=dtype)
                    arr['X'], arr['Y'] = tx[mask], ty[mask]
                    arr['Z'] = chunk.z[mask]

                    for d, dim_name in active_dims:
                        # getattr with a default zero array prevents errors if a subsequent 
                        # file in the same subtile is missing a dimension (e.g. colors)
                        default_zeros = np.zeros_like(chunk.z)
                        arr[dim_name] = getattr(chunk, d, default_zeros)[mask]

                    all_chunks.append(arr)

        if not all_chunks:
            return f"Empty {tile_data['idx']}"

        full_array = np.concatenate(all_chunks)

        # 4) PDAL pipeline
        pipeline_json =[
            {"type": "filters.outlier", "method": "statistical", "mean_k": 8, "multiplier": 2.5},
            {"type": "filters.hag_delaunay", "allow_extrapolation": True},
            {"type": "filters.range", "limits": "HeightAboveGround[0:100],Classification![7:7]"},
            {"type": "filters.relaxationdartthrowing", "count": 10000000},
            {
                "type": "writers.copc",
                "filename": out_path,
                "forward": "all",
                "a_srs": f"EPSG:{output_epsg}+7837",
                "extra_dims": "HeightAboveGround=float32",
                "offset_x": "auto",
                "offset_y": "auto",
                "offset_z": "auto", 
                "scale_x": 0.001,
                "scale_y": 0.001,
                "scale_z": 0.001     
            }
        ]

        pipeline = pdal.Pipeline(json.dumps(pipeline_json), arrays=[full_array])
        pipeline.execute()
        return f"Success {tile_data['idx']}: {out_path}"

    except Exception as e:
        return f"Error {tile_data['idx']}: {str(e)}"


def main():
    parser = argparse.ArgumentParser(
        description="Parallel PDAL reprojection for FORCE 30km cube into 1km subtiles"
    )
    parser.add_argument("--force_gpkg", required=True,
                        help="GPKG with FORCE cubes (EPSG:25832 or EPSG:3035)")
    parser.add_argument("--force_layer", default=None,
                        help="Layer name (default: first layer)")
    parser.add_argument("--force_cube_id", required=True,
                        help="ID of FORCE cube to process, e.g. X0065_Y0051")

    parser.add_argument("--input_laz_dir", required=True,
                        help="Directory containing source LAZ tiles (EPSG:25832)")
    parser.add_argument("--output_laz_dir", required=True,
                        help="Directory for output COPC LAZ files")
    parser.add_argument("--header_cache", required=True,
                        help="Pickle header cache for this force cube")

    parser.add_argument("--start_id", type=int, default=0,
                        help="Start subtile index [0..899]")
    parser.add_argument("--end_id", type=int, default=899,
                        help="End subtile index [0..899]")
    parser.add_argument("--num_workers", type=int, default=1,
                        help="Number of parallel processes")
    parser.add_argument("--src_epsg", type=int, default=25832)
    parser.add_argument("--dst_epsg", type=int, default=3035)

    args = parser.parse_args()

    # --- Step 1: read target FORCE cube (single feature) ---
    ds = ogr.Open(args.force_gpkg)
    if ds is None:
        raise RuntimeError(f"Cannot open {args.force_gpkg}")

    if args.force_layer is None:
        layer = ds.GetLayer(0)
    else:
        layer = ds.GetLayerByName(args.force_layer)
        if layer is None:
            raise RuntimeError(f"Layer {args.force_layer} not found in {args.force_gpkg}")

    layer.SetAttributeFilter(f"id = '{args.force_cube_id}'")
    feature = layer.GetNextFeature()
    if feature is None:
        raise RuntimeError(f"Cube id {args.force_cube_id} not found in {args.force_gpkg}")

    # Read native geometry from the Geopackage
    polygon_raw = load_wkb(bytes(feature.GetGeometryRef().ExportToWkb(ogr.wkbNDR)))

    # Dynamically read the EPSG of the Geopackage
    spatial_ref = layer.GetSpatialRef()
    gpkg_epsg = args.src_epsg  # Default fallback if unknown
    if spatial_ref is not None:
        auth_code = spatial_ref.GetAuthorityCode(None)
        if auth_code:
            gpkg_epsg = int(auth_code)

    # Force project the GPKG geometry to the target EPSG (3035) to get perfectly square FORCE cubes
    if gpkg_epsg != args.dst_epsg:
        print(f"Projecting FORCE cube {args.force_cube_id} geometry from EPSG:{gpkg_epsg} to EPSG:{args.dst_epsg}...")
        gpkg_to_dst = Transformer.from_crs(f"EPSG:{gpkg_epsg}", f"EPSG:{args.dst_epsg}", always_xy=True).transform
        polygon = shapely_transform(gpkg_to_dst, polygon_raw)
    else:
        polygon = polygon_raw

    # Perfectly align bounds
    minx, miny, maxx, maxy = polygon.bounds

    print(f"FORCE cube {args.force_cube_id} bounds (Aligned to EPSG:{args.dst_epsg}): "
          f"{(minx, miny, maxx, maxy)}")

    # --- Step 2: subdivide 30 km tile into 900×1 km subtiles ---
    subtiles =[]
    for i in range(30):
        for j in range(30):
            subtiles.append({
                "idx": (i, j),
                "polygon": box(minx + i*1000, miny + j*1000,
                               minx + (i+1)*1000, miny + (j+1)*1000),
            })

    # --- Step 3: load header cache ---
    with open(args.header_cache, "rb") as f:
        cache = pickle.load(f)
    raw_headers = cache["tiles"] if isinstance(cache, dict) else cache

    input_polys = [h["polygon"] for h in raw_headers]
    tree = STRtree(input_polys)

    # Transformer to temporarily project the 1km subtiles back to src (25832) to query the spatial index
    project_to_src = Transformer.from_crs(
        f"EPSG:{args.dst_epsg}", f"EPSG:{args.src_epsg}", always_xy=True
    ).transform

    # --- Step 4: build work queue (map subtiles to relevant source files) ---
    work_queue =[]
    for i in range(args.start_id, args.end_id + 1):
        if i >= len(subtiles):
            break
        st = subtiles[i]

        # Project 3035 subtile polygon back to EPSG:25832 specifically for intersecting with LAZ headers
        st_src = shapely_transform(project_to_src, st["polygon"])
        
        indices = tree.query(st_src)
        relevant_files = [
            raw_headers[idx]["filename"]
            for idx in indices
            if st_src.intersects(input_polys[idx])
        ]
        if relevant_files:
            st["files"] = relevant_files
            work_queue.append(st)

    print(f"Node processing {len(work_queue)} subtiles with {args.num_workers} workers.")

    os.makedirs(args.output_laz_dir, exist_ok=True)

    # --- Step 5: run workers ---
    with ProcessPoolExecutor(max_workers=args.num_workers) as executor:
        futures = {
            executor.submit(process_single_subtile,
                            item,
                            args.src_epsg,
                            args.dst_epsg,
                            args.output_laz_dir,
                            args.force_cube_id): item
            for item in work_queue
        }
        for future in as_completed(futures):
            print(future.result())

if __name__ == "__main__":
    main()
