#!/usr/bin/env python3
"""
Generate a spatial header cache (.pkl) for LAZ/LAS tiles intersecting a FORCE cube,
using a CSV with columns: id,left,right,top,bottom.
"""

import os
import sys
import pickle
import argparse
import time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

import laspy
from shapely.geometry import box


def build_force_cube_polygon_from_csv(force_csv, force_cube_id):
    import csv
    with open(force_csv, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row["id"].strip() == force_cube_id:
                left   = float(row["left"])
                right  = float(row["right"])
                top    = float(row["top"])
                bottom = float(row["bottom"])
                return box(left, bottom, right, top)  # (minx, miny, maxx, maxy)
    sys.exit(f"ERROR: Force cube '{force_cube_id}' not found in {force_csv}")


def find_laz_files(input_dir):
    patterns = ["**/*.laz", "**/*.las", "**/*.LAZ", "**/*.LAS"]
    all_files = []
    for pat in patterns:
        all_files.extend(Path(input_dir).glob(pat))
    return sorted(set(all_files))


def read_header(filepath):
    filepath = str(filepath)
    try:
        with laspy.open(filepath) as lf:
            hdr = lf.header
            poly = box(hdr.mins[0], hdr.mins[1], hdr.maxs[0], hdr.maxs[1])
            return {
                "filename": os.path.abspath(filepath),
                "polygon": poly,
                "minx": float(hdr.mins[0]),
                "miny": float(hdr.mins[1]),
                "maxx": float(hdr.maxs[0]),
                "maxy": float(hdr.maxs[1]),
                "minz": float(hdr.mins[2]),
                "maxz": float(hdr.maxs[2]),
                "point_count": hdr.point_count,
            }
    except Exception as e:
        return {"filename": filepath, "error": str(e)}


def main():
    parser = argparse.ArgumentParser(
        description="Generate header cache for LAZ/LAS tiles intersecting a FORCE cube (CSV definition)"
    )
    parser.add_argument("--input", required=True,
                        help="Root directory with .laz/.las files (recursive)")
    parser.add_argument("--force_csv", required=True,
                        help="CSV with cube extents: id,left,right,top,bottom")
    parser.add_argument("--force_cube_id", required=True,
                        help="Target cube ID, e.g. X0065_Y0051")
    parser.add_argument("--output", required=True,
                        help="Output dir or .pkl path. Dir → header_cache_<id>.pkl")
    parser.add_argument("--workers", type=int, default=8,
                        help="Parallel workers (default: 8)")
    args = parser.parse_args()

    t0 = time.time()
    cube_poly = build_force_cube_polygon_from_csv(args.force_csv, args.force_cube_id)
    print(f"[1/3] Cube {args.force_cube_id} bounds: {cube_poly.bounds}")

    laz_files = find_laz_files(args.input)
    print(f"[2/3] Found {len(laz_files)} LAZ/LAS files under {args.input}")

    results = []
    n_workers = min(args.workers, max(len(laz_files), 1))
    with ProcessPoolExecutor(max_workers=n_workers) as ex:
        futures = {ex.submit(read_header, f): f for f in laz_files}
        for i, fut in enumerate(as_completed(futures), 1):
            r = fut.result()
            if "polygon" in r and r["polygon"].intersects(cube_poly):
                results.append(r)
            if i % 200 == 0 or i == len(futures):
                print(f"    checked {i}/{len(futures)} headers...")

    errors = [r for r in results if "error" in r]
    successes = [r for r in results if "error" not in r]

    if errors:
        print(f"WARNING: {len(errors)} files had read errors (showing first 5):")
        for e in errors[:5]:
            print(f"  {e['filename']}: {e['error']}")

    out_path = args.output
    if os.path.isdir(out_path) or not out_path.endswith(".pkl"):
        os.makedirs(out_path, exist_ok=True)
        out_path = os.path.join(out_path, f"header_cache_{args.force_cube_id}.pkl")

    cache = {
        "force_cube_id": args.force_cube_id,
        "force_cube_bounds": cube_poly.bounds,
        "n_tiles": len(successes),
        "tiles": successes,
    }

    with open(out_path, "wb") as p:
        pickle.dump(cache, p, protocol=pickle.HIGHEST_PROTOCOL)

    print(f"[3/3] Saved {len(successes)} tile headers to {out_path}")
    print(f"       Elapsed: {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
