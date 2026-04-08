"""
Module: generate_tiles.py
Functions to create, reproject, export and index subtiles of the forceCube main tiles
Autor: Paul Magdon
"""

from shapely.geometry import box, polygon
from shapely.wkb import loads as load_wkb
from pyproj import Transformer
from shapely.ops import transform
from shapely.strtree import STRtree
from dataclasses import dataclass
from typing import List
from osgeo import ogr, osr
import os

# ====== Dataklasse für Tiles ======

@dataclass
class Tile:
    idx_x: int
    idx_y: int
    polygon: box

# ====== Tile-Generierung ======

def generate_target_tiles(
    origin_x: float,
    origin_y: float,
    count_x: int,
    count_y: int,
    size: float = 1000.0
    ) -> List[Tile]:
    """
    Erstellt ein Raster aus Polygon-Kacheln ab origin_x/origin_y.
    
    Args:
        origin_x: links-unten X (float, im Ziel-SRS)
        origin_y: links-unten Y (float, im Ziel-SRS)
        count_x: Zahl Kacheln in X-Richtung
        count_y: Zahl Kacheln in Y-Richtung
        size: Kantenlänge jeder Kachel (Meter)
    Returns:
        List[Tile] mit Attributen idx_x, idx_y, polygon (Shapely)
    """
    tiles = []
    for ix in range(count_x):
        for iy in range(count_y):
            minx = origin_x + ix * size
            miny = origin_y + iy * size
            maxx = minx + size
            maxy = miny + size
            tiles.append(Tile(idx_x=ix, idx_y=iy, polygon=box(minx, miny, maxx, maxy)))
    return tiles

# ====== Reprojection =======

def reproject_tiles(tiles: List[Tile], from_epsg: int, to_epsg: int) -> List[Tile]:
    transformer = Transformer.from_crs(from_epsg, to_epsg, always_xy=True)
    return [
        Tile(
            idx_x=tile.idx_x,
            idx_y=tile.idx_y,
            polygon=transform(transformer.transform, tile.polygon)
        )
        for tile in tiles
    ]

# ====== Export ======

def export_tiles_ogr(
        tiles: List[Tile],
        out_gpkg: str,
        layername: str = 'tiles',
        epsg: int = 3035
    ):
    # Existierende Datei ggf. löschen (Achtung!)
    if os.path.exists(out_gpkg):
        os.remove(out_gpkg)
    # Treiber wählen 
    drv = ogr.GetDriverByName('GPKG') 
    ds = drv.CreateDataSource(out_gpkg)
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(epsg)
    layer = ds.CreateLayer(layername, srs=srs, geom_type=ogr.wkbPolygon)
    # Attribute für Indizes
    layer.CreateField(ogr.FieldDefn("idx_x", ogr.OFTInteger))
    layer.CreateField(ogr.FieldDefn("idx_y", ogr.OFTInteger))
    # Tiles hinzufügen:
    for tile in tiles:
        ring = ogr.Geometry(ogr.wkbLinearRing)
        coords = list(tile.polygon.exterior.coords)
        for x, y in coords:
            ring.AddPoint(float(x), float(y))
        ogr_poly = ogr.Geometry(ogr.wkbPolygon)
        ogr_poly.AddGeometry(ring)
        feat = ogr.Feature(layer.GetLayerDefn())
        feat.SetGeometry(ogr_poly)
        feat.SetField('idx_x', tile.idx_x)
        feat.SetField('idx_y', tile.idx_y)
        layer.CreateFeature(feat)
        feat = None       # Speicher freigeben
        ogr_poly = None
        ring = None
    ds = None  # Datei flushen/schließen

# ====== Spatial Index ======

def build_spatial_index(tiles: List[Tile]):
    polys = [t.polygon for t in tiles]
    index = STRtree(polys)
    mapping = {id(poly): t for t, poly in zip(tiles, polys)}
    return index, mapping

# ====== Beispiel-Nutzung ======
if __name__ == "__main__":

    
    # Import of force tile
    forceTilebb ="/mnt/data/raw/th/project_tests/fc_x0065_y0051.gpkg"

    # Extract origin of force tile
    ds = ogr.Open(forceTilebb)
    layer = ds.GetLayer(0)  
    feature = layer.GetNextFeature()
    geom = feature.GetGeometryRef()
    wkb_geom = geom.ExportToWkb(ogr.wkbNDR) 
    polygon = load_wkb(bytes(wkb_geom))  # Shapely-Objekt
    minx, miny, maxx, maxy = polygon.bounds
    print(f"x_orig = {minx}, y_orig = {miny}")  
    # Parameter
    origin_x, origin_y = minx, miny
    n_x, n_y = 30, 30
    size = 1000.0
    epsg = 3035

    # Tiles generieren
    tiles = generate_target_tiles(origin_x, origin_y, n_x, n_y, size)
    # Export als GPKG
    out_gpkg = "tiles_3035.gpkg"
    export_tiles_ogr(tiles, out_gpkg, layername="tiles", epsg=epsg)
    # Für Reprojektion z.B.:
    tiles_utm = reproject_tiles(tiles, 3035, 25832)
    export_tiles_ogr(tiles_utm, "tiles_utm32.gpkg", layername="tiles", epsg=25832)