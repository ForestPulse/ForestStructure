#!/usr/bin/env python

import argparse
import geopandas as gpd
from pathlib import Path

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--tindex_gpkg", required=True)
    p.add_argument("--tindex_layer", required=True)
    p.add_argument("--tindex_id_field", required=True)
    p.add_argument("--force_gpkg", required=True)
    p.add_argument("--force_layer", required=True)
    p.add_argument("--force_id_field", required=True)
    p.add_argument("--target_id", required=True)
    p.add_argument("--output", required=True)
    a = p.parse_args()

    tindex = gpd.read_file(a.tindex_gpkg, layer=a.tindex_layer)
    force  = gpd.read_file(a.force_gpkg,  layer=a.force_layer)

    if tindex.crs != force.crs:
        force = force.to_crs(tindex.crs)

    fc = force[force[a.force_id_field] == a.target_id]
    if fc.empty:
        raise SystemExit(f"Target {a.target_id} not found.")
    geom = fc.geometry.iloc[0]

    cand = tindex.iloc[list(tindex.sindex.intersection(geom.bounds))]
    hits = cand[cand.intersects(geom)][a.tindex_id_field].astype(str)

    Path(a.output).write_text("\n".join(hits), encoding="utf-8")
    print(f"{len(hits)} tiles intersect {a.target_id}")
    print(f"Written: {a.output}")

if __name__ == "__main__":
    main()
