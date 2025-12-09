## PDAL Tindex and ForceCube Tile Intersection Workflow
This documentation describes the workflow for generating a PDAL
**tindex** (tile index) in **EPSG:3035** from LiDAR LAZ tiles and
performing spatial intersection with **ForceCube** grid tiles to
determine which LiDAR tiles fall inside a selected ForceCube cell.

### 1. Generating a Tindex

The input LAZ tiles are stored in **EPSG:25832** (ETRS89 / UTM 32N).\
The resulting GeoPackage tindex is written in **EPSG:3035**, ensuring
compatibility with European equal-area grid systems and seamless
integration with ForceCube.

### Command

``` bash
pdal tindex create " /path/to/write/Thu_2020_25_tindex_epsg3035.gpkg" \
                    "/input/LAZtiles/th_data/Thuringia_2020-2025_tiles/*.laz" \
                    --ogrdriver GPKG \
                    --tindex_name tile \
                    --fast_boundary \
                    --lyr_name "Thu_2020_25_tindex_epsg3035" \
                    --a_srs "EPSG:25832" \
                    --t_srs "EPSG:3035"
```

  `--ogrdriver GPKG`          Writes the tile index to an OGR GeoPackage
                              file.

  `--tindex_name tile`        Creates a column named **tile**, containing
                              the full path or name of each LAZ file.

  `--lyr_name`                Sets the name of the vector layer inside
                              the GeoPackage.

  `--fast_boundary`           Uses fast bounding boxes rather than exact
                              tile hulls (recommended for many-tile
                              workflows).

  `--a_srs EPSG:25832`        Assigns the CRS of the input LAZ data
                              (ETRS89 / UTM 32N).

  `--t_srs EPSG:3035`         Reprojects the tile polygons to the target
                              CRS of the output GPKG.
                              
------------------------------------------------------------------------

### 2. Intersecting Tindex Tiles with ForceCube Grid Cells
After generating the tindex, the next step is to identify all LiDAR
tiles that fall inside a selected ForceCube tile (e.g., `X0065_Y0051`).

### Example Command
``` bash
/path/of/this/file/tilesClip.py \
      --tindex_gpkg "D:\ForceCube\Thu_2014_19_tindex_epsg3035.gpkg" \
      --tindex_layer "Thu_2014_19_tindex_epsg3035" \
      --tindex_id_field "tile" \
      --force_gpkg "D:\ForceCube\DE_force-cube.gpkg" \
      --force_layer "epsg3035" \
      --force_id_field "id" \
      --target_id "X0065_Y0051" \
      --output "D:\ForceCube\tiles_in_X0065_Y0051.txt"
```

### Workflow Summary

    Tindex (EPSG:3035)        ForceCube (EPSG:3035)
            ▼                          ▼
       tilesClip.py  → geometric intersection → list of intersecting LiDAR tiles

Both datasets in **EPSG:3035** ensures accurate, CRS-consistent spatial
operations.



