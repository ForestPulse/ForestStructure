# NFI Buffersampling | LF-RLP 2025

#-------Setup----------|

import geopandas as gpd
import pandas as pd
import numpy as np
from rasterio.mask import mask as rio_mask
from shapely.strtree import STRtree
from shapely.geometry import box
import rasterio
import math
import subprocess
from pathlib import Path
from tqdm import tqdm
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# | standard deviation of masurements for trees

SIGMA_GNSS  =  5.0                        # X,Y coordinates
SIGMA_S_REL = 0.0025 / 10                 # [m/m] distance: 2,5 mm per 10 m 
SIGMA_ALPHA = 2 * math.pi / 200           # [rad] angle: 2 gon 

# | input

trees = gpd.read_file("")                # Input nfi-file in gpkg format
height_col = "M_Hoe"                     # declare name of height specific column

tile_dir = Path("")                      # directory of your raster-tiles
YEARS = {"2021", "2022"}                 # years to sample from should match the aquisition time of nfi data (BWI4=2021/2022)
tile_paths = [
    p for p in tile_dir.glob("*.tif")
    if p.stem.split("_")[-1] in YEARS      
]

#--------functions-------|

def assign_layers(group):

    h = group[height_col]

    # | fewer than 3 trees --> main layer
    if h.notna().sum() < 3:
        group["Schicht"] = "Hauptschicht" # main-layer
        return group

    # | quantile per plot

    plot_max=  np.max(h)

    def classify(val):
        if pd.isna(val):
            return np.nan
        elif val <= (1/3)*plot_max :
            return "Unterschicht"     # lower-layer
        elif val <= (2/3)*plot_max:
            return "Zwischenschicht"  # intermediate-layer
        else:
            return "Hauptschicht"     # main-layer

    group["Schicht"] = h.apply(classify)

    return group

#-------preprocessing--------|

trees[height_col] = pd.to_numeric(trees[height_col], errors="coerce")
trees["M_Hoe"] = trees["M_Hoe"] / 10

# | only living trees in the main layer (other filters can optionally be applied)

trees = trees.groupby("tnr_enr", group_keys=False).apply(assign_layers)

trees["DIST"] = trees["Hori"] / 100

trees_valid = trees[trees["Pk"] == 1 | 0 ]                            # 0 = new tree, 1 = tree in previous nfi
trees_valid = trees_valid[trees_valid["tot"] == 0]                    # attribute for dead trees
trees_valid = trees_valid[trees_valid["Schicht"] == 'Hauptschicht']

trees_valid.head(n = 50)

# ---------buffer generation---------|

# | calculate individual tree positions

trees_valid["azimut_rad"] = trees_valid["Azi"] * np.pi / 200

trees_valid["tree_x"] = trees_valid["soll_x_32n"] + trees_valid["DIST"] * np.sin(trees_valid["azimut_rad"]) # fill in X-coordinate, distance and angle(rad)
trees_valid["tree_y"] = trees_valid["soll_y_32n"] + trees_valid["DIST"] * np.cos(trees_valid["azimut_rad"]) # fill in Y-coordinate, distance and angle(rad)

trees_valid["geometry"] = gpd.points_from_xy(trees_valid.tree_x, trees_valid.tree_y)

trees_points = gpd.GeoDataFrame(trees_valid, geometry="geometry", crs = trees.crs)

s = trees_valid["DIST"]                  # distance [m]
a = trees_valid["azimut_rad"]            # angle [rad]

sigma_s = (s * SIGMA_S_REL) + 0.005                # absolute distance error [m]

#  | partial deviation for new coordinates 
#  | x_new = x_st + s·sin(α)  →  ∂x/∂x_st=1, ∂x/∂s=sin(α), ∂x/∂α= s·cos(α)
#  | y_new = y_st + s·cos(α)  →  ∂y/∂y_st=1, ∂y/∂s=cos(α), ∂y/∂α=-s·sin(α)

sigma_x2 = (SIGMA_GNSS**2
            + np.sin(a)**2 * sigma_s**2
            + s**2 * np.cos(a)**2 * SIGMA_ALPHA**2)

sigma_y2 = (SIGMA_GNSS**2
            + np.cos(a)**2 * sigma_s**2
            + s**2 * np.sin(a)**2 * SIGMA_ALPHA**2)

#  | mean point error (isotropic buffer radius)
#  | σ_x² + σ_y² = 2·σ_GNSS² + σ_s² + s²·σ_α²
trees_points["sigma_pos"] = np.sqrt(sigma_x2 + sigma_y2)

#  | save
trees_points["sigma_gnss_val"] = np.sqrt(2 * SIGMA_GNSS**2)
trees_points["sigma_dist_val"] = sigma_s
trees_points["sigma_angle_val"] = s * SIGMA_ALPHA

#  | buffer with r = σ_pos
trees_buffer = trees_points.copy()
trees_buffer["geometry"] = trees_points.geometry.buffer(trees_points["sigma_pos"])

#-----------sampling-----------| 

max_vals = []

for _, row in tqdm(trees_buffer.iterrows(), total=len(trees_buffer)):
    geom = row.geometry
    
   
    candidate_idxs = tree.query(geom, predicate="intersects")
    
    pixel_max = np.nan
    
    for ci in candidate_idxs:
        with rasterio.open(tile_paths[ci]) as src:
            try:
                
                out_image, _ = rio_mask(
                    src,
                    [geom],
                    crop=True,
                    nodata=np.nan,
                    all_touched=False  # just pixels with centroid in buffer
                )
                tile_max = np.nanmax(out_image)
                
                if np.isnan(pixel_max) or tile_max > pixel_max:
                    pixel_max = tile_max
            except Exception:
                pass  
    
    max_vals.append(pixel_max)

trees_buffer["nDOM_max"] = max_vals

#----------outlier removal-----------|

# | MAE

df = trees_buffer[["M_Hoe", "nDOM_max"]].dropna().copy()
df["diff"] = df["M_Hoe"] - df["nDOM_max"]        
df["abs_diff"] = df["diff"].abs()

n_before = len(df)

# | Z-Score for outlier removal

Z_THRESHOLD = 2

z_diff = np.abs(stats.zscore(df["diff"]))        
mask_clean = z_diff < Z_THRESHOLD
df_clean = df[mask_clean].copy()

n_after   = len(df_clean)
n_removed = n_before - n_after
print(f"objects before outlier removal:          {n_before}")
print(f"objects removed (|z_diff| ≥ {Z_THRESHOLD}): {n_removed}  ({n_removed/n_before*100:.1f} %)")
print(f"objects before outlier removal:         {n_after}")

#  Stats for clean data

std_diff  = df_clean["abs_diff"].std()
mean_diff = df_clean["abs_diff"].mean()
print(f"\nabsolute difference |M_Hoe – nDOM_max|:")
print(f"  Mean : {mean_diff:.3f} m")
print(f"  Std.-dev.  : {std_diff:.3f} m")

r, p_value = stats.pearsonr(df_clean["M_Hoe"], df_clean["nDOM_max"])
print(f"\nPearson-Correlation r = {r:.4f}  (p = {p_value:.2e})")

# ---------plotting----------|

fig, ax = plt.subplots(figsize=(7, 7))


sc = ax.scatter(
    df_clean["M_Hoe"], df_clean["nDOM_max"],
    c=df_clean["abs_diff"],
    cmap="RdYlGn_r",
    alpha=0.7,
    edgecolors="white",
    linewidths=0.3,
    s=40,
    zorder=3
)
plt.colorbar(sc, ax=ax, label="Absolute difference [m]")

# 1:1-Linie
lim_min = min(df_clean["M_Hoe"].min(), df_clean["nDOM_max"].min()) - 1
lim_max = max(df_clean["M_Hoe"].max(), df_clean["nDOM_max"].max()) + 1
ax.plot([lim_min, lim_max], [lim_min, lim_max],
        color="gray", linestyle="--", linewidth=1, zorder=2)

# Annotations
ax.text(0.05, 0.95,
        f"r = {r:.3f}\np = {n_after}\n"
        f"Ø |Δ| = {mean_diff:.2f} m\nσ |Δ| = {std_diff:.2f} m\n"
        f"Outlier removed: {n_removed}",
        transform=ax.transAxes, fontsize=9,
        verticalalignment="top",
        bbox=dict(boxstyle="round,pad=0.4", facecolor="white", alpha=0.8))

ax.set_xlim(lim_min, lim_max)
ax.set_ylim(lim_min, lim_max)
ax.set_xlabel("measured height M_Hoe [m]", fontsize=11)
ax.set_ylabel("nDOM Max-value [m]", fontsize=11)
ax.set_title("measured vs. nDOM", fontsize=13, fontweight="bold")
ax.legend(fontsize=9)
ax.set_aspect("equal")
ax.grid(True, linestyle=":", alpha=0.5)

plt.tight_layout()
plt.savefig("scatterplot_hoehe_chm.png", dpi=150, bbox_inches="tight")
plt.show()

