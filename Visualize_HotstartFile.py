#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
process_hotstart.py
-------------------
Visualises the MRSWAT hotstart file (mr-hotstart.txt) as a longitudinal
flow-profile view of the lower Mississippi River.

Figure layout (two subplots, stacked):
    Top    — salinity  colored flow field, elevation on y-axis
    Bottom — velocity  colored flow field, elevation on y-axis

The flow domain is split into two layers (surface and bottom) whose
boundaries are the water surface, the layer interface, and the river bed.
Each layer is filled by its scalar value using pcolormesh on a regular
(river-mile × elevation) grid built via vectorised NumPy operations.
"""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
script_dir = Path(__file__).resolve().parent       # …/Code Storage/
master_dir = script_dir.parent                     # working directory
hotstart_path = master_dir / "mr-hotstart.txt"

# ---------------------------------------------------------------------------
# Load the hotstart matrix
# ---------------------------------------------------------------------------
# Each line contains 11 whitespace-delimited values.
data = np.loadtxt(hotstart_path, skiprows=1)

# Column mapping (0-indexed):
#  0  Time (elapsed days)
#  1  River mile
#  2  Water-surface elevation (ft)
#  3  Layer interface elevation (ft)
#  4  Bed elevation (ft)
#  5  Surface-layer salinity (ppt)
#  6  Bottom-layer salinity (ppt)
#  7  Surface-layer velocity (ft/s)
#  8  Bottom-layer velocity (ft/s)
#  9  Surface-layer discharge (cfs)
# 10  Bottom-layer discharge (cfs)

river_mile = data[:, 1]
wsel       = data[:, 2]   # water surface elevation (ft)
interface  = data[:, 3]   # layer interface elevation (ft)
bed        = data[:, 4]   # bed elevation (ft)
surf_salt  = data[:, 5]
bot_salt   = data[:, 6]
surf_vel   = data[:, 7]
bot_vel    = data[:, 8]

# ---------------------------------------------------------------------------
# Read hotstart metadata for the title (if available)
# ---------------------------------------------------------------------------
meta_path = master_dir / "hotstart_metadata.txt"
title_date = ""
if meta_path.exists():
    with open(meta_path, "r") as fh:
        for line in fh:
            if line.startswith("Hotstart actual date:"):
                title_date = line.split(":", 1)[1].strip()
                break

# ---------------------------------------------------------------------------
# Build 2-D scalar grids on a (elevation × river-mile) mesh
# ---------------------------------------------------------------------------
N_ELEV = 400   # vertical resolution of the colour grid

elev_min = bed.min()  - 2.0   # a bit below deepest bed
elev_max = wsel.max() + 0.5   # a bit above highest water surface
elevs = np.linspace(elev_min, elev_max, N_ELEV)

# Broadcast arrays for vectorised layer assignment
# E  shape (N_ELEV, 1),  others shape (1, n_x)
E   = elevs[:, np.newaxis]
B   = bed[np.newaxis, :]
I   = interface[np.newaxis, :]
W   = wsel[np.newaxis, :]

in_bottom  = (E >= B) & (E <= I)
in_surface = (E >  I) & (E <= W)

# Salinity grid  —  NaN outside the flow domain
salt_grid = np.where(
    in_surface, surf_salt[np.newaxis, :],
    np.where(in_bottom, bot_salt[np.newaxis, :], np.nan),
)
salt_grid = np.ma.masked_invalid(salt_grid)

# Velocity grid  —  NaN outside the flow domain
vel_grid = np.where(
    in_surface, surf_vel[np.newaxis, :],
    np.where(in_bottom, bot_vel[np.newaxis, :], np.nan),
)
vel_grid = np.ma.masked_invalid(vel_grid)

# Coordinate grids for pcolormesh  (shape: N_ELEV × n_x)
X, Y = np.meshgrid(river_mile, elevs)

# ---------------------------------------------------------------------------
# Velocity colour scale — symmetric diverging around zero
# ---------------------------------------------------------------------------
vel_abs_max = np.nanmax(np.abs(vel_grid))

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
fig, (ax_salt, ax_vel) = plt.subplots(
    2, 1, figsize=(16, 9), sharex=True,
)

suptitle = "MRSWAT Hotstart — Longitudinal Flow Profile"
if title_date:
    suptitle += f"  ({title_date})"
fig.suptitle(suptitle, fontsize=14, fontweight="bold")

# ── Helper: draw bed fill and boundary lines on an axis ──────────────────
def _add_boundaries(ax):
    """Fill the riverbed and overlay water-surface / interface lines."""
    ax.fill_between(river_mile, elev_min, bed,
                    color="#8B7355", zorder=3, label="Riverbed")
    ax.plot(river_mile, wsel,       color="steelblue",
            lw=1.0, zorder=4, label="Water surface")
    ax.plot(river_mile, interface,  color="white",
            lw=0.8, ls="--", zorder=4, label="Layer interface")
    ax.set_ylabel("Elevation (ft)", fontsize=10)
    ax.set_ylim(elev_min, elev_max)
    ax.grid(True, linestyle="--", alpha=0.25, zorder=0)

# ── Salinity ─────────────────────────────────────────────────────────────
im_salt = ax_salt.pcolormesh(
    X, Y, salt_grid,
    cmap="YlOrRd", shading="nearest", zorder=2,
)
_add_boundaries(ax_salt)
cb_salt = fig.colorbar(im_salt, ax=ax_salt, pad=0.01, fraction=0.025)
cb_salt.set_label("Salinity (ppt)", fontsize=9)
ax_salt.set_title("Salinity", fontsize=11, loc="left")
ax_salt.legend(loc="upper right", fontsize=8, framealpha=0.85)

# ── Velocity ─────────────────────────────────────────────────────────────
im_vel = ax_vel.pcolormesh(
    X, Y, vel_grid,
    cmap="RdBu_r", shading="nearest", zorder=2,
    vmin=-vel_abs_max, vmax=vel_abs_max,
)
_add_boundaries(ax_vel)
cb_vel = fig.colorbar(im_vel, ax=ax_vel, pad=0.01, fraction=0.025)
cb_vel.set_label("Velocity (ft/s)", fontsize=9)
ax_vel.set_title("Velocity", fontsize=11, loc="left")
ax_vel.set_xlabel("River Mile", fontsize=10)
ax_vel.legend(loc="upper right", fontsize=8, framealpha=0.85)

fig.subplots_adjust(top=0.93, hspace=0.18, right=0.95)

# ---------------------------------------------------------------------------
# Save and display
# ---------------------------------------------------------------------------
out_path = master_dir / "Output" / "hotstart_profile.png"
out_path.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(out_path, dpi=150, bbox_inches="tight")
print(f"  ✓ Figure saved: {out_path}")
plt.show()
