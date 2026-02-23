"""
Generate a synthetic hot-start (restart) file for the MRSWAT model.

Reads the first time-block of mr-full-output.txt to obtain the spatial
grid (river miles) and bed elevations, then constructs a synthetic
initial-condition matrix with a user-defined intrusion-toe position.

Columns in the output hot-start file match the full-output layout:
  1  Time                        → 9.9999  (synthetic flag)
  2  River-mile                  → from template
  3  Wsel (ft)                   → 1.06 (default surface-water elev.)
  4  Layer-Interface-el (ft)     → computed
  5  Bed-el (ft)                 → from template
  6  Surface-salt (ppt)          → computed
  7  Bottom-salt (ppt)           → computed
  8  Surface-vel (ft/sec)        → 0.0
  9  Bottom-vel (ft/sec)         → 0.0
 10  Surface-Q (cfs)             → 0.0
 11  Bottom-Q (cfs)              → 0.0
"""

import numpy as np
import sys
import os

# ── file paths ───────────────────────────────────────────────────────
full_output_file = "mr-full-output.txt"
hotstart_file    = "mr-hotstart.txt"

# ── constants ────────────────────────────────────────────────────────
SYNTHETIC_TIME = 9.9999          # flag value for all rows
DEFAULT_WSEL   = 1.06            # default surface-water elevation (ft)

# =====================================================================
# 1.  Parse the FIRST time-block from mr-full-output.txt
# =====================================================================
print(f"Reading template from {full_output_file} ...")

river_miles = []
bed_els     = []

with open(full_output_file, "r") as fh:
    found_first_block = False

    for raw_line in fh:
        line = raw_line.strip()
        if not line:
            continue

        # Detect "Time = ..." header lines
        if line.startswith("Time"):
            if "=" in line:
                if found_first_block:
                    # We've hit the second time block → done
                    break
                found_first_block = True
                continue
            else:
                # Part of the column-header text → skip
                continue

        # Skip remaining header lines (non-numeric leading words)
        try:
            float(line.split()[0])
        except (ValueError, IndexError):
            continue

        # Each data row is a single line with 11 whitespace-separated values:
        #  time  river-mile  wsel  interface  bed  Ssurf  Sbot  Usurf  Ubot  Qsurf  Qbot
        vals = line.split()
        if len(vals) >= 5:
            river_miles.append(float(vals[1]))   # col 2: river mile
            bed_els.append(float(vals[4]))        # col 5: bed elevation

river_miles = np.array(river_miles)
bed_els     = np.array(bed_els)
n_rows      = len(river_miles)

print(f"  Parsed {n_rows} spatial nodes.")
print(f"  River-mile range: {river_miles.min():.4f}  to  {river_miles.max():.4f}")

# =====================================================================
# 2.  User inputs
# =====================================================================
try:
    toe_input = input(
        "\nEnter the river-mile location of the Intrusion Toe "
        "(9.0-ppt isohaline): "
    )
    toe_rm = float(toe_input)
except (EOFError, ValueError) as exc:
    sys.exit(f"Invalid toe location: {exc}")

try:
    surf_input = input(
        "Enter the downstream distance (river miles) over which the\n"
        "  interface rises to the water surface [default 20]: "
    )
    surface_dist = float(surf_input) if surf_input.strip() else 20.0
except (EOFError, ValueError):
    surface_dist = 20.0

print(f"\n  Toe location           : RM {toe_rm:.4f}")
print(f"  Interface-surface dist : {surface_dist:.1f} RM downstream of toe")
print(f"  Surface WSEL           : {DEFAULT_WSEL:.4f} ft")

# =====================================================================
# 3.  Compute derived columns
# =====================================================================
wsel       = np.full(n_rows, DEFAULT_WSEL)
interface  = np.copy(bed_els)          # start at bed everywhere
surf_salt  = np.zeros(n_rows)
bot_salt   = np.zeros(n_rows)

# --- downstream-most river mile in the domain -------------------------
rm_min = river_miles.min()

# --- find the bed elevation at (or nearest to) the toe ----------------
idx_toe = int(np.argmin(np.abs(river_miles - toe_rm)))
bed_at_toe = bed_els[idx_toe]
interface_at_toe = bed_at_toe + 1.0    # 1 m above bed at the toe

# River mile where interface reaches the water surface (downstream)
rm_surface = toe_rm - surface_dist

for i in range(n_rows):
    rm = river_miles[i]

    # ── interface elevation ──────────────────────────────────────────
    if rm >= toe_rm:
        # Upstream of (or at) the toe → interface sits on the bed
        interface[i] = bed_els[i]
    elif rm >= rm_surface:
        # Between toe and surface point → linear interpolation
        frac = (toe_rm - rm) / surface_dist   # 0 at toe, 1 at surface pt
        interface[i] = interface_at_toe + frac * (DEFAULT_WSEL - interface_at_toe)
    else:
        # Downstream of surface point → interface at water surface
        interface[i] = DEFAULT_WSEL

    # ── salinity profiles (linear from toe to downstream boundary) ───
    if rm >= toe_rm:
        surf_salt[i] = 0.0
        bot_salt[i]  = 0.0
    else:
        dist_from_toe = toe_rm - rm
        total_dist    = toe_rm - rm_min
        if total_dist > 0:
            frac_salt = dist_from_toe / total_dist
        else:
            frac_salt = 0.0
        surf_salt[i] = 3.0  * frac_salt
        bot_salt[i]  = 35.0 * frac_salt

# =====================================================================
# 4.  Generate Salinity Profile Plot
# =====================================================================
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

print("\nGenerating salinity depth profile plot...")

# Invert river miles for plotting if they are not already descending
if river_miles[0] < river_miles[-1]:
    river_miles = river_miles[::-1]
    bed_els = bed_els[::-1]
    wsel = wsel[::-1]
    interface = interface[::-1]
    surf_salt = surf_salt[::-1]
    bot_salt = bot_salt[::-1]

# Define the vertical grid for the plot
y_res = 200  # Number of vertical points
min_elev = np.floor(bed_els.min()) - 5  # Add some buffer
max_elev = np.ceil(wsel.max()) + 5    # Add some buffer
y_grid = np.linspace(min_elev, max_elev, y_res)

# Create 2D arrays for plotting
x_grid_2d, y_grid_2d = np.meshgrid(river_miles, y_grid)
salinity_2d = np.full_like(x_grid_2d, np.nan) # Initialize with NaN

# Populate the salinity grid based on the two-layer system
for i in range(n_rows):
    rm_idx = i
    bed_elev = bed_els[rm_idx]
    interface_elev = interface[rm_idx]
    wsel_elev = wsel[rm_idx]
    
    # Assign salinity based on layer for each vertical point
    for j in range(y_res):
        y_pos = y_grid[j]
        if y_pos >= bed_elev and y_pos <= wsel_elev:
            if y_pos < interface_elev:
                salinity_2d[j, i] = bot_salt[rm_idx]  # Bottom layer
            else:
                salinity_2d[j, i] = surf_salt[rm_idx] # Top layer

# Create the plot
fig, ax = plt.subplots(figsize=(18, 6))

# Set up a colormap that emphasizes the 5-15 ppt range
colors = ["#1f77b4", "#2ca02c", "#ff7f0e", "#d62728"]  # Blue, Green, Orange, Red
nodes = [0.0, 5.0/35.0, 15.0/35.0, 1.0]
salinity_cmap = mcolors.LinearSegmentedColormap.from_list("salinity_cmap", list(zip(nodes, colors)))
norm = mcolors.Normalize(vmin=0, vmax=35)

# Plot the colored salinity field
im = ax.pcolormesh(x_grid_2d, y_grid_2d, salinity_2d, cmap=salinity_cmap, norm=norm, shading='gouraud')
cbar = fig.colorbar(im, ax=ax)
cbar.set_label('Salinity (ppt)')

# Plot contour lines for 0.25 and 9 ppt isohalines
ax.contour(x_grid_2d, y_grid_2d, salinity_2d, levels=[0.25], colors=['#8c564b'], linewidths=1.5, linestyles='--')

# Plot bed, water surface, and interface
ax.fill_between(river_miles, min_elev, bed_els, color='saddlebrown', alpha=0.8)
ax.plot(river_miles, wsel, color='cyan', linewidth=2, label='Surface Water Elevation')
ax.plot(river_miles, interface, color='white', linestyle=':', linewidth=2, label='Layer Interface')

# Find and plot the upstream extent of 0.25 ppt at the surface
surf_salt_mask = surf_salt >= 0.25
if np.any(surf_salt_mask):
    rm_surf_025_idx = np.where(surf_salt_mask)[0][-1]
    rm_surf_025 = river_miles[rm_surf_025_idx]
    wsel_at_point = wsel[rm_surf_025_idx]
    ax.plot(rm_surf_025, wsel_at_point, 'o', markerfacecolor='none', markeredgecolor='yellow', markersize=10, markeredgewidth=2, label='Upstream Extent (0.25 ppt surface)')

# The 9 ppt isohaline at the bed is the user-defined `toe_rm`
bed_at_toe_plot = bed_els[np.argmin(np.abs(river_miles - toe_rm))]
ax.plot(toe_rm, bed_at_toe_plot, 'o', markerfacecolor='none', markeredgecolor='lime', markersize=10, markeredgewidth=2, label='Upstream Extent (9 ppt bed)')

# Final plot setup
ax.set_xlabel("Distance from Head of Passes (river miles)", fontsize=12)
ax.set_ylabel("Elevation (ft)", fontsize=12)
ax.set_title("Longitudinal Salinity Depth Profile", fontsize=14)
ax.legend(loc='lower left')
ax.set_xlim(river_miles.max(), river_miles.min()) # Flipped axis
ax.set_ylim(min_elev, max_elev)
plt.grid(True, linestyle=':', alpha=0.4)
plt.tight_layout()

# Save the plot
plot_filename = "Salinity_Depth_Profile.png"
plt.savefig(plot_filename, dpi=300)
print(f"  Plot saved to {plot_filename}")

# =====================================================================
# 5.  Write the hot-start file
# =====================================================================
print(f"\nWriting synthetic hot-start file: {hotstart_file} ...")

with open(hotstart_file, "w") as out:
    # Single header line – Fortran skips exactly one line with Read(106,*)
    out.write("Time River-mile Wsel(ft) Layer-Interface-el(ft) Bed-el(ft) "
              "Surface-salt(ppt) Bottom-salt(ppt) "
              "Surface-vel(ft/sec) Bottom-vel(ft/sec) "
              "Surface-Q(cfs) Bottom-Q(cfs)\n")

    # Data rows: one line per row, 11 values
    # Format matches Fortran output: 9(f10.4,1x) , 2(f14.4,1x)
    for i in range(n_rows):
        out.write(
            f"{SYNTHETIC_TIME:10.4f} "
            f"{river_miles[i]:10.4f} "
            f"{wsel[i]:10.4f} "
            f"{interface[i]:10.4f} "
            f"{bed_els[i]:10.4f} "
            f"{surf_salt[i]:10.4f} "
            f"{bot_salt[i]:10.4f} "
            f"{0.0:10.4f} "
            f"{0.0:10.4f} "
            f"{0.0:14.4f} "
            f"{0.0:14.4f}\n"
        )

print(f"  {n_rows} rows written.  Done.")
