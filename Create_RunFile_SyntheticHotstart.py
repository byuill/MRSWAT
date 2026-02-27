#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Create_RunFile_SyntheticHotstart.py
------------------------------------
Modifies the MRSWAT 'main.inp' run file in the master folder (one level above
the directory containing this script) to configure a hot-start simulation
using a **synthetic** initial condition rather than results from a previous
simulation.

The synthetic hotstart is constructed from:
  • The spatial grid (river miles) and bed elevations read from the first
    time-block of mr-full-output.txt in the Template folder.
  • A user-defined position of the upstream extent of the 9-ppt isohaline
    (the intrusion toe).
  • A transition zone (default 20 river miles downstream of the toe) over
    which the layer interface rises from the bed to near the water surface.
  • Downstream of the transition zone the interface stays at least 10 ft
    below the water surface elevation.

The script also sets main.inp for restart mode, identical to the real
hotstart script:
  • Total simulation time              → from Date_Information.txt
  • Frequency of Model Output          → 0.125 days
  • Restart conditions flag            → 1
  • Restart file name                  → 'mr-hotstart.txt'
  • Restart time                       → days between Ref_Date and Sim_Date

All other lines in main.inp are left unchanged.  Paths are resolved
dynamically from __file__ so the script works regardless of where the
project is located.
"""

from pathlib import Path
from datetime import datetime
import re
import shutil
import numpy as np

# ===========================================================================
# USER-CONFIGURABLE PARAMETERS
# ===========================================================================
TOE_RIVER_MILE       = 40.0     # River-mile of the 9-ppt isohaline (intrusion toe)
INTERFACE_RISE_DIST  = 20.0     # River miles downstream of toe for interface to
                                #   rise from bed to near the water surface
MIN_DEPTH_BELOW_WSEL = 10.0     # Minimum depth (ft) the interface must stay
                                #   below the water surface once it has risen
DEFAULT_WSEL         = 1.06     # Default surface-water elevation (ft NAVD88)
OCEAN_SALINITY       = 35.0     # ppt — downstream boundary salinity
SYNTHETIC_TIME       = 9.9999   # Flag value written in the Time column

# ---------------------------------------------------------------------------
# Paths — dynamic, relative to this script's location
# ---------------------------------------------------------------------------
script_dir  = Path(__file__).resolve().parent          # …/Code Storage/
master_dir  = script_dir.parent                        # …/MRSWAT_YYYYMMDD_XX/
inp_path    = master_dir / "main.inp"

# ---------------------------------------------------------------------------
# Override TOE_RIVER_MILE from hs_qc_result.txt if present
# (written by HS_QC.py when an observation mismatch triggers this script)
# ---------------------------------------------------------------------------
qc_result_path = master_dir / "hs_qc_result.txt"
if qc_result_path.exists():
    try:
        with open(qc_result_path, "r") as _fh:
            for _line in _fh:
                if _line.startswith("Observed_Toe_RM:"):
                    TOE_RIVER_MILE = float(_line.split(":", 1)[1].strip())
                    break
        print(f"  ✓ TOE_RIVER_MILE overridden by hs_qc_result.txt: "
              f"RM {TOE_RIVER_MILE}")
    except Exception as _e:
        print(f"  ⚠ Could not read hs_qc_result.txt ({_e}); "
              f"using default TOE_RIVER_MILE = {TOE_RIVER_MILE}")

# ---------------------------------------------------------------------------
# Read dates and elapsed days from Date_Information.txt
# ---------------------------------------------------------------------------
date_info_path = master_dir / "Date_Information.txt"
try:
    with open(date_info_path, "r") as _fh:
        _di_lines = [ln.rstrip("\n") for ln in _fh.readlines()]
    _headers      = _di_lines[0].split("\t")
    _dates_row    = _di_lines[1].split("\t")
    _elapsed_row  = _di_lines[2].split("\t")

    _ref_col      = _headers.index("Ref_Date")
    _final_col    = _headers.index("Final_Date")
    _sim_col      = _headers.index("Sim_Date")

    TOTAL_SIM_DAYS = int(_elapsed_row[_final_col])
    RESTART_DAYS   = int(_elapsed_row[_sim_col])

    CURRENT_REF_DATE = datetime.strptime(_dates_row[_ref_col], "%Y%m%d")
    CURRENT_SIM_DATE = datetime.strptime(_dates_row[_sim_col], "%Y%m%d")

    print(f"  ✓ Read Date_Information.txt:")
    print(f"      total simulation days = {TOTAL_SIM_DAYS}")
    print(f"      restart time (days)   = {RESTART_DAYS}")
except Exception as _e:
    TOTAL_SIM_DAYS = 210
    RESTART_DAYS   = 180
    CURRENT_REF_DATE = None
    CURRENT_SIM_DATE = None
    print(f"  ⚠ Could not read Date_Information.txt ({_e}); "
          f"using fallback values: total={TOTAL_SIM_DAYS}, restart={RESTART_DAYS}")

# ===========================================================================
# BUILD SYNTHETIC HOTSTART FILE
# ===========================================================================
print("\n  ── Creating synthetic hotstart file ──")
print(f"      Toe location (9-ppt isohaline) : RM {TOE_RIVER_MILE}")
print(f"      Interface rise distance         : {INTERFACE_RISE_DIST} RM")
print(f"      Min depth below WSEL            : {MIN_DEPTH_BELOW_WSEL} ft")
print(f"      Default WSEL                    : {DEFAULT_WSEL} ft")

# ── 1. Parse the first time-block from mr-full-output.txt ────────────────
#    We need the spatial grid (river miles) and bed elevations to build the
#    synthetic condition on the actual model geometry.
# -------------------------------------------------------------------------
full_output_path = master_dir / "mr-full-output.txt"
if not full_output_path.exists():
    raise FileNotFoundError(
        f"Full output file not found: {full_output_path}\n"
        f"The model must have been run at least once (even a cold start) "
        f"to produce this file, OR copy a template version."
    )

river_miles = []
bed_els     = []

print(f"  … Reading spatial grid from {full_output_path.name} …")
with open(full_output_path, "r") as fh:
    found_first_block = False
    for raw_line in fh:
        line = raw_line.strip()
        if not line:
            continue
        # Detect "Time = ..." header lines
        if line.startswith("Time"):
            if "=" in line:
                if found_first_block:
                    break                           # hit second block → done
                found_first_block = True
                continue
            else:
                continue                            # column header text
        # Skip non-numeric lines
        try:
            float(line.split()[0])
        except (ValueError, IndexError):
            continue
        vals = line.split()
        if len(vals) >= 5:
            river_miles.append(float(vals[1]))      # col 2: river mile
            bed_els.append(float(vals[4]))           # col 5: bed elevation (ft)

river_miles = np.array(river_miles)
bed_els     = np.array(bed_els)
n_rows      = len(river_miles)
print(f"  ✓ Parsed {n_rows} spatial nodes  "
      f"(RM {river_miles.min():.2f} to {river_miles.max():.2f})")

# ── 2. Compute synthetic fields ─────────────────────────────────────────
wsel       = np.full(n_rows, DEFAULT_WSEL)
interface  = np.copy(bed_els)                       # default: interface at bed
surf_salt  = np.zeros(n_rows)
bot_salt   = np.zeros(n_rows)

rm_min     = river_miles.min()                      # downstream boundary

# Bed elevation at (nearest to) the toe
idx_toe      = int(np.argmin(np.abs(river_miles - TOE_RIVER_MILE)))
bed_at_toe   = bed_els[idx_toe]
intf_at_toe  = bed_at_toe + 1.0                     # 1 ft above bed at the toe

# River-mile where the interface has fully risen (downstream of toe)
rm_surface   = TOE_RIVER_MILE - INTERFACE_RISE_DIST

# Target elevation once the interface has "surfaced":
# WSEL minus the minimum clearance (10 ft)
intf_near_surface = DEFAULT_WSEL - MIN_DEPTH_BELOW_WSEL

for i in range(n_rows):
    rm = river_miles[i]

    # ── Layer interface elevation ────────────────────────────────────
    if rm >= TOE_RIVER_MILE:
        # Upstream of (or at) the toe → interface sits on the bed
        interface[i] = bed_els[i]
    elif rm >= rm_surface:
        # Transition zone: linear interpolation from bed+1 ft at the toe
        # up to (WSEL − 10 m) at rm_surface
        frac = (TOE_RIVER_MILE - rm) / INTERFACE_RISE_DIST   # 0→1
        interface[i] = intf_at_toe + frac * (intf_near_surface - intf_at_toe)
    else:
        # Downstream of the transition → interface stays ≥ 10 ft below WSEL
        interface[i] = intf_near_surface

    # Enforce: interface can never drop below the bed
    if interface[i] < bed_els[i]:
        interface[i] = bed_els[i]

    # ── Salinity ─────────────────────────────────────────────────────
    if rm >= TOE_RIVER_MILE:
        # Fresh water upstream of the toe
        surf_salt[i] = 0.0
        bot_salt[i]  = 0.0
    else:
        # Linear ramp from 0 at the toe to ocean salinity at the
        # downstream boundary
        dist_from_toe = TOE_RIVER_MILE - rm
        total_dist    = TOE_RIVER_MILE - rm_min
        frac_salt     = dist_from_toe / total_dist if total_dist > 0 else 0.0
        # Surface layer: light salinity (capped at ~9 ppt at toe interface)
        surf_salt[i] = 3.0 * frac_salt
        # Bottom layer: full ocean gradient
        bot_salt[i]  = OCEAN_SALINITY * frac_salt

# ── 3. Write the synthetic hotstart file ─────────────────────────────────
hotstart_path = master_dir / "mr-hotstart.txt"

# Back up any existing file
if hotstart_path.exists():
    old_path = master_dir / "mr-hotstart_old.txt"
    shutil.move(str(hotstart_path), str(old_path))
    print(f"  ✓ Renamed existing mr-hotstart.txt → mr-hotstart_old.txt")

HOTSTART_HEADER = ("Time River-mile Wsel(ft) Layer-Interface-el(ft) "
                   "Bed-el(ft) Surface-salt(ppt) Bottom-salt(ppt) "
                   "Surface-vel(ft/sec) Bottom-vel(ft/sec) "
                   "Surface-Q(cfs) Bottom-Q(cfs)\n")

with open(hotstart_path, "w") as out:
    out.write(HOTSTART_HEADER)
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

print(f"  ✓ Wrote synthetic mr-hotstart.txt "
      f"(1 header + {n_rows} data lines)")

# ── 4. Write bookkeeping metadata ────────────────────────────────────────
metadata_path = master_dir / "hotstart_metadata.txt"
with open(metadata_path, "w") as fh:
    fh.write(f"Hotstart type: SYNTHETIC\n")
    fh.write(f"Toe river mile (9-ppt isohaline): {TOE_RIVER_MILE}\n")
    fh.write(f"Interface rise distance (RM): {INTERFACE_RISE_DIST}\n")
    fh.write(f"Min depth below WSEL (ft): {MIN_DEPTH_BELOW_WSEL}\n")
    fh.write(f"Default WSEL (ft): {DEFAULT_WSEL}\n")
    fh.write(f"Spatial nodes: {n_rows}\n")
print(f"  ✓ Wrote hotstart_metadata.txt")

# ===========================================================================
# MODIFY main.inp FOR HOT-START CONFIGURATION
# ===========================================================================

# ---------------------------------------------------------------------------
# Parameter labels to match (matched with str.strip() for robustness)
# ---------------------------------------------------------------------------
LABEL_SIM_DAYS      = "Total Time for Model Simulation (starting from time 0), DAYS"
LABEL_OUTPUT_FREQ   = "Frequency of Model Output, DAYS"
LABEL_RESTART_FLAG  = "Enter a 0 for no restart conditions, enter a 1 for restart conditions"
LABEL_RESTART_FILE  = "Enter the restart file name"
LABEL_RESTART_TIME  = "Enter the time to restart (DAYS)"

# ---------------------------------------------------------------------------
# Read the file
# ---------------------------------------------------------------------------
with open(inp_path, "r") as fh:
    lines = fh.readlines()

# ---------------------------------------------------------------------------
# Process lines — label / value pairs
# ---------------------------------------------------------------------------
output_lines = []
i = 0

while i < len(lines):
    raw = lines[i]
    label = raw.strip()

    # ── Total simulation time ────────────────────────────────────────
    if label == LABEL_SIM_DAYS:
        output_lines.append(raw)
        i += 1
        if i < len(lines):
            output_lines.append(f"{TOTAL_SIM_DAYS}\n")
            i += 1
        continue

    # ── Output frequency ─────────────────────────────────────────────
    if label == LABEL_OUTPUT_FREQ:
        output_lines.append(raw)
        i += 1
        if i < len(lines):
            output_lines.append("0.125\n")
            i += 1
        continue

    # ── Restart conditions flag ──────────────────────────────────────
    if label == LABEL_RESTART_FLAG:
        output_lines.append(raw)
        i += 1
        if i < len(lines):
            output_lines.append("1\n")
            i += 1

        # Insert restart file / time lines if missing
        next_label = lines[i].strip() if i < len(lines) else ""
        if next_label != LABEL_RESTART_FILE:
            output_lines.append(f"{LABEL_RESTART_FILE}\n")
            output_lines.append("mr-hotstart.txt\n")
            output_lines.append(f"{LABEL_RESTART_TIME}\n")
            output_lines.append(f"{RESTART_DAYS}\n")
        continue

    # ── Restart file name (already present) ──────────────────────────
    if label == LABEL_RESTART_FILE:
        output_lines.append(raw)
        i += 1
        if i < len(lines):
            output_lines.append("mr-hotstart.txt\n")
            i += 1
        continue

    # ── Restart time (already present) ───────────────────────────────
    if label == LABEL_RESTART_TIME:
        output_lines.append(raw)
        i += 1
        if i < len(lines):
            output_lines.append(f"{RESTART_DAYS}\n")
            i += 1
        continue

    # ── All other lines — pass through unchanged ─────────────────────
    output_lines.append(raw)
    i += 1

# ---------------------------------------------------------------------------
# Write the modified file back
# ---------------------------------------------------------------------------
with open(inp_path, "w") as fh:
    fh.writelines(output_lines)

print(f"\n  ✓ Hot-start run file written: {inp_path}")
print("    Changes applied:")
print(f"      Total simulation days → {TOTAL_SIM_DAYS}")
print(f"      Output frequency      → 0.125")
print(f"      Restart flag          → 1")
print(f"      Restart file          → mr-hotstart.txt")
print(f"      Restart time (days)   → {RESTART_DAYS}")
