#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Create_RunFile_Coldstart.py
----------------------------
Modifies the MRSWAT 'main.inp' run file in the master folder (one level above
the directory containing this script) to configure a cold-start simulation:

  • Total simulation time              → days between Ref_Date and Final_Date
                                         read from Date_Information.txt
  • Restart conditions flag            → 0  (no restart / cold start)
  • Restart file name label + value    → deleted (lines shifted up)
  • Restart time label + value         → deleted (lines shifted up)

All other lines are left unchanged.  Paths are resolved dynamically from
__file__ so the script works regardless of where the project is located.
"""

from pathlib import Path

# ---------------------------------------------------------------------------
# Paths — dynamic, relative to this script's location
# ---------------------------------------------------------------------------
script_dir  = Path(__file__).resolve().parent          # …/Code Storage/
master_dir  = script_dir.parent                        # …/MRSWAT_20260224_00/
inp_path    = master_dir / "main.inp"

# ---------------------------------------------------------------------------
# Read total simulation days from Date_Information.txt
# Layout: row 1 = headers (Ref_Date, Sim_Date, Final_Date)
#         row 2 = dates as YYYYMMDD
#         row 3 = elapsed days after reference date
# The Final_Date elapsed value in row 3 is the total simulation duration.
# ---------------------------------------------------------------------------
date_info_path = master_dir / "Date_Information.txt"
try:
    with open(date_info_path, "r") as _fh:
        _di_lines = [ln.rstrip("\n") for ln in _fh.readlines()]
    # row 3 (index 2) contains elapsed-day integers
    _elapsed_row  = _di_lines[2].split("\t")
    _headers      = _di_lines[0].split("\t")
    _final_col    = _headers.index("Final_Date")
    TOTAL_SIM_DAYS = int(_elapsed_row[_final_col])
    print(f"  ✓ Read Date_Information.txt: total simulation days = {TOTAL_SIM_DAYS}")
except Exception as _e:
    TOTAL_SIM_DAYS = 210   # safe fallback if file is missing or malformed
    print(f"  ⚠ Could not read Date_Information.txt ({_e}); "
          f"using fallback value of {TOTAL_SIM_DAYS} days")

# ---------------------------------------------------------------------------
# Parameter labels to match (matched with str.strip() for robustness)
# ---------------------------------------------------------------------------
LABEL_SIM_DAYS      = "Total Time for Model Simulation (starting from time 0), DAYS"
LABEL_RESTART_FLAG  = "Enter a 0 for no restart conditions, enter a 1 for restart conditions"
LABEL_RESTART_FILE  = "Enter the restart file name"
LABEL_RESTART_TIME  = "Enter the time to restart (DAYS)"
LABEL_OUTPUT_FREQ   = "Frequency of Model Output, DAYS"

# ---------------------------------------------------------------------------
# Read the file
# ---------------------------------------------------------------------------
with open(inp_path, "r") as fh:
    lines = fh.readlines()

# ---------------------------------------------------------------------------
# Process lines — the file alternates: label line, value line, label line, …
# We walk index-by-index so we can look ahead for the value line.
# ---------------------------------------------------------------------------
output_lines = []
i = 0
while i < len(lines):
    raw = lines[i]
    label = raw.strip()

    # ── Total simulation time ────────────────────────────────────────────
    if label == LABEL_SIM_DAYS:
        output_lines.append(raw)                      # keep the label unchanged
        i += 1
        if i < len(lines):
            output_lines.append(f"{TOTAL_SIM_DAYS}\n")  # value from Date_Information.txt
            i += 1
        continue

    # ── Output frequency ─────────────────────────────────────────────────
    if label == LABEL_OUTPUT_FREQ:
        output_lines.append(raw)
        i += 1
        if i < len(lines):
            output_lines.append("1.0\n")  # set output frequency to 1.0 days
            i += 1
        continue

    # ── Restart conditions flag ──────────────────────────────────────────
    if label == LABEL_RESTART_FLAG:
        output_lines.append(raw)
        i += 1
        if i < len(lines):
            output_lines.append("0\n")    # 0 = cold start
            i += 1
        continue

    # ── Restart file name — delete label AND value lines ─────────────────
    if label == LABEL_RESTART_FILE:
        i += 1          # skip label line
        if i < len(lines):
            i += 1      # skip value line
        continue

    # ── Restart time — delete label AND value lines ──────────────────────
    if label == LABEL_RESTART_TIME:
        i += 1          # skip label line
        if i < len(lines):
            i += 1      # skip value line
        continue

    # ── All other lines — pass through unchanged ─────────────────────────
    output_lines.append(raw)
    i += 1

# ---------------------------------------------------------------------------
# Write the modified file back
# ---------------------------------------------------------------------------
with open(inp_path, "w") as fh:
    fh.writelines(output_lines)

print(f"  ✓ Cold-start run file written: {inp_path}")
print("    Changes applied:")
print(f"      {LABEL_SIM_DAYS!r:>8} → {TOTAL_SIM_DAYS}")
print(f"      {LABEL_OUTPUT_FREQ!r:>8} → 1.0")
print(f"      {LABEL_RESTART_FLAG!r:>8} → 0")
print(f"      {LABEL_RESTART_FILE!r:>8} → deleted (label + value)")
print(f"      {LABEL_RESTART_TIME!r:>8} → deleted (label + value)")
