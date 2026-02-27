#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Create_RunFile_Hotstart.py
----------------------------
Modifies the MRSWAT 'main.inp' run file in the master folder (one level above
the directory containing this script) to configure a hot-start simulation:

  • Total simulation time              → days between Ref_Date and Final_Date
                                         read from Date_Information.txt
  • Frequency of Model Output          → 0.125 days
  • Restart conditions flag            → 1  (restart / hot start)
  • Restart file name                  → 'mr-hotstart.txt'
  • Restart time                       → days between Ref_Date and Sim_Date

If the restart file name and restart time entries are missing from main.inp
(as is the case for cold-start versions), they are inserted immediately after
the restart conditions flag value line.

Optionally (CREATE_HOTSTART_FILE = True), the script also:
  • Identifies the immediately-previous simulation folder in the same parent
    directory, using the naming convention MRSWAT_YYYYMMDD_## to determine
    chronological order.
  • Reads mr-full-output.txt from that previous simulation and extracts the
    results matrix whose simulated time is closest to the start of the current
    simulation day (Sim_Date).
  • Renames any pre-existing mr-hotstart.txt → mr-hotstart_old.txt.
  • Writes the extracted matrix as mr-hotstart.txt (no headings).
  • Writes hotstart_metadata.txt with the actual calendar date of the hotstart.

All other lines in main.inp are left unchanged.  Paths are resolved
dynamically from __file__ so the script works regardless of where the project
is located.
"""

from pathlib import Path
from datetime import datetime, timedelta
import re
import shutil

# ===========================================================================
# USER SWITCHES
# ===========================================================================
CREATE_HOTSTART_FILE = True   # Set True to build a new mr-hotstart.txt from
                              # the previous simulation's full output

# ---------------------------------------------------------------------------
# Paths — dynamic, relative to this script's location
# ---------------------------------------------------------------------------
script_dir  = Path(__file__).resolve().parent          # …/Code Storage/
master_dir  = script_dir.parent                        # …/MRSWAT_20260224_XX/
inp_path    = master_dir / "main.inp"

# ---------------------------------------------------------------------------
# Read dates and elapsed days from Date_Information.txt
# Layout: row 1 = headers (Ref_Date, Sim_Date, Final_Date)
#         row 2 = dates as YYYYMMDD
#         row 3 = elapsed days after reference date
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

    TOTAL_SIM_DAYS = int(_elapsed_row[_final_col])   # Ref_Date → Final_Date
    RESTART_DAYS   = int(_elapsed_row[_sim_col])      # Ref_Date → Sim_Date

    # Also store actual dates for hotstart creation
    CURRENT_REF_DATE = datetime.strptime(_dates_row[_ref_col], "%Y%m%d")
    CURRENT_SIM_DATE = datetime.strptime(_dates_row[_sim_col], "%Y%m%d")

    print(f"  ✓ Read Date_Information.txt:")
    print(f"      total simulation days = {TOTAL_SIM_DAYS}")
    print(f"      restart time (days)   = {RESTART_DAYS}")
except Exception as _e:
    TOTAL_SIM_DAYS = 210   # safe fallback if file is missing or malformed
    RESTART_DAYS   = 180
    CURRENT_REF_DATE = None
    CURRENT_SIM_DATE = None
    print(f"  ⚠ Could not read Date_Information.txt ({_e}); "
          f"using fallback values: total={TOTAL_SIM_DAYS}, restart={RESTART_DAYS}")

# ===========================================================================
# HOTSTART FILE CREATION  (runs only when the switch is True)
# ===========================================================================
if CREATE_HOTSTART_FILE:
    print("\n  ── Creating hotstart file ──")

    # ── 1. Parse current folder name ─────────────────────────────────────
    current_folder_name = master_dir.name           # e.g. MRSWAT_20260224_02
    folder_re = re.compile(r"^MRSWAT_(\d{8})_(\d{2})$")
    m_cur = folder_re.match(current_folder_name)
    if not m_cur:
        raise ValueError(
            f"Current folder '{current_folder_name}' does not match "
            f"the expected MRSWAT_YYYYMMDD_## naming convention."
        )
    cur_date_str = m_cur.group(1)
    cur_seq      = int(m_cur.group(2))

    # ── 2. Inventory all simulation folders in the same parent directory ──
    parent_dir   = master_dir.parent
    sim_folders  = []                             # (date_str, seq, Path)
    for d in parent_dir.iterdir():
        if d.is_dir():
            m = folder_re.match(d.name)
            if m:
                sim_folders.append((m.group(1), int(m.group(2)), d))
    sim_folders.sort(key=lambda x: (x[0], x[1]))  # chronological order

    # ── 3. Find the folder immediately before the current one ────────────
    current_idx = None
    for idx, (ds, seq, _) in enumerate(sim_folders):
        if ds == cur_date_str and seq == cur_seq:
            current_idx = idx
            break
    if current_idx is None or current_idx == 0:
        raise FileNotFoundError(
            "Cannot identify a previous simulation folder. "
            "The current folder appears to be the earliest one."
        )
    prev_folder = sim_folders[current_idx - 1][2]
    print(f"  ✓ Previous simulation folder: {prev_folder.name}")

    # ── 4. Read the previous simulation's reference date ─────────────────
    prev_date_info_path = prev_folder / "Date_Information.txt"
    with open(prev_date_info_path, "r") as fh:
        prev_di_lines = [ln.rstrip("\n") for ln in fh.readlines()]
    prev_headers   = prev_di_lines[0].split("\t")
    prev_dates_row = prev_di_lines[1].split("\t")
    prev_ref_col   = prev_headers.index("Ref_Date")
    prev_ref_date  = datetime.strptime(prev_dates_row[prev_ref_col], "%Y%m%d")
    print(f"  ✓ Previous simulation reference date: "
          f"{prev_ref_date.strftime('%Y-%m-%d')}")

    # ── 5. Compute the target elapsed days ───────────────────────────────
    #    "closest in time to the start of the current day" →
    #    Sim_Date of this simulation expressed as elapsed days from the
    #    *previous* simulation's reference date.
    target_elapsed = (CURRENT_SIM_DATE - prev_ref_date).days
    print(f"  ✓ Target elapsed days (in prev sim reference): {target_elapsed}")

    # ── 6. Scan the previous simulation's full output for Time = markers ─
    full_output_path = prev_folder / "mr-full-output.txt"
    if not full_output_path.exists():
        raise FileNotFoundError(
            f"Full output file not found: {full_output_path}"
        )

    time_re = re.compile(r"^Time\s*=\s*([\d.]+)")
    time_entries = []                             # (elapsed_days, line_number)
    print(f"  … Scanning {full_output_path.name} for time markers …")
    with open(full_output_path, "r") as fh:
        for line_num, line in enumerate(fh):
            tm = time_re.match(line.strip())
            if tm:
                time_entries.append((float(tm.group(1)), line_num))

    if not time_entries:
        raise ValueError("No 'Time = …' markers found in the full output file.")

    # ── 7. Find the time step closest to the target ──────────────────────
    best_idx = min(
        range(len(time_entries)),
        key=lambda k: abs(time_entries[k][0] - target_elapsed),
    )
    best_time, best_line = time_entries[best_idx]
    hotstart_actual_date = prev_ref_date + timedelta(days=best_time)
    print(f"  ✓ Closest time step: {best_time} elapsed days  "
          f"(actual date: {hotstart_actual_date.strftime('%Y-%m-%d')})")

    # ── 8. Determine the line range for the data block ───────────────────
    data_start = best_line + 1                    # first data line after Time =
    if best_idx + 1 < len(time_entries):
        data_end = time_entries[best_idx + 1][1]  # exclusive (next Time = line)
    else:
        data_end = None                           # read to end of file

    # ── 9. Second pass — extract the data lines ─────────────────────────
    data_lines = []
    with open(full_output_path, "r") as fh:
        for line_num, line in enumerate(fh):
            if line_num < data_start:
                continue
            if data_end is not None and line_num >= data_end:
                break
            data_lines.append(line)
    print(f"  ✓ Extracted {len(data_lines)} data lines from results matrix")

    # ── 10. Rename any existing hotstart file in the current folder ──────
    hotstart_path = master_dir / "mr-hotstart.txt"
    if hotstart_path.exists():
        old_path = master_dir / "mr-hotstart_old.txt"
        shutil.move(str(hotstart_path), str(old_path))
        print(f"  ✓ Renamed existing mr-hotstart.txt → mr-hotstart_old.txt")

    # ── 11. Write the new hotstart file (with header expected by Fortran) ──
    HOTSTART_HEADER = ("Time River-mile Wsel(ft) Layer-Interface-el(ft) "
                       "Bed-el(ft) Surface-salt(ppt) Bottom-salt(ppt) "
                       "Surface-vel(ft/sec) Bottom-vel(ft/sec) "
                       "Surface-Q(cfs) Bottom-Q(cfs)\n")
    with open(hotstart_path, "w") as fh:
        fh.write(HOTSTART_HEADER)
        fh.writelines(data_lines)
    print(f"  ✓ Wrote new mr-hotstart.txt (1 header + {len(data_lines)} data lines)")

    # ── 12. Write bookkeeping metadata ───────────────────────────────────
    metadata_path = master_dir / "hotstart_metadata.txt"
    with open(metadata_path, "w") as fh:
        fh.write(f"Hotstart actual date: "
                 f"{hotstart_actual_date.strftime('%Y-%m-%d')}\n")
        fh.write(f"Elapsed days in previous simulation: {best_time}\n")
        fh.write(f"Source file: {full_output_path}\n")
        fh.write(f"Previous simulation folder: {prev_folder.name}\n")
        fh.write(f"Previous simulation reference date: "
                 f"{prev_ref_date.strftime('%Y-%m-%d')}\n")
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
# Process lines — the file alternates: label line, value line, label line, …
# We walk index-by-index so we can look ahead for the value line.
# ---------------------------------------------------------------------------
output_lines = []
found_restart_file = False
found_restart_time = False
i = 0

while i < len(lines):
    raw = lines[i]
    label = raw.strip()

    # ── Total simulation time ────────────────────────────────────────────
    if label == LABEL_SIM_DAYS:
        output_lines.append(raw)                        # keep the label unchanged
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
            output_lines.append("0.125\n")               # hot-start output frequency
            i += 1
        continue

    # ── Restart conditions flag ──────────────────────────────────────────
    if label == LABEL_RESTART_FLAG:
        output_lines.append(raw)
        i += 1
        if i < len(lines):
            output_lines.append("1\n")                    # 1 = hot start
            i += 1

        # Check whether the next lines already contain the restart file /
        # restart time entries.  If not, insert them here.
        next_label = lines[i].strip() if i < len(lines) else ""
        if next_label != LABEL_RESTART_FILE:
            # Lines are missing — insert the four restart lines
            output_lines.append(f"{LABEL_RESTART_FILE}\n")
            output_lines.append("mr-hotstart.txt\n")
            output_lines.append(f"{LABEL_RESTART_TIME}\n")
            output_lines.append(f"{RESTART_DAYS}\n")
        continue

    # ── Restart file name (already present) ──────────────────────────────
    if label == LABEL_RESTART_FILE:
        found_restart_file = True
        output_lines.append(raw)
        i += 1
        if i < len(lines):
            output_lines.append("mr-hotstart.txt\n")     # ensure correct filename
            i += 1
        continue

    # ── Restart time (already present) ───────────────────────────────────
    if label == LABEL_RESTART_TIME:
        found_restart_time = True
        output_lines.append(raw)
        i += 1
        if i < len(lines):
            output_lines.append(f"{RESTART_DAYS}\n")      # value from Date_Information.txt
            i += 1
        continue

    # ── All other lines — pass through unchanged ─────────────────────────
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
