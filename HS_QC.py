#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
HS_QC.py — Hotstart Quality-Control Check
==========================================

Compares the most recently observed 9-ppt isohaline toe location (from a
spreadsheet) against the toe location embedded in the mr-hotstart.txt file
that initialised the current simulation.  If the two disagree by more than
a configurable river-mile threshold, the script signals the calling .bat
file to build a corrective synthetic hotstart.

Exit codes
----------
  0  No qualifying observation within the look-back window, **or** the
     observed and hotstart toe locations agree within the threshold.
     No further action is required.

  2  A recent observation exists **and** the mismatch exceeds the
     threshold.  The file ``hs_qc_result.txt`` is written to the
     simulation folder so that Create_RunFile_SyntheticHotstart.py
     can consume it automatically.

Expected spreadsheet layout  (MRSWAT_Data.xlsx → 'Observations')
-----------------------------------------------------------------
  Column A : observation date  (MM/DD/YYYY, YYYYMMDD, or any
             format pandas can parse — e.g. YYYY-MM-DD)
  Column B : 9-ppt isohaline river-mile location
"""

# ╔═════════════════════════════════════════════════════════════════════════╗
# ║  IMPORTS                                                              ║
# ╚═════════════════════════════════════════════════════════════════════════╝
import sys
from pathlib import Path
from datetime import datetime, timedelta

import numpy as np
import pandas as pd

# ╔═════════════════════════════════════════════════════════════════════════╗
# ║  USER CONFIGURATION  —  edit these values as needed                   ║
# ║                                                                       ║
# ║  These are the ONLY parameters a user should normally need to change. ║
# ╚═════════════════════════════════════════════════════════════════════════╝

# River-mile difference that triggers a corrective synthetic-hotstart run.
# If |observed_toe − hotstart_toe| > THRESHOLD_RM  →  exit code 2.
THRESHOLD_RM = 15.0                   # [river miles]

# How many calendar days backward from today to search for a qualifying
# observation.  Observations older than this are ignored.
LOOKBACK_DAYS = 3                     # [days]

# Salinity threshold (parts per thousand) used to identify the salt-wedge
# toe in both the observations sheet and the hotstart file.
SALT_THRESHOLD = 9.0                  # [ppt]

# Name of the Excel workbook (must reside in the workspace root).
XLSX_FILENAME = "MRSWAT_Data.xlsx"    # workbook name

# Name of the sheet inside the workbook that holds observation data.
SHEET_NAME = "Observations"           # sheet tab name

# Name of the hotstart file inside the simulation folder.
HOTSTART_FILENAME = "mr-hotstart.txt" # hotstart data file

# ╔═════════════════════════════════════════════════════════════════════════╗
# ║  DERIVED PATHS  —  computed automatically from script location        ║
# ╚═════════════════════════════════════════════════════════════════════════╝
# Folder layout assumed:
#   <workspace_root>/
#       MRSWAT_Data.xlsx
#       <SIMFOLDER>/
#           mr-hotstart.txt
#           Code Storage/          ← this script lives here
#               HS_QC.py

script_dir  = Path(__file__).resolve().parent   # …/<SIMFOLDER>/Code Storage/
master_dir  = script_dir.parent                 # …/<SIMFOLDER>/
basedir     = master_dir.parent                 # workspace root

xlsx_path     = basedir    / XLSX_FILENAME      # full path to workbook
hotstart_path = master_dir / HOTSTART_FILENAME  # full path to hotstart file

# ╔═════════════════════════════════════════════════════════════════════════╗
# ║  HELPER FUNCTIONS                                                     ║
# ╚═════════════════════════════════════════════════════════════════════════╝

def _parse_date(val):
    """Convert *val* to a ``pd.Timestamp``, handling multiple date formats.

    The function attempts the following conversions **in order**, returning
    the first that succeeds:

    1. **YYYYMMDD** — 8-digit integer *or* string (e.g. ``20260224``).
    2. **MM/DD/YYYY** — US-style slash-separated string (e.g. ``02/24/2026``).
       1- or 2-digit month/day values are accepted.
    3. **General fallback** — ``pd.to_datetime`` handles ISO 8601, Excel
       datetime objects, and many other representations.

    Returns ``pd.NaT`` when none of the above succeed.
    """
    # -- Guard: None / NaN --------------------------------------------------
    if val is None or (isinstance(val, float) and np.isnan(val)):
        return pd.NaT

    # -- 1. YYYYMMDD as integer ---------------------------------------------
    if isinstance(val, (int, np.integer)):
        s = str(int(val))
        if len(s) == 8:
            try:
                return pd.Timestamp(datetime.strptime(s, "%Y%m%d"))
            except ValueError:
                pass

    # -- 1b / 2. String-based formats ---------------------------------------
    elif isinstance(val, str):
        stripped = val.strip()

        # Pure 8-digit string → YYYYMMDD
        if stripped.isdigit() and len(stripped) == 8:
            try:
                return pd.Timestamp(datetime.strptime(stripped, "%Y%m%d"))
            except ValueError:
                pass

        # MM/DD/YYYY (1- or 2-digit month/day)
        try:
            return pd.Timestamp(datetime.strptime(stripped, "%m/%d/%Y"))
        except ValueError:
            pass

    # -- 3. General pandas fallback (YYYY-MM-DD, datetime objects, etc.) ----
    try:
        return pd.to_datetime(val)
    except Exception:
        return pd.NaT


def _print_banner():
    """Print a detailed description of inputs, assumptions, actions,
    and outputs to the terminal so the operator knows what the script
    will do before any work begins.
    """
    print("=" * 72)
    print("  HS_QC.py — Hotstart Quality-Control Check")
    print("=" * 72)
    print()
    print("  INPUTS")
    print("  ------")
    print(f"  1) Observations spreadsheet : {xlsx_path}")
    print(f"     Sheet '{SHEET_NAME}', Column A = date, Column B = 9-ppt")
    print( "     isohaline toe location (river mile).")
    print( "     Recognised date formats: MM/DD/YYYY, YYYYMMDD, YYYY-MM-DD,")
    print( "     or any format pandas can infer automatically.")
    print(f"  2) Hotstart file            : {hotstart_path}")
    print( "     Columns: col 1 = river mile, col 6 = bottom salinity (ppt).")
    print()
    print("  ASSUMPTIONS")
    print("  -----------")
    print(f"  • Qualifying observations must fall within the last "
          f"{LOOKBACK_DAYS} day(s)")
    print( "    relative to the current system date.")
    print(f"  • The 9-ppt isohaline toe in the hotstart is the highest")
    print(f"    river mile where bottom salinity >= {SALT_THRESHOLD} ppt.")
    print( "  • If no salinity exceeds the threshold the toe defaults to the")
    print( "    minimum (most downstream) river mile in the hotstart.")
    print(f"  • Agreement threshold: {THRESHOLD_RM} river miles.")
    print()
    print("  ACTIONS")
    print("  -------")
    print(f"  1. Read {XLSX_FILENAME} → '{SHEET_NAME}' sheet.")
    print( "  2. Parse observation dates (MM/DD/YYYY, YYYYMMDD, or general).")
    print(f"  3. Filter to observations within the last {LOOKBACK_DAYS} day(s).")
    print(f"  4. Parse {HOTSTART_FILENAME} to locate the 9-ppt toe (river mile).")
    print( "  5. Compare observed toe vs. hotstart toe.")
    print()
    print("  OUTPUTS")
    print("  -------")
    print( "  • Exit code 0 : No recent observation found, OR difference is")
    print(f"                   within {THRESHOLD_RM} RM.  No action needed.")
    print( "  • Exit code 2 : Difference exceeds threshold.  Writes")
    print( "                   hs_qc_result.txt with the observed toe RM")
    print( "                   so Create_RunFile_SyntheticHotstart.py can")
    print( "                   build a corrective run automatically.")
    print("=" * 72)
    print()
    print(f"  Workspace root  : {basedir}")
    print(f"  Simulation dir  : {master_dir}")
    print(f"  Script dir      : {script_dir}")
    print()


# ╔═════════════════════════════════════════════════════════════════════════╗
# ║  MAIN EXECUTION                                                       ║
# ╚═════════════════════════════════════════════════════════════════════════╝

# -- Print the description banner -------------------------------------------
_print_banner()

# ---------------------------------------------------------------------------
# Step 1 — Read the Excel observations
# ---------------------------------------------------------------------------
if not xlsx_path.exists():
    print(f"  ⚠ {XLSX_FILENAME} not found at {xlsx_path}.  Skipping QC.")
    sys.exit(0)

try:
    df = pd.read_excel(
        xlsx_path,
        sheet_name=SHEET_NAME,
        header=0,
        usecols=[0, 1],               # Column A = date, Column B = toe RM
    )
    df.columns = ["Date", "Toe_RM"]
    df = df.dropna(subset=["Date", "Toe_RM"])   # drop rows missing either value
    df["Date"] = df["Date"].apply(_parse_date)   # flexible date conversion
    df = df.dropna(subset=["Date"])              # drop rows whose date didn't parse
    print(f"  ✓ Parsed {len(df)} observation(s) with valid dates.")
except Exception as exc:
    print(f"  ⚠ Could not read '{SHEET_NAME}' sheet: {exc}.  Skipping QC.")
    sys.exit(0)

# ---------------------------------------------------------------------------
# Step 2 — Filter to the look-back window and pick the most recent row
# ---------------------------------------------------------------------------
today  = datetime.now()
cutoff = today - timedelta(days=LOOKBACK_DAYS)
recent = df[df["Date"] >= cutoff].copy()

if recent.empty:
    print(f"  ✓ No observations within the last {LOOKBACK_DAYS} days.  "
          f"No QC action needed.")
    sys.exit(0)

latest_row = recent.sort_values("Date", ascending=False).iloc[0]
obs_date   = latest_row["Date"]
obs_toe_rm = float(latest_row["Toe_RM"])
print(f"  ✓ Recent observation found: {obs_date.strftime('%Y-%m-%d')}  "
      f"Observed toe RM = {obs_toe_rm:.2f}")

# ---------------------------------------------------------------------------
# Step 3 — Locate the 9-ppt toe in the hotstart file
#           Toe = highest river mile where bottom salinity >= SALT_THRESHOLD
# ---------------------------------------------------------------------------
if not hotstart_path.exists():
    print(f"  ⚠ {HOTSTART_FILENAME} not found at {hotstart_path}.  Skipping QC.")
    sys.exit(0)

try:
    data = np.loadtxt(hotstart_path, skiprows=1)
except Exception as exc:
    print(f"  ⚠ Could not parse {HOTSTART_FILENAME}: {exc}.  Skipping QC.")
    sys.exit(0)

# Column layout (0-based): col 1 = river mile, col 6 = bottom salinity (ppt)
river_miles = data[:, 1]
bot_salt    = data[:, 6]

salt_mask = bot_salt >= SALT_THRESHOLD
if salt_mask.any():
    hs_toe_rm = float(river_miles[salt_mask].max())
else:
    # No salinity above threshold → place toe at the downstream end
    hs_toe_rm = float(river_miles.min())

print(f"  ✓ Hotstart 9-ppt toe location : RM {hs_toe_rm:.2f}")
print(f"  ✓ Observed 9-ppt toe location : RM {obs_toe_rm:.2f}")
diff = abs(obs_toe_rm - hs_toe_rm)
print(f"  ✓ Difference                  : {diff:.2f} RM  "
      f"(threshold = {THRESHOLD_RM} RM)")

# ---------------------------------------------------------------------------
# Step 4 — Evaluate and act
# ---------------------------------------------------------------------------
if diff <= THRESHOLD_RM:
    print(f"  ✓ Within threshold.  Hotstart is consistent with observations.")
    sys.exit(0)

# Difference exceeds threshold — write result file and signal the .bat script
print(f"\n  ⚠ Difference exceeds threshold ({diff:.2f} > {THRESHOLD_RM} RM).")
print(f"    A new synthetic hotstart using the observed toe "
      f"(RM {obs_toe_rm:.2f}) is required.")

result_path = master_dir / "hs_qc_result.txt"
with open(result_path, "w") as fh:
    fh.write(f"Observed_Toe_RM: {obs_toe_rm}\n")
    fh.write(f"Hotstart_Toe_RM: {hs_toe_rm}\n")
    fh.write(f"Difference_RM: {diff:.4f}\n")
    fh.write(f"Observation_Date: {obs_date.strftime('%Y-%m-%d')}\n")
    fh.write(f"Threshold_RM: {THRESHOLD_RM}\n")

print(f"  ✓ Wrote hs_qc_result.txt to: {result_path}")
sys.exit(2)
