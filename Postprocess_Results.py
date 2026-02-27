#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
=============================================================================
 Visualize_Results.py  –  MRSWAT Post-Processing Visualization Suite
=============================================================================
 Purpose
 -------
   Loads MRSWAT model output files from the master simulation folder (one
   level above the directory containing this script) and produces a
   configurable set of plots and animations:

     1. **Summary Figure** – 8.5 × 11 in portrait page with:
            • Time-series of salt-wedge intrusion distance, TMDL threshold
              line, upstream discharge, and shaded forecast region.
            • Scatter plot of discharge vs intrusion distance coloured by
              downstream stage relative to its temporal mean.
            • 10-day and 28-day forecast summary text.
     2. **Satellite Salinity Map** – Esri World Imagery basemap overlaid
            with bottom-layer salinity markers at integer river miles from
            RM −19 (SW Pass outlet) to RM 100.
     3. **Salinity Animation** – MP4 (or GIF fallback) showing the time-
            evolving longitudinal salinity cross-section with boundary-
            condition time-series and a moving time marker.

 Configuration
 -------------
   Toggle each output on/off with the boolean switches in SECTION 1.

 Data Sources  (read from the master folder)
 ------------
   • mr-toe-location-output.txt  – toe-location time series
   • mr-full-output.txt          – full longitudinal output (all time steps)
   • Discharges.txt              – upstream discharge boundary condition
   • Stages.txt                  – downstream tidal stage boundary condition
   • Date_Information.txt        – reference date and simulation date
   • river_miles.csv             – USACE river-mile marker lat/lon lookup
                                   (expected in the same folder as this script)

 Outputs  (saved to <master_folder>/Output/)
 -------
   • MRSWAT_<START>_<END>.jpg     – summary figure  (if MAKE_SUMMARY_FIGURE)
   • Map.jpg                      – satellite salinity map (if MAKE_SALINITY_MAP)
   • MRSWAT_<START>_<END>.mp4/.gif – animation (if MAKE_ANIMATION)

 Physical System Overview
 -------------------------
   MRSWAT (Mississippi River Salt Water Intrusion model) is a 1-D, two-layer,
   cross-sectionally averaged hydrodynamic and salinity transport model for
   the tidal freshwater and estuarine reach of the Mississippi River.

   Two-layer estuarine structure
   ------------------------------
   The lower Mississippi operates as a partially-mixed salt-wedge estuary.
   When river discharge is sufficiently low, denser saline Gulf water intrudes
   upstream along the bottom of the channel beneath the lighter, seaward-
   flowing freshwater layer above.  This creates two distinct layers separated
   by a sharp halocline (the *interface*):

     ┌─────────────────────────────────────┐  ← Water Surface Elevation (WSEL)
     │   UPPER LAYER  –  fresh to brackish │    (driven by river stage and tides)
     │   (s_surf ≈ 0–5 ppt)               │
     ├─────────────────────────────────────┤  ← Interface Elevation
     │   LOWER LAYER  –  saline            │    (salt wedge / halocline)
     │   (s_bot ≈ 5–35 ppt)               │
     └─────────────────────────────────────┘  ← Channel Bed Elevation

   Salt wedge dynamics
   --------------------
     • High river discharge → strong seaward barotropic pressure gradient →
       salt wedge is pushed downstream toward Head of Passes (RM 0) or offshore.
     • Low river discharge → reduced flushing → saline Gulf water migrates
       upstream, threatening intakes for municipal water supply and industrial use.
     • Tidal forcing modulates the interface on a 12.4-hour (semi-diurnal) cycle,
       causing the wedge to surge upstream on the flood tide and retreat on ebb.

   Colormap – salinity classification (USGS / Venice system)
   -----------------------------------------------------------
     Blue   (0 ppt)   – Fresh water
     Green  (~5 ppt)  – Oligohaline  (1–5 ppt boundary)
     Orange (~15 ppt) – Mesohaline / Polyhaline
     Red    (35 ppt)  – Marine / Euhaline

   Output file columns (mr-full-output.txt, 11 columns per row)
   --------------------------------------------------------------
     Index  Variable
     ─────  ─────────────────────────────────────────────────────
       0    Node index
       1    River Mile (RM) from Head of Passes
       2    Water Surface Elevation – WSEL (ft, NAVD 88)
       3    Interface Elevation (ft, NAVD 88)
       4    Bed Elevation (ft, NAVD 88)
       5    Surface-layer salinity – s_surf (ppt)
       6    Bottom-layer salinity  – s_bot  (ppt)
       7-10 Additional model state variables (not used here)
=============================================================================
"""

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 0 — IMPORTS
# ═══════════════════════════════════════════════════════════════════════════════
import os
import sys
import csv as _csv
import datetime as dt
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import matplotlib
# Prefer non-interactive backend when running in batch; allow override
if os.environ.get("MPLBACKEND") is None:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.dates as mdates
import matplotlib.colors as mcolors
import matplotlib.patheffects
import matplotlib.gridspec as gridspec
import matplotlib.animation as animation


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 1 — USER-CONFIGURABLE SWITCHES
# ═══════════════════════════════════════════════════════════════════════════════
# Set each flag to True/False to enable/disable individual outputs.

MAKE_SUMMARY_FIGURE    = False  # Time-series + scatter summary (JPG)
MAKE_SALINITY_MAP      = False  # Satellite basemap salinity map (JPG)
MAKE_MAP_ANIMATION     = False  # Animated satellite salinity map (MP4/GIF)
MAKE_ANIMATION         = False  # Longitudinal salinity animation (MP4/GIF)
MAKE_CALIBRATION_ANIMATION = True  # Longitudinal salinity + observed CTD overlay (MP4/GIF)
MAKE_COMBINED_ANIMATION= False   # Combined 4-quadrant animation (MP4/GIF)
RECORD_TO_EXCEL        = False   # Record isohaline forecast to MRSWAT_Data.xlsx

# ── Animation parameters ─────────────────────────────────────────────────────
MAX_RIVER_MILE = 90.0          # Animation spatial domain (RM from Head of Passes)
FPS            = 10            # Frames per second in the output video

# ── Alternative toe definition ────────────────────────────────────────────────
#   When ALT_TOE_ENABLED is True the code computes an alternative salt-wedge
#   toe from the full model output and shows it alongside the main toe (which
#   is always read from mr-toe-location-output.txt).
#
#   The alternative toe is still defined as the 9 ppt isohaline, but uses a
#   connectivity rule based on the full model output rather than the simple
#   upstream-most location reported by the Fortran model.
#
#   A bottom salinity layer is considered *continuous* (connected) when every
#   node satisfies both:
#     • bottom salinity  ≥  ALT_TOE_SALINITY  (ppt)  ← continuity threshold,
#                                                        NOT the toe isohaline
#     • bottom-layer thickness  >  ALT_TOE_MIN_THICKNESS  (ft)
#   Gaps shorter than ALT_TOE_MAX_GAP (RM) do not break the layer.
#   The reported alternative toe is then the most upstream node within that
#   connected layer where bottom salinity ≥ 9 ppt.
#
#   Salinity pools are saline segments that are *not* part of the connected
#   layer.  A pool is reported only if its max bottom-layer thickness
#   exceeds ALT_TOE_MIN_POOL_DEPTH (ft).

ALT_TOE_ENABLED       = False  # Set True to show alternative toe on figures
ALT_TOE_SALINITY      = 5.0   # ppt – min salinity for layer CONTINUITY (not the toe isohaline – that is always 9 ppt)
ALT_TOE_MIN_THICKNESS = 1.0   # ft  – min bottom-layer thickness for connectivity
ALT_TOE_MAX_GAP       = 2.0   # RM  – max gap before the layer is broken
ALT_TOE_MIN_POOL_DEPTH = 5.0  # ft  – only report pools thicker than this


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 2 — DYNAMIC PATHS
# ═══════════════════════════════════════════════════════════════════════════════
SCRIPT_DIR = Path(__file__).resolve().parent          # …/Code Storage/
MASTER_DIR = SCRIPT_DIR.parent                        # …/dummy_sim/ (or dated run)
OUTPUT_DIR = MASTER_DIR / "Output"
OUTPUT_DIR.mkdir(exist_ok=True)

# Input files in the master folder
TOE_FILE       = MASTER_DIR / "mr-toe-location-output.txt"
FULL_OUT_FILE  = MASTER_DIR / "mr-full-output.txt"
DISCHARGE_FILE = MASTER_DIR / "Discharges.txt"
STAGE_FILE     = MASTER_DIR / "Stages.txt"
DATE_FILE      = MASTER_DIR / "Date_Information.txt"

# River-mile CSV lives alongside this script
RIVER_MILES_CSV = SCRIPT_DIR / "river_miles.csv"

# Network path to vessel-based CTD field observations
FIELD_DATA_ROOT = Path(
    r"\\mvd\mvn\Data_EDR\02 Stream Gaging Section"
    r"\Gaging Stations\Gages_RandomProjects\Saltwater Wedge\Field Data"
)


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 3 — SHARED UTILITY FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════════════

def load_reference_and_sim_dates():
    """
    Read reference and simulation dates from Date_Information.txt.

    File format (tab-separated, two columns):
        Ref_Date    Sim_Date
        YYYYMMDD    YYYYMMDD   ← row 2

    Returns
    -------
    ref_date : datetime.datetime  – day-0 of model elapsed time
    sim_date : datetime.date      – simulation start date
    """
    date_info = np.loadtxt(str(DATE_FILE), dtype=str)
    if date_info.shape[0] < 2 or date_info.shape[1] < 2:
        raise SystemExit(f"{DATE_FILE} must have at least two rows and two columns.")
    ref_str = date_info[1, 0].strip()
    sim_str = date_info[1, 1].strip()
    ref_date = dt.datetime.strptime(ref_str, "%Y%m%d")
    sim_date = dt.datetime.strptime(sim_str, "%Y%m%d").date()
    return ref_date, sim_date


def load_two_col_ts(filepath):
    """
    Read a two-column boundary-condition time series from a space/tab-
    delimited text file with a single header row.

    Used for:
      Discharges.txt – upstream open boundary condition.
          Column 0: elapsed time (days since reference date)
          Column 1: river discharge (cfs).
          Discharge is the dominant control on salt-wedge position.

      Stages.txt – downstream open boundary condition.
          Column 0: elapsed time (days since reference date)
          Column 1: water surface elevation / tidal stage (ft, NAVD 88).
          The ~0.5–1.5 ft semi-diurnal tidal range modulates the salt wedge.

    Returns
    -------
    times  : np.ndarray  – elapsed-day timestamps
    values : np.ndarray  – corresponding values
    """
    times, values = [], []
    with open(filepath, "r") as f:
        lines = f.readlines()
    for line in lines[1:]:
        parts = line.strip().split()
        if len(parts) >= 2:
            try:
                times.append(float(parts[0]))
                values.append(float(parts[1]))
            except ValueError:
                continue
    return np.array(times), np.array(values)


def load_toe_data():
    """
    Load toe-location output (skip header row).

    Columns:
        0 – elapsed days since reference
        1 – River Mile of the toe (Toe-Location)
        2 – River Mile of the 0.25-ppt-at-surface TMDL threshold
    """
    data = np.loadtxt(str(TOE_FILE), skiprows=1)
    return data[:, 0], data[:, 1], data[:, 2]


def compute_connected_toe(data, sal_threshold=None, min_thickness=None,
                          max_gap=None, min_pool_depth=None):
    """
    Identify an alternative salt-wedge toe as the upstream extent of a
    **continuous** significant salinity layer connected to the downstream
    boundary, and classify disconnected saline segments as pools.

    Connectivity criterion (both must hold at every node):
      • bottom-layer thickness  > *min_thickness*  (ft)
      • bottom-layer salinity   ≥ *sal_threshold*  (ppt)  ← continuity threshold

    Note: *sal_threshold* controls which nodes count as "connected" to the
    downstream boundary.  The reported toe is always the upstream-most node
    within that connected layer where bottom salinity ≥ 9 ppt — identical
    isohaline to the main toe from ``mr-toe-location-output.txt``.
    Gaps shorter than *max_gap* RM are bridged (treated as continuous).

    Salinity pools are qualifying segments that are not absorbed into the
    connected layer.  Only pools whose **max** bottom-layer thickness
    exceeds *min_pool_depth* are returned.

    Parameters
    ----------
    data : np.ndarray, shape (N, 11)
        Single time-step snapshot from ``read_mrswat_output()``.
    sal_threshold : float, optional
        Minimum bottom salinity for the continuous layer (default:
        ``ALT_TOE_SALINITY``).
    min_thickness : float, optional
        Minimum bottom-layer thickness (default: ``ALT_TOE_MIN_THICKNESS``).
    max_gap : float, optional
        Maximum gap that is bridged (default: ``ALT_TOE_MAX_GAP``).
    min_pool_depth : float, optional
        Report only pools thicker than this (default:
        ``ALT_TOE_MIN_POOL_DEPTH``).

    Returns
    -------
    toe_rm : float or np.nan
        River mile of the alternative connected-layer 9-ppt toe.
    pools : list of dict
        Each pool dict contains:
          'lo'  : float – downstream RM of the pool
          'hi'  : float – upstream RM of the pool
          'avg_thickness' : float – mean bottom-layer thickness (ft)
          'max_thickness' : float – max bottom-layer thickness (ft)
    """
    if sal_threshold is None:
        sal_threshold = ALT_TOE_SALINITY
    if min_thickness is None:
        min_thickness = ALT_TOE_MIN_THICKNESS
    if max_gap is None:
        max_gap = ALT_TOE_MAX_GAP
    if min_pool_depth is None:
        min_pool_depth = ALT_TOE_MIN_POOL_DEPTH

    if len(data) == 0:
        return np.nan, []

    idx       = np.argsort(data[:, 1])          # sort downstream → upstream
    rm        = data[idx, 1]
    bot_sal   = data[idx, 6]
    intf_elev = data[idx, 3]
    bed_elev  = data[idx, 4]
    bot_depth = intf_elev - bed_elev            # bottom-layer thickness (ft)

    # ── Build contiguous saline segments ──────────────────────────────────
    segments = []       # list of lists of node indices
    current_seg = []
    for i in range(len(rm)):
        if bot_depth[i] > min_thickness and bot_sal[i] >= sal_threshold:
            current_seg.append(i)
        else:
            if current_seg:
                segments.append(current_seg[:])
                current_seg = []
    if current_seg:
        segments.append(current_seg[:])

    if not segments:
        return np.nan, []

    # ── Merge segments across small gaps ─────────────────────────────────
    #    Segment 0 seeds the connected layer.  Subsequent segments whose
    #    gap to the connected layer is < max_gap are absorbed.
    connected_nodes = list(segments[0])
    connected_hi    = rm[segments[0][-1]]
    pools_raw = []

    for seg in segments[1:]:
        seg_lo = rm[seg[0]]
        gap    = seg_lo - connected_hi

        if gap < max_gap:
            connected_nodes.extend(seg)
            connected_hi = rm[seg[-1]]
        else:
            pools_raw.append(seg)

    # ── Identify alternative toe within the connected layer ──────────────
    toe_rm = np.nan
    for i in connected_nodes:
        if bot_sal[i] >= 9.0:
            toe_rm = rm[i]              # last (most upstream) wins

    # ── Build pool dicts, filtering by depth ─────────────────────────────
    pools = []
    for seg in pools_raw:
        thicknesses = [bot_depth[j] for j in seg]
        max_th = float(np.max(thicknesses))
        if max_th > min_pool_depth:
            pools.append({
                'lo':  float(rm[seg[0]]),
                'hi':  float(rm[seg[-1]]),
                'avg_thickness': float(np.mean(thicknesses)),
                'max_thickness': max_th,
            })

    return toe_rm, pools


def _format_pool_text(pools):
    """
    Return a human-readable string summarising salinity-pool locations and
    thicknesses.

    Parameters
    ----------
    pools : list of dict
        Each dict has keys 'lo', 'hi', 'avg_thickness', 'max_thickness'.

    Returns
    -------
    str – e.g. ``"RM 65.2–67.8 (avg 12.3 ft, max 18.1 ft), RM 72.1 (8.7 ft)"``
          or ``"none"`` when the list is empty.
    """
    if not pools:
        return "none"
    parts = []
    for p in pools:
        lo, hi = p['lo'], p['hi']
        if abs(hi - lo) < 0.1:
            loc = f"RM {lo:.1f}"
        else:
            loc = f"RM {lo:.1f}\u2013{hi:.1f}"
        parts.append(f"{loc} (avg {p['avg_thickness']:.1f} ft, "
                     f"max {p['max_thickness']:.1f} ft)")
    return ", ".join(parts)


def compute_alt_toe_timeseries():
    """
    Compute the alternative toe time series and per-snapshot pool lists
    from the full model output using ``compute_connected_toe``.

    Returns
    -------
    snapshots   : list of dict         – raw model snapshots (for reuse)
    times       : np.ndarray           – elapsed days since reference date
    alt_locs    : np.ndarray           – alternative toe RM per snapshot
    all_pools   : list of list of dict – pool dicts for each snapshot
    """
    snapshots = read_mrswat_output()
    times, alt_locs, all_pools = [], [], []
    for snap in snapshots:
        alt_rm, pools = compute_connected_toe(snap["data"])
        times.append(snap["time"])
        alt_locs.append(alt_rm)
        all_pools.append(pools)
    return (snapshots, np.array(times), np.array(alt_locs), all_pools)


def load_stage_data():
    """
    Read downstream stage, tolerating incomplete rows.
    """
    rows = []
    with open(str(STAGE_FILE), "r") as f:
        for line in f:
            parts = line.split()
            if len(parts) == 2:
                try:
                    rows.append([float(parts[0]), float(parts[1])])
                except ValueError:
                    pass
    arr = np.array(rows)
    return arr[:, 0], arr[:, 1]


def read_mrswat_output():
    """
    Parse the MRSWAT full longitudinal output file into a list of time
    snapshots.

    File structure
    ──────────────
    The file is organised as a sequence of time blocks.  Each block begins
    with a header ``Time = <fractional_days>`` followed by one row per
    computational node (11 whitespace-delimited floats):

      Col  0  – Node index
      Col  1  – River Mile (RM)
      Col  2  – WSEL (ft, NAVD 88)
      Col  3  – Interface Elevation (ft, NAVD 88)
      Col  4  – Bed Elevation (ft, NAVD 88)
      Col  5  – Surface-layer salinity (ppt)
      Col  6  – Bottom-layer salinity  (ppt)
      Col 7-10 – Additional variables (not used)

    Returns
    -------
    list of dict  –  [{'time': float, 'data': np.ndarray (N, 11)}, …]
    """
    filename = str(FULL_OUT_FILE)
    if not os.path.exists(filename):
        print(f"Error: {filename} not found.")
        sys.exit(1)
    print(f"Reading {filename} …")

    snapshots    = []
    current_time = None
    current_data = []

    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("Time"):
                if current_time is not None and current_data:
                    snapshots.append({"time": current_time,
                                      "data": np.array(current_data)})
                try:
                    parts = line.split("=")
                    if len(parts) > 1:
                        current_time = float(parts[1])
                    current_data = []
                except ValueError:
                    pass
                continue
            try:
                parts = line.split()
                float(parts[0])
                if len(parts) == 11:
                    current_data.append([float(x) for x in parts])
            except (ValueError, IndexError):
                continue

    if current_time is not None and current_data:
        snapshots.append({"time": current_time,
                          "data": np.array(current_data)})
    print(f"Parsed {len(snapshots)} time steps.")
    return snapshots


def interpolate_profile(rm, bed, wsel, interface, s_surf, s_bot, y_grid):
    """
    Map the MRSWAT two-layer salinity solution onto a regular 2-D Cartesian
    grid (river-mile × elevation) for pseudo-colour rendering.

    MRSWAT is a *layer-averaged* model: each node stores one salinity for
    the upper layer (s_surf) and one for the lower layer (s_bot).  We paint
    every grid cell within the water column with the salinity of whichever
    layer it belongs to:

        elevation > interface  →  upper layer  →  s_surf
        elevation ≤ interface  →  lower layer  →  s_bot

    Cells above the water surface or below the bed are NaN (transparent).
    'gouraud' shading then gives smooth horizontal interpolation.

    Returns
    -------
    X, Y : (M, N) arrays  – river-mile and elevation meshes
    Z    : (M, N) array   – salinity (ppt); NaN outside water column
    """
    X, Y = np.meshgrid(rm, y_grid)
    Z    = np.full_like(X, np.nan)
    for i in range(len(rm)):
        elevs    = Y[:, i]
        in_water = (elevs >= bed[i]) & (elevs <= wsel[i])
        Z[in_water & (elevs <= interface[i]), i] = s_bot[i]
        Z[in_water & (elevs >  interface[i]), i] = s_surf[i]
    return X, Y, Z


def build_salinity_cmap():
    """
    Custom salinity colormap anchored to standard estuarine zones
    (Venice / USGS classification, 0–35 ppt range):
      0 ppt → blue (fresh)  |  5 ppt → green (oligohaline)
      15 ppt → orange        |  35 ppt → red  (marine)
    """
    colors = ["#1f77b4", "#2ca02c", "#ff7f0e", "#d62728"]
    nodes  = [0.0, 5.0 / 35.0, 15.0 / 35.0, 1.0]
    cmap = mcolors.LinearSegmentedColormap.from_list(
        "salinity_cmap", list(zip(nodes, colors)))
    norm = mcolors.Normalize(vmin=0, vmax=35)
    return cmap, norm


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 4 — RIVER-MILE GEOMETRY  (shared by map and animation)
# ═══════════════════════════════════════════════════════════════════════════════

def load_river_mile_waypoints():
    """
    Load USACE river-mile marker coordinates from the companion CSV.

    Two USACE series cover the model domain (RM −19 … RM 100):
      • MISSISSIPPI-LO  (RIVER_CODE "MI") — whole miles RM 0 … 100
      • SOUTHWEST PASS   (RIVER_CODE "SW") — SW-Pass miles 0 … 20
        (SW-Pass mile N ↔ model RM −N)

    Returns
    -------
    np.ndarray, shape (K, 3) – columns: model_RM, latitude, longitude
    """
    mi_rows, sw_rows = [], []
    with open(str(RIVER_MILES_CSV), "r") as cf:
        reader = _csv.DictReader(cf)
        for row in reader:
            code = row["RIVER_CODE"].strip()
            name = row["RIVER_NAME"].strip()
            try:
                mile = float(row["MILE"])
                lat  = float(row["LATITUDE1"])
                lon  = float(row["LONGITUDE1"])
            except (ValueError, KeyError):
                continue
            if code == "MI" and name == "MISSISSIPPI-LO" and 0 <= mile <= 100:
                mi_rows.append((mile, lat, lon))
            elif code == "SW" and 0 < mile <= 20:
                sw_rows.append((-mile, lat, lon))
    all_rows = mi_rows + sw_rows
    all_rows.sort(key=lambda r: r[0])
    return np.array(all_rows)


def rm_to_latlon(rm_values, waypoints):
    """Linearly interpolate (lat, lon) for arbitrary river-mile values."""
    lats = np.interp(rm_values, waypoints[:, 0], waypoints[:, 1])
    lons = np.interp(rm_values, waypoints[:, 0], waypoints[:, 2])
    return lats, lons


def latlon_to_web_mercator(lat, lon):
    """
    WGS-84 (lat, lon) → Web Mercator (EPSG:3857) without pyproj.
    """
    x = lon * 20037508.34 / 180.0
    y = np.log(np.tan((90.0 + lat) * np.pi / 360.0)) * 20037508.34 / np.pi
    return x, y


# ── Legacy hardcoded waypoints (kept for reference only) ─────────────
# Originally embedded in plot_model_results.py; see _OLD_RM_WAYPOINTS in
# that file.  Not used here — the CSV table above is more maintainable.


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 5 — FFMPEG RESOLVER  (for animation output)
# ═══════════════════════════════════════════════════════════════════════════════

def _resolve_ffmpeg():
    """
    Return the path to a usable ffmpeg executable, or None.

    Resolution order (first match wins):
      1. System PATH (covers system-wide installs and conda environments)
      2. WinGet store (Gyan.FFmpeg via ``winget install Gyan.FFmpeg``)
      3. imageio-ffmpeg bundled binary
    """
    import shutil as _shutil
    import glob as _glob

    exe = _shutil.which("ffmpeg")
    if exe:
        return exe

    winget_pattern = os.path.join(
        os.environ.get("LOCALAPPDATA", ""),
        "Microsoft", "WinGet", "Packages",
        "Gyan.FFmpeg_*", "**", "bin", "ffmpeg.exe",
    )
    matches = _glob.glob(winget_pattern, recursive=True)
    if matches:
        return matches[0]

    try:
        import imageio_ffmpeg as iioff
        return iioff.get_ffmpeg_exe()
    except Exception:
        pass
    return None


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 6 — SUMMARY FIGURE
# ═══════════════════════════════════════════════════════════════════════════════

def make_summary_figure(ref_date, sim_date):
    """
    Produce the 8.5 × 11 in summary figure:
      Panel 1 – intrusion distance + TMDL + discharge time-series
      Panel 2 – discharge vs intrusion distance scatter
    Saves to <OUTPUT_DIR>/MRSWAT_<START>_<END>.jpg
    """
    print("\n" + "=" * 70)
    print("  Generating Summary Figure …")
    print("=" * 70)

    # ── Load data ────────────────────────────────────────────────────────
    toe_time, toe_loc, tmdl_loc = load_toe_data()
    flow_data = np.loadtxt(str(DISCHARGE_FILE), skiprows=1)
    flow_time, flow_q = flow_data[:, 0], flow_data[:, 1]
    stage_time, stage_wsel = load_stage_data()
    stage_mean = np.mean(stage_wsel)

    # Alt toe (optional)
    alt_toe_loc = None
    all_pools   = None
    if ALT_TOE_ENABLED:
        _, alt_times, alt_locs, all_pools = compute_alt_toe_timeseries()
        # Interpolate onto the same time vector as the main toe
        alt_toe_loc = np.interp(toe_time, alt_times, alt_locs)

    today = dt.datetime.now().date()
    if sim_date != today:
        offset_days = (today - sim_date).days
        print(f"  WARNING: simulation is {abs(offset_days)} day(s) from the "
              f"current date (offset = {offset_days}).")

    # Convert elapsed days → calendar dates
    toe_dates  = np.array([ref_date + dt.timedelta(days=float(d)) for d in toe_time])
    flow_dates = np.array([ref_date + dt.timedelta(days=float(d)) for d in flow_time])

    now = dt.datetime.now()
    forecast_mask  = toe_dates > now
    forecast_start = toe_dates[np.argmax(forecast_mask)] if np.any(forecast_mask) else None

    # ── Figure layout ────────────────────────────────────────────────────
    fig = plt.figure(figsize=(8.5, 11))
    page_w, page_h = 8.5, 11.0
    margin = 0.25
    plot_h = 2.0
    left   = 1.00 / page_w
    right  = 0.90 / page_w
    gap    = 0.60
    ax_w   = 1.0 - left - right
    ax_h   = plot_h / page_h

    top1 = 1.0 - margin / page_h - ax_h
    ax1  = fig.add_axes([left, top1, ax_w, ax_h])

    top2 = top1 - (gap + 0.5) / page_h - ax_h
    ax2  = fig.add_axes([left, top2, ax_w, ax_h])

    # ── PLOT 1 – Time series ─────────────────────────────────────────────
    color_toe  = "#1f77b4"
    color_tmdl = "#d62728"
    color_flow = "#2ca02c"
    color_alt  = "#ff7f0e"          # orange for alternative toe

    ln1 = ax1.plot(toe_dates, toe_loc,  color=color_toe,  lw=1.2,
                   label="intrusion distance (9 ppt)")
    ln2 = ax1.plot(toe_dates, tmdl_loc, color=color_tmdl, lw=1.2,
                   label="TMDL threshold (0.25 ppt)")
    lns_extra = []
    if alt_toe_loc is not None:
        _ln_alt = ax1.plot(toe_dates, alt_toe_loc, color=color_alt, lw=1.2,
                           ls=":", label="alt toe (connected)")
        lns_extra = _ln_alt
    ax1.set_xlabel("Date", fontsize=10)
    ax1.set_ylabel("River Mile", fontsize=10, color="black")
    ax1.tick_params(axis="y", labelcolor="black")

    ax1r = ax1.twinx()
    ln3  = ax1r.plot(flow_dates, flow_q / 1000.0, color=color_flow,
                     lw=1.0, ls="--", alpha=0.8, label="Discharge")
    ax1r.set_ylabel("Discharge (×1 000 cfs)", fontsize=10, color=color_flow)
    ax1r.tick_params(axis="y", labelcolor=color_flow)

    lns  = ln1 + ln2 + lns_extra + ln3
    labs = [l.get_label() for l in lns]
    ax1.legend(lns, labs, loc="upper center", bbox_to_anchor=(0.5, -0.18),
               ncol=2, fontsize=8, framealpha=0.9)
    ax1.set_title("Saltwater Intrusion – Time Series", fontsize=12,
                  fontweight="bold")
    ax1.grid(True, lw=0.3, alpha=0.6)

    ax1.xaxis_date()
    ax1.xaxis.set_major_locator(mdates.AutoDateLocator())
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b-%d"))
    fig.autofmt_xdate(rotation=30)

    # Shade forecast region
    from matplotlib.patches import Patch
    if forecast_start is not None:
        ax1.axvspan(mdates.date2num(forecast_start),
                    mdates.date2num(toe_dates[-1]),
                    color="lightblue", alpha=0.40, zorder=0)
        forecast_patch = Patch(facecolor="lightblue", alpha=0.60, label="Forecast")
        ax1.legend(lns + [forecast_patch],
                   [l.get_label() for l in lns] + ["Forecast"],
                   loc="upper center", bbox_to_anchor=(0.5, -0.18),
                   ncol=2, fontsize=8, framealpha=0.9)

    # Sill line at RM 64
    ax1.axhline(64.0, color="k", ls="--", lw=1.0)

    # ── PLOT 2 – Discharge vs Intrusion Distance (scatter) ───────────────
    q_matched     = np.interp(toe_time, flow_time, flow_q)
    stage_matched = np.interp(toe_time, stage_time, stage_wsel)
    colors = np.where(stage_matched > stage_mean, "#FF0000", "#0000FF")
    base_s = 18

    ax2.scatter(q_matched / 1000.0, toe_loc, c=colors, s=base_s,
                alpha=0.8, edgecolors="none")

    if len(toe_dates) > 0:
        diffs   = np.array([abs((d - now).total_seconds()) for d in toe_dates])
        idx_now = int(np.argmin(diffs))
        idx_last = len(toe_dates) - 1
        hl_s = base_s * 2
        ax2.scatter([q_matched[idx_now] / 1000.0], [toe_loc[idx_now]],
                    facecolors="none", edgecolors="k", s=hl_s, lw=1.2,
                    zorder=5, label="current date")
        ax2.annotate("Now",
                     xy=(q_matched[idx_now] / 1000.0, toe_loc[idx_now]),
                     xytext=(4, 4), textcoords="offset points",
                     fontsize=8, color="k", fontweight="bold", zorder=6)
        ax2.scatter([q_matched[idx_last] / 1000.0], [toe_loc[idx_last]],
                    facecolors="none", edgecolors="b", s=hl_s, lw=1.2,
                    zorder=5, label="last date")
        ax2.annotate("Last",
                     xy=(q_matched[idx_last] / 1000.0, toe_loc[idx_last]),
                     xytext=(4, 4), textcoords="offset points",
                     fontsize=8, color="b", fontweight="bold", zorder=6)

    ax2.set_xlabel("Discharge (×1 000 cfs)", fontsize=10)
    ax2.set_ylabel("Intrusion Distance (River Mile)", fontsize=10)
    ax2.set_title("Discharge vs Intrusion Distance", fontsize=12,
                  fontweight="bold")
    ax2.axhline(64.0, color="k", ls="--", lw=1.0)
    ax2.grid(True, lw=0.3, alpha=0.6)

    from matplotlib.lines import Line2D
    legend_el = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#FF0000",
               markersize=8, label=f"Stage > mean ({stage_mean:+.2f} ft)"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#0000FF",
               markersize=8, label=f"Stage ≤ mean ({stage_mean:+.2f} ft)"),
    ]
    ax2.legend(handles=legend_el, loc="upper left",
               bbox_to_anchor=(1.02, 1.0), borderaxespad=0,
               fontsize=8, framealpha=0.9)

    # ── 10-day and 28-day forecast values ────────────────────────────────
    def _fc_val(days_ahead):
        target = now + dt.timedelta(days=days_ahead)
        if target < toe_dates[0] or target > toe_dates[-1]:
            return None, None
        td = (target - ref_date).total_seconds() / 86400.0
        return np.interp(td, toe_time, toe_loc), np.interp(td, toe_time, tmdl_loc)

    toe_10, tmdl_10 = _fc_val(10)
    toe_28, tmdl_28 = _fc_val(28)

    l10 = (f"  10-day forecast  – Toe: {toe_10:.2f} RM  |  0.25-ppt: {tmdl_10:.2f} RM"
           if toe_10 is not None
           else "  10-day forecast  – not available in the current result time series.")
    l28 = (f"  28-day forecast  – Toe: {toe_28:.2f} RM  |  0.25-ppt: {tmdl_28:.2f} RM"
           if toe_28 is not None
           else "  28-day forecast  – not available in the current result time series.")

    # Alt toe + pool text (only when alt toe is enabled)
    alt_lines = ""
    if ALT_TOE_ENABLED and alt_toe_loc is not None and all_pools is not None:
        def _alt_at_day(days_ahead):
            target = now + dt.timedelta(days=days_ahead)
            td = (target - ref_date).total_seconds() / 86400.0
            if td < alt_times[0] or td > alt_times[-1]:
                return None
            return float(np.interp(td, alt_times, alt_locs))

        def _pools_at_day(days_ahead):
            target = now + dt.timedelta(days=days_ahead)
            td = (target - ref_date).total_seconds() / 86400.0
            if td < alt_times[0] or td > alt_times[-1]:
                return []
            idx = int(np.argmin(np.abs(alt_times - td)))
            return all_pools[idx] if idx < len(all_pools) else []

        alt_10 = _alt_at_day(10)
        alt_28 = _alt_at_day(28)
        pools_10 = _pools_at_day(10)
        pools_28 = _pools_at_day(28)
        print(f"  Pools at forecast day 10 : {_format_pool_text(pools_10)}")
        print(f"  Pools at forecast day 28 : {_format_pool_text(pools_28)}")
        a10 = f"  10-day alt toe   – {alt_10:.2f} RM" if alt_10 is not None else ""
        a28 = f"  28-day alt toe   – {alt_28:.2f} RM" if alt_28 is not None else ""
        p10 = f"    Pools (day 10): {_format_pool_text(pools_10)}"
        p28 = f"    Pools (day 28): {_format_pool_text(pools_28)}"
        alt_lines = "\n".join(s for s in [a10, p10, a28, p28] if s)
        if alt_lines:
            alt_lines = "\n" + alt_lines

    fc_text = "Forecast summary:\n" + l10 + "\n" + l28 + alt_lines
    print(fc_text)

    fig.text(0.5, top2 - 0.55 / page_h, fc_text,
             ha="center", va="top", fontsize=9, family="monospace",
             bbox=dict(boxstyle="round,pad=0.4", facecolor="#f0f0f0",
                       edgecolor="gray", lw=0.8))

    # ── Save ──────────────────────────────────────────────────────────────
    fig_start = (ref_date + dt.timedelta(days=float(toe_time[0]))).strftime("%Y%m%d")
    fig_end   = (ref_date + dt.timedelta(days=float(toe_time[-1]))).strftime("%Y%m%d")
    fig_name  = OUTPUT_DIR / f"MRSWAT_{fig_start}_{fig_end}.jpg"
    plt.savefig(str(fig_name), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ Summary figure saved: {fig_name}")


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 7 — SATELLITE SALINITY MAP
# ═══════════════════════════════════════════════════════════════════════════════

def make_salinity_map(ref_date):
    """
    Satellite basemap with bottom-layer salinity markers (RM −19 to RM 100)
    for the snapshot closest to the current system date.
    Saves to <OUTPUT_DIR>/Map.jpg
    """
    print("\n" + "=" * 70)
    print("  Generating Satellite Salinity Map …")
    print("=" * 70)

    try:
        import contextily as ctx
    except ImportError:
        print("  WARNING: contextily not installed – skipping salinity map.")
        return

    now        = dt.datetime.now()
    target_day = (now - ref_date).total_seconds() / 86400.0

    # Stream-parse for the single closest snapshot
    best_time, best_data = None, None
    snap_time, snap_data = None, []
    best_diff = 1e30
    passed_target = False

    print(f"  Scanning {FULL_OUT_FILE.name} for snapshot closest to day {target_day:.3f} …")
    with open(str(FULL_OUT_FILE), "r") as fo:
        for line in fo:
            ls = line.strip()
            if not ls:
                continue
            if ls.startswith("Time"):
                if snap_time is not None and snap_data:
                    diff = abs(snap_time - target_day)
                    if diff < best_diff:
                        best_diff = diff
                        best_time = snap_time
                        best_data = snap_data[:]
                    if snap_time >= target_day:
                        passed_target = True
                    if passed_target and diff > best_diff:
                        break
                try:
                    snap_time = float(ls.split("=")[1])
                except (IndexError, ValueError):
                    snap_time = None
                snap_data = []
                continue
            try:
                pts = ls.split()
                float(pts[0])
                if len(pts) == 11:
                    snap_data.append([float(p) for p in pts])
            except (ValueError, IndexError):
                continue

    if snap_time is not None and snap_data:
        diff = abs(snap_time - target_day)
        if diff < best_diff:
            best_time = snap_time
            best_data = snap_data

    snap_arr = np.array(best_data)
    snap_date_str = (ref_date + dt.timedelta(days=best_time)).strftime("%b %d %Y")
    print(f"  Using snapshot at time = {best_time:.4f} days ({snap_date_str})")

    MAP_RM_MIN, MAP_RM_MAX = -19.0, 100.0
    mask    = (snap_arr[:, 1] >= MAP_RM_MIN) & (snap_arr[:, 1] <= MAP_RM_MAX)
    map_rm    = snap_arr[mask, 1]
    map_sal   = snap_arr[mask, 6]
    map_surf  = snap_arr[mask, 5]

    marker_rms = np.arange(int(np.ceil(MAP_RM_MIN)), int(np.floor(MAP_RM_MAX)) + 1)
    marker_sal  = np.interp(marker_rms, map_rm, map_sal)
    marker_surf = np.interp(marker_rms, map_rm, map_surf)

    waypoints = load_river_mile_waypoints()
    marker_lat, marker_lon = rm_to_latlon(marker_rms.astype(float), waypoints)
    marker_x, marker_y     = latlon_to_web_mercator(marker_lat, marker_lon)

    cl_lat, cl_lon = rm_to_latlon(map_rm, waypoints)
    cl_x, cl_y     = latlon_to_web_mercator(cl_lat, cl_lon)

    cmap, norm = build_salinity_cmap()

    fig_map, ax_map = plt.subplots(figsize=(10, 14))
    ax_map.plot(cl_x, cl_y, color="white", lw=1.0, alpha=0.40, zorder=1)

    sc = ax_map.scatter(marker_x, marker_y, c=marker_sal,
                        cmap=cmap, norm=norm, s=30,
                        edgecolors="k", linewidths=0.3, zorder=3)
    cbar = fig_map.colorbar(sc, ax=ax_map,
                            label="Bottom-Layer Salinity (ppt)",
                            shrink=0.55, pad=0.02)
    cbar.ax.yaxis.label.set_color("white")
    cbar.ax.tick_params(colors="white")

    for i, rm_val in enumerate(marker_rms):
        if rm_val % 10 == 0:
            ax_map.annotate(
                f"RM {rm_val:g}",
                xy=(marker_x[i], marker_y[i]),
                xytext=(6, 4), textcoords="offset points",
                fontsize=7, fontweight="bold", color="white",
                path_effects=[matplotlib.patheffects.withStroke(
                    linewidth=2, foreground="black")],
                zorder=4)

    # ── Infrastructure labels at key river miles ──────────────────────────
    INFRA_LABELS = {
        19: "Boothville",
        49: "Port Sulphur/PALH",
        76: "Belle Chasse",
        81: "Dalcour",
        88: "St. Bernard",
        96: "Algiers",
    }
    for rm_val, label_text in INFRA_LABELS.items():
        idx_arr = np.where(marker_rms == rm_val)[0]
        if len(idx_arr) == 0:
            continue
        idx = idx_arr[0]
        ax_map.annotate(
            label_text,
            xy=(marker_x[idx], marker_y[idx]),
            xytext=(8, -10), textcoords="offset points",
            fontsize=7.5, fontweight="bold", color="yellow",
            path_effects=[matplotlib.patheffects.withStroke(
                linewidth=2, foreground="black")],
            zorder=6)

    # ── Salinity threshold rings ──────────────────────────────────────────
    from matplotlib.patches import Circle

    # ── Main toe from output file ────────────────────────────────────────
    toe_time_f, toe_loc_f, _ = load_toe_data()
    _main_toe_rm = float(np.interp(best_time, toe_time_f, toe_loc_f))

    # White ring (solid): main 9-ppt toe from mr-toe-location-output.txt
    _idx_9 = None
    if not np.isnan(_main_toe_rm):
        _idx_9 = int(np.argmin(np.abs(marker_rms - _main_toe_rm)))
        circ_white = Circle((marker_x[_idx_9], marker_y[_idx_9]),
                            radius=(marker_x.max() - marker_x.min()) * 0.006,
                            transform=ax_map.transData,
                            facecolor='none', edgecolor='white', linewidth=1.5,
                            zorder=6)
        ax_map.add_patch(circ_white)
        ax_map.annotate(
            "9 ppt (toe)",
            xy=(marker_x[_idx_9], marker_y[_idx_9]),
            xytext=(8, 4), textcoords="offset points",
            fontsize=7.5, fontweight="bold", color="white",
            path_effects=[matplotlib.patheffects.withStroke(
                linewidth=2, foreground="black")],
            zorder=7)

    # ── Alternative toe + salinity pools (conditional) ────────────────────
    _pools = []
    if ALT_TOE_ENABLED:
        _alt_toe_rm, _pools = compute_connected_toe(snap_arr)
        # Orange ring (dashed): alternative connected toe
        if not np.isnan(_alt_toe_rm):
            _idx_alt = int(np.argmin(np.abs(marker_rms - _alt_toe_rm)))
            circ_alt = Circle((marker_x[_idx_alt], marker_y[_idx_alt]),
                              radius=(marker_x.max() - marker_x.min()) * 0.006,
                              transform=ax_map.transData,
                              facecolor='none', edgecolor='orange', linewidth=1.5,
                              linestyle='--', zorder=6)
            ax_map.add_patch(circ_alt)
            ax_map.annotate(
                "alt toe",
                xy=(marker_x[_idx_alt], marker_y[_idx_alt]),
                xytext=(8, -12), textcoords="offset points",
                fontsize=7, fontweight="bold", color="orange",
                path_effects=[matplotlib.patheffects.withStroke(
                    linewidth=2, foreground="black")],
                zorder=7)

        # Yellow dashed rings: salinity pools (disconnected from main toe)
        for _pool in _pools:
            _p_center = (_pool['lo'] + _pool['hi']) / 2.0
            _p_idx = int(np.argmin(np.abs(marker_rms - _p_center)))
            circ_pool = Circle((marker_x[_p_idx], marker_y[_p_idx]),
                               radius=(marker_x.max() - marker_x.min()) * 0.006,
                               transform=ax_map.transData,
                               facecolor='none', edgecolor='yellow', linewidth=1.5,
                               linestyle='--', zorder=6)
            ax_map.add_patch(circ_pool)
            if abs(_pool['hi'] - _pool['lo']) < 0.1:
                _pool_label = f"Pool RM {_pool['lo']:.0f}"
            else:
                _pool_label = f"Pool RM {_pool['lo']:.0f}\u2013{_pool['hi']:.0f}"
            ax_map.annotate(
                _pool_label,
                xy=(marker_x[_p_idx], marker_y[_p_idx]),
                xytext=(8, -12), textcoords="offset points",
                fontsize=7, fontweight="bold", color="yellow",
                path_effects=[matplotlib.patheffects.withStroke(
                    linewidth=2, foreground="black")],
                zorder=7)

        print(f"  Salinity pools: {_format_pool_text(_pools)}")

    # Red ring (dashed): upstream-most marker where SURFACE salinity > 0.25 ppt
    _mask_025 = marker_surf > 0.25
    if np.any(_mask_025):
        _idx_025 = int(np.max(np.where(_mask_025)[0]))
        circ_red = Circle((marker_x[_idx_025], marker_y[_idx_025]),
                          radius=(marker_x.max() - marker_x.min()) * 0.006,
                          transform=ax_map.transData,
                          facecolor='none', edgecolor='red', linewidth=2.0,
                          linestyle='--', zorder=7)
        ax_map.add_patch(circ_red)
        # place red label to the left of the marker
        ax_map.annotate(
            f"0.25 ppt",
            xy=(marker_x[_idx_025], marker_y[_idx_025]),
            xytext=(-65, 6), textcoords="offset points",
            fontsize=7.5, fontweight="bold", color="red",
            path_effects=[matplotlib.patheffects.withStroke(
                linewidth=2, foreground="black")],
            zorder=8)

    pad_x = (marker_x.max() - marker_x.min()) * 0.08
    pad_y = (marker_y.max() - marker_y.min()) * 0.05
    ax_map.set_xlim(marker_x.min() - pad_x, marker_x.max() + pad_x)
    ax_map.set_ylim(marker_y.min() - pad_y, marker_y.max() + pad_y)

    try:
        ctx.add_basemap(ax_map, source=ctx.providers.Esri.WorldImagery,
                        zoom="auto", attribution=False)
        print("  Satellite basemap tiles loaded.")
    except Exception as e:
        print(f"  WARNING: Could not fetch satellite tiles ({e}).")

    ax_map.set_title(
        f"Mississippi River Salinity  ·  {snap_date_str}\n"
        f"RM {int(np.ceil(MAP_RM_MIN))} (SW Pass outlet) to RM {int(MAP_RM_MAX)}",
        fontsize=13, fontweight="bold", color="white",
        path_effects=[matplotlib.patheffects.withStroke(
            linewidth=2, foreground="black")])
    ax_map.set_axis_off()
    fig_map.tight_layout()

    # Name the map file using the snapshot date: MRSWAT_MapYYYYMMDD.jpg
    try:
        snap_dt = ref_date + dt.timedelta(days=best_time)
        map_name = f"MRSWAT_Map{snap_dt.strftime('%Y%m%d')}.jpg"
    except Exception:
        map_name = "MRSWAT_Map.jpg"

    map_path = OUTPUT_DIR / map_name
    fig_map.savefig(str(map_path), dpi=300, bbox_inches="tight",
                    facecolor="black", edgecolor="none")
    plt.close(fig_map)
    print(f"  ✓ Salinity map saved: {map_path}")


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 7B — ANIMATED SATELLITE SALINITY MAP
# ═══════════════════════════════════════════════════════════════════════════════

def make_salinity_map_animation(ref_date):
    """
    Animated satellite basemap showing bottom-layer salinity markers, 9 ppt
    and 0.25 ppt threshold circles, and infrastructure labels evolving
    through all model output time-steps.

    Static elements (drawn once):
      • Esri World Imagery basemap tiles
      • River centerline
      • RM labels (every 10 miles)
      • Infrastructure location labels
      • Colorbar

    Dynamic elements (updated per frame):
      • Scatter marker colours (bottom-layer salinity at integer RMs)
      • White solid circle  – upstream-most RM with bottom salinity > 9 ppt
      • Red dashed circle – upstream-most RM with surface salinity > 0.25 ppt
      • Associated text labels for both circles
      • Title date string

    Saves to <OUTPUT_DIR>/MRSWAT_Map_Animation.mp4 (or .gif fallback).
    """
    print("\n" + "=" * 70)
    print("  Generating Animated Satellite Salinity Map …")
    print("=" * 70)

    try:
        import contextily as ctx
    except ImportError:
        print("  WARNING: contextily not installed – skipping map animation.")
        return

    from matplotlib.patches import Circle as MplCircle

    # Resolve ffmpeg
    ffmpeg_exe = _resolve_ffmpeg()
    if ffmpeg_exe:
        matplotlib.rcParams["animation.ffmpeg_path"] = ffmpeg_exe
        print(f"  ffmpeg located: {ffmpeg_exe}")
    else:
        print("  WARNING: ffmpeg not found – GIF fallback will be used.")

    # ── Load all snapshots ────────────────────────────────────────────────
    snapshots = read_mrswat_output()
    if not snapshots:
        print("  No data found – skipping map animation.")
        return

    # ── Static geometry ─────────────────────────────────────────────────
    MAP_RM_MIN, MAP_RM_MAX = -19.0, 100.0
    marker_rms  = np.arange(int(np.ceil(MAP_RM_MIN)), int(np.floor(MAP_RM_MAX)) + 1)
    waypoints   = load_river_mile_waypoints()
    marker_lat, marker_lon = rm_to_latlon(marker_rms.astype(float), waypoints)
    marker_x, marker_y     = latlon_to_web_mercator(marker_lat, marker_lon)
    circ_radius = (marker_x.max() - marker_x.min()) * 0.006

    # Centerline from the first snapshot
    d0   = snapshots[0]["data"]
    m0   = (d0[:, 1] >= MAP_RM_MIN) & (d0[:, 1] <= MAP_RM_MAX)
    cl_lat, cl_lon = rm_to_latlon(d0[m0, 1], waypoints)
    cl_x, cl_y     = latlon_to_web_mercator(cl_lat, cl_lon)

    cmap, norm = build_salinity_cmap()

    # ── Build figure & static elements ─────────────────────────────────
    fig_anim, ax_anim = plt.subplots(figsize=(10, 14))
    ax_anim.plot(cl_x, cl_y, color="white", lw=1.0, alpha=0.40, zorder=1)

    # Initial scatter (colours will be overwritten every frame)
    sc = ax_anim.scatter(marker_x, marker_y, c=np.zeros(len(marker_rms)),
                         cmap=cmap, norm=norm, s=30,
                         edgecolors="k", linewidths=0.3, zorder=3)
    cbar = fig_anim.colorbar(sc, ax=ax_anim,
                             label="Bottom-Layer Salinity (ppt)",
                             shrink=0.55, pad=0.02)
    cbar.ax.yaxis.label.set_color("white")
    cbar.ax.tick_params(colors="white")

    # RM labels every 10 miles
    for i, rm_val in enumerate(marker_rms):
        if rm_val % 10 == 0:
            ax_anim.annotate(
                f"RM {rm_val:g}",
                xy=(marker_x[i], marker_y[i]),
                xytext=(6, 4), textcoords="offset points",
                fontsize=7, fontweight="bold", color="white",
                path_effects=[matplotlib.patheffects.withStroke(
                    linewidth=2, foreground="black")],
                zorder=4)

    # Infrastructure labels
    INFRA_LABELS = {
        19: "Boothville",
        49: "Port Sulphur/PALH",
        76: "Belle Chasse",
        81: "Dalcour",
        88: "St. Bernard",
        96: "Algiers",
    }
    for rm_val, label_text in INFRA_LABELS.items():
        idx_arr = np.where(marker_rms == rm_val)[0]
        if len(idx_arr) == 0:
            continue
        idx = idx_arr[0]
        ax_anim.annotate(
            label_text,
            xy=(marker_x[idx], marker_y[idx]),
            xytext=(8, -10), textcoords="offset points",
            fontsize=7.5, fontweight="bold", color="yellow",
            path_effects=[matplotlib.patheffects.withStroke(
                linewidth=2, foreground="black")],
            zorder=6)

    # Axis limits and basemap
    pad_x = (marker_x.max() - marker_x.min()) * 0.08
    pad_y = (marker_y.max() - marker_y.min()) * 0.05
    ax_anim.set_xlim(marker_x.min() - pad_x, marker_x.max() + pad_x)
    ax_anim.set_ylim(marker_y.min() - pad_y, marker_y.max() + pad_y)

    try:
        ctx.add_basemap(ax_anim, source=ctx.providers.Esri.WorldImagery,
                        zoom="auto", attribution=False)
        print("  Satellite basemap tiles loaded.")
    except Exception as e:
        print(f"  WARNING: Could not fetch satellite tiles ({e}).")

    title_obj = ax_anim.set_title(
        "", fontsize=13, fontweight="bold", color="white",
        path_effects=[matplotlib.patheffects.withStroke(
            linewidth=2, foreground="black")])
    ax_anim.set_axis_off()
    fig_anim.tight_layout()

    # ── Load main toe from output file ──────────────────────────────────
    toe_time_f, toe_loc_f, _ = load_toe_data()

    # ── Mutable artist holders for threshold circles / labels ───────────
    _artists = {"circ_white": None, "circ_red": None, "circ_alt": None,
                "label_9": None, "label_025": None, "label_alt": None,
                "pool_artists": []}

    # ── Per-frame callback ───────────────────────────────────────────────
    def _update_map(frame):
        snap  = snapshots[frame]
        d     = snap["data"]
        t_day = snap["time"]

        # Update title with date only (no time)
        snap_date = (ref_date + dt.timedelta(days=t_day)).strftime("%b %d %Y")
        title_obj.set_text(
            f"Mississippi River Salinity  \u00b7  {snap_date}\n"
            f"RM {int(np.ceil(MAP_RM_MIN))} (SW Pass outlet) to "
            f"RM {int(MAP_RM_MAX)}")

        # Interpolate salinity at integer river miles
        mask      = (d[:, 1] >= MAP_RM_MIN) & (d[:, 1] <= MAP_RM_MAX)
        rm_snap   = d[mask, 1]
        sal_bot   = d[mask, 6]
        sal_surf  = d[mask, 5]
        m_sal     = np.interp(marker_rms, rm_snap, sal_bot)
        m_surf    = np.interp(marker_rms, rm_snap, sal_surf)

        # Update scatter colours
        sc.set_array(m_sal)

        # -- Remove previous dynamic artists -------------------------------
        for key in ("circ_white", "circ_red", "circ_alt",
                    "label_9", "label_025", "label_alt"):
            if _artists[key] is not None:
                _artists[key].remove()
                _artists[key] = None
        for _pa in _artists["pool_artists"]:
            _pa.remove()
        _artists["pool_artists"] = []

        # -- White circle: main 9-ppt toe from file -----------------------
        _main_toe = float(np.interp(t_day, toe_time_f, toe_loc_f))
        if not np.isnan(_main_toe):
            idx_9 = int(np.argmin(np.abs(marker_rms - _main_toe)))
            circ_w = MplCircle(
                (marker_x[idx_9], marker_y[idx_9]),
                radius=circ_radius, transform=ax_anim.transData,
                facecolor="none", edgecolor="white", linewidth=1.5,
                zorder=6)
            ax_anim.add_patch(circ_w)
            _artists["circ_white"] = circ_w
            _artists["label_9"] = ax_anim.annotate(
                "9 ppt",
                xy=(marker_x[idx_9], marker_y[idx_9]),
                xytext=(8, 4), textcoords="offset points",
                fontsize=7.5, fontweight="bold", color="white",
                path_effects=[matplotlib.patheffects.withStroke(
                    linewidth=2, foreground="black")],
                zorder=7)

        # -- Alt toe + pool circles (conditional) --------------------------
        if ALT_TOE_ENABLED:
            _alt_rm, _pools = compute_connected_toe(d)
            # Orange dashed circle: alternative connected toe
            if not np.isnan(_alt_rm):
                idx_alt = int(np.argmin(np.abs(marker_rms - _alt_rm)))
                circ_a = MplCircle(
                    (marker_x[idx_alt], marker_y[idx_alt]),
                    radius=circ_radius, transform=ax_anim.transData,
                    facecolor="none", edgecolor="orange", linewidth=1.5,
                    linestyle="--", zorder=6)
                ax_anim.add_patch(circ_a)
                _artists["circ_alt"] = circ_a
                _artists["label_alt"] = ax_anim.annotate(
                    "alt toe",
                    xy=(marker_x[idx_alt], marker_y[idx_alt]),
                    xytext=(8, -12), textcoords="offset points",
                    fontsize=7, fontweight="bold", color="orange",
                    path_effects=[matplotlib.patheffects.withStroke(
                        linewidth=2, foreground="black")],
                    zorder=7)

            # Yellow dashed circles: salinity pools
            for _pool in _pools:
                _p_center = (_pool['lo'] + _pool['hi']) / 2.0
                _p_idx = int(np.argmin(np.abs(marker_rms - _p_center)))
                circ_p = MplCircle(
                    (marker_x[_p_idx], marker_y[_p_idx]),
                    radius=circ_radius, transform=ax_anim.transData,
                    facecolor="none", edgecolor="yellow", linewidth=1.5,
                    linestyle="--", zorder=6)
                ax_anim.add_patch(circ_p)
                _artists["pool_artists"].append(circ_p)

        # -- Red dashed circle: surface salinity > 0.25 ppt ----------------
        mask_025 = m_surf > 0.25
        if np.any(mask_025):
            idx_025 = int(np.max(np.where(mask_025)[0]))
            circ_r = MplCircle(
                (marker_x[idx_025], marker_y[idx_025]),
                radius=circ_radius, transform=ax_anim.transData,
                facecolor="none", edgecolor="red", linewidth=2.0,
                linestyle="--", zorder=7)
            ax_anim.add_patch(circ_r)
            _artists["circ_red"] = circ_r
            _artists["label_025"] = ax_anim.annotate(
                "0.25 ppt",
                xy=(marker_x[idx_025], marker_y[idx_025]),
                xytext=(-65, 6), textcoords="offset points",
                fontsize=7.5, fontweight="bold", color="red",
                path_effects=[matplotlib.patheffects.withStroke(
                    linewidth=2, foreground="black")],
                zorder=8)

        return (sc, title_obj)

    # ── Build and save animation ─────────────────────────────────────
    print(f"  Generating map animation with {len(snapshots)} frames …")
    ani = animation.FuncAnimation(fig_anim, _update_map,
                                  frames=len(snapshots),
                                  interval=200, blit=False)

    start_dt = ref_date + dt.timedelta(
        days=int(np.floor(snapshots[0]["time"])))
    end_dt = ref_date + dt.timedelta(
        days=int(np.floor(snapshots[-1]["time"])))
    out_base = (f"MRSWAT_Map_{start_dt.strftime('%Y%m%d')}"
                f"_{end_dt.strftime('%Y%m%d')}")

    mp4_path = OUTPUT_DIR / f"{out_base}.mp4"
    try:
        ani.save(str(mp4_path), writer="ffmpeg", fps=FPS, dpi=200,
                 savefig_kwargs={"facecolor": "black", "edgecolor": "none"})
        print(f"  ✓ Map animation saved: {mp4_path}")
    except Exception as e:
        print(f"  MP4 save failed ({e}); trying GIF fallback …")
        gif_path = OUTPUT_DIR / f"{out_base}.gif"
        try:
            ani.save(str(gif_path), writer="pillow", fps=FPS,
                     savefig_kwargs={"facecolor": "black", "edgecolor": "none"})
            print(f"  ✓ Map animation saved: {gif_path}")
        except Exception as e2:
            print(f"  GIF save also failed ({e2}).")
    plt.close(fig_anim)


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 8 — SALINITY ANIMATION
# ═══════════════════════════════════════════════════════════════════════════════

def make_animation(ref_date):
    """
    Produce an MP4 (or GIF fallback) animation of the time-evolving
    longitudinal salinity cross-section with boundary-condition time series.

    Animation layout
    -----------------
      TOP PANEL   – Pseudo-2-D salinity field via bilinear interpolation onto
                    a vertical grid, overlain by bed profile (saddlebrown),
                    water surface (cyan), and interface (white dashed).
      BOTTOM PANEL – Discharge (steelblue, left axis) and tidal stage
                     (orange, right axis) with a moving crimson marker.

    Saves to <OUTPUT_DIR>/MRSWAT_<START>_<END>.mp4 (or .gif).
    """
    print("\n" + "=" * 70)
    print("  Generating Salinity Animation …")
    print("=" * 70)

    # Resolve ffmpeg
    ffmpeg_exe = _resolve_ffmpeg()
    if ffmpeg_exe:
        matplotlib.rcParams["animation.ffmpeg_path"] = ffmpeg_exe
        print(f"  ffmpeg located: {ffmpeg_exe}")
    else:
        print("  WARNING: ffmpeg not found – GIF fallback will be used.")

    # ── Load data ────────────────────────────────────────────────────────
    disc_times, discharges = load_two_col_ts(str(DISCHARGE_FILE))
    stage_times, stages    = load_two_col_ts(str(STAGE_FILE))
    snapshots = read_mrswat_output()
    if not snapshots:
        print("  No data found – skipping animation.")
        return

    all_snap_times = [s["time"] for s in snapshots]
    t_min, t_max   = min(all_snap_times), max(all_snap_times)
    span_hw = max((t_max - t_min) * 0.008, 0.05)

    start_dt = ref_date + timedelta(days=int(np.floor(t_min)))
    end_dt   = ref_date + timedelta(days=int(np.floor(t_max)))
    out_base = f"MRSWAT_{start_dt.strftime('%Y%m%d')}_{end_dt.strftime('%Y%m%d')}"

    # ── Figure layout ────────────────────────────────────────────────────
    fig = plt.figure(figsize=(14, 10))
    gs  = gridspec.GridSpec(2, 1, height_ratios=[2, 1], hspace=0.42)
    ax     = fig.add_subplot(gs[0])
    ax_ts  = fig.add_subplot(gs[1])
    ax_ts2 = ax_ts.twinx()

    # ── Vertical limits (global scan) ────────────────────────────────────
    all_bed, all_wsel = [], []
    for s in snapshots:
        d = s["data"]
        m = d[:, 1] <= MAX_RIVER_MILE
        if np.any(m):
            all_bed.extend(d[m, 4])
            all_wsel.extend(d[m, 2])
    min_elev = np.floor(np.min(all_bed))  - 5
    max_elev = np.ceil(np.max(all_wsel))  + 5
    y_grid   = np.linspace(min_elev, max_elev, 150)

    cmap, norm = build_salinity_cmap()

    # ── Initial plot elements ────────────────────────────────────────────
    dummy_X, dummy_Y = np.meshgrid([0, 1], y_grid)
    mesh = ax.pcolormesh(dummy_X, dummy_Y, np.zeros_like(dummy_X),
                         cmap=cmap, norm=norm, shading="gouraud")
    fig.colorbar(mesh, ax=ax, label="Salinity (ppt)")

    line_bed,  = ax.plot([], [], color="saddlebrown", lw=2,
                         label="Bed", zorder=3)
    line_wsel, = ax.plot([], [], color="cyan", lw=2,
                         label="Water Surface", zorder=4)
    line_intf, = ax.plot([], [], color="white", lw=1.5, ls="--",
                         label="Interface", zorder=5)
    poly_bed   = ax.fill_between([], [], color="saddlebrown")

    ax.set_xlim(0, MAX_RIVER_MILE)
    ax.set_ylim(min_elev, max_elev)
    ax.set_xlabel("River Mile (from Head of Passes)")
    ax.set_ylabel("Elevation (ft, NAVD 88)")
    ax.legend(loc="upper right")
    ax.grid(True, alpha=0.3)
    title_text = ax.set_title("", fontsize=13, fontweight="bold")

    # ── Bottom panel (static) ────────────────────────────────────────────
    # Discharge – primary hydrodynamic forcing controlling wedge position.
    # Typical low-flow: ~150 000–200 000 cfs; major floods >1 200 000 cfs.
    ax_ts.plot(disc_times, discharges, color="steelblue", lw=1.5,
               label="Discharge (cfs)")
    # Tidal stage at or near Head of Passes; ~0.5–1.5 ft semi-diurnal range.
    ax_ts2.plot(stage_times, stages, color="darkorange", lw=1.0,
                label="Stage (ft)")
    ax_ts.set_xlim(t_min, t_max)
    ax_ts.set_xlabel("Elapsed Days Since Reference Date")
    ax_ts.set_ylabel("Discharge (cfs)", color="steelblue")
    ax_ts2.set_ylabel("Downstream Stage (ft)", color="darkorange")
    ax_ts.tick_params(axis="y",  labelcolor="steelblue")
    ax_ts2.tick_params(axis="y", labelcolor="darkorange")
    ax_ts.grid(True, alpha=0.3)
    ax_ts.set_title("Open Boundary Conditions: River Discharge & Tidal Stage",
                    fontsize=11)
    l1, lb1 = ax_ts.get_legend_handles_labels()
    l2, lb2 = ax_ts2.get_legend_handles_labels()
    ax_ts.legend(l1 + l2, lb1 + lb2, loc="upper right", fontsize=9)

    vline = ax_ts.axvline(x=t_min, color="crimson", lw=1.5,
                          alpha=0.90, zorder=6)
    vspan = [ax_ts.axvspan(t_min, t_min + span_hw,
                           color="crimson", alpha=0.20, zorder=5)]

    prev_int_day = [None]
    cached_title = [""]

    # ── Per-frame callback ───────────────────────────────────────────────
    def update(frame):
        """Redraw all animated elements for a single time-step frame."""
        nonlocal mesh, poly_bed

        s = snapshots[frame]
        d = s["data"]
        t = s["time"]

        int_day = int(np.floor(t))
        if int_day != prev_int_day[0]:
            cd = ref_date + timedelta(days=int_day)
            cached_title[0] = (f"MRSWAT Salinity Profile  ·  "
                               f"{cd.strftime('%b %d %Y')}  ·  "
                               f"Day {int_day} since reference")
            prev_int_day[0] = int_day
        title_text.set_text(cached_title[0])

        msk = d[:, 1] <= MAX_RIVER_MILE
        if np.any(msk):
            rm_i = d[msk, 1]
            idx  = np.argsort(rm_i)
            rm_i   = rm_i[idx]
            bed_i  = d[msk, 4][idx]
            wsel_i = d[msk, 2][idx]
            intf_i = d[msk, 3][idx]
            ss_i   = d[msk, 5][idx]
            sb_i   = d[msk, 6][idx]

            X_new, Y_new, Z_new = interpolate_profile(
                rm_i, bed_i, wsel_i, intf_i, ss_i, sb_i, y_grid)
            mesh.remove()
            mesh = ax.pcolormesh(X_new, Y_new, Z_new,
                                 cmap=cmap, norm=norm,
                                 shading="gouraud", zorder=1)
            line_bed.set_data(rm_i, bed_i)
            line_wsel.set_data(rm_i, wsel_i)
            line_intf.set_data(rm_i, intf_i)
            poly_bed.remove()
            poly_bed = ax.fill_between(rm_i, min_elev, bed_i,
                                       color="saddlebrown", alpha=1.0, zorder=2)

        vline.set_xdata([t, t])
        vspan[0].remove()
        vspan[0] = ax_ts.axvspan(max(t_min, t - span_hw),
                                 min(t_max, t + span_hw),
                                 color="crimson", alpha=0.20, zorder=5)
        return (line_bed, line_wsel, line_intf, mesh, poly_bed,
                title_text, vline)

    # ── Build and save ───────────────────────────────────────────────────
    print(f"  Generating animation with {len(snapshots)} frames …")
    ani = animation.FuncAnimation(fig, update, frames=len(snapshots),
                                  interval=100, blit=False)
    mp4_path = OUTPUT_DIR / f"{out_base}.mp4"
    try:
        ani.save(str(mp4_path), writer="ffmpeg", fps=FPS, dpi=150)
        print(f"  ✓ Animation saved: {mp4_path}")
    except Exception as e:
        print(f"  MP4 save failed ({e}); trying GIF fallback …")
        gif_path = OUTPUT_DIR / f"{out_base}.gif"
        try:
            ani.save(str(gif_path), writer="pillow", fps=FPS)
            print(f"  ✓ Animation saved: {gif_path}")
        except Exception as e2:
            print(f"  GIF save also failed ({e2}); displaying interactively.")
            plt.show()
            return
    plt.close(fig)


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 8B — CALIBRATION ANIMATION  (Longitudinal salinity + observed CTD)
# ═══════════════════════════════════════════════════════════════════════════════

# ---------- CTD field-data helpers (needed only by calibration animation) -----

def _haversine_ft(lat1, lon1, lat2, lon2):
    """Great-circle distance between two points in **feet**."""
    R = 20_902_231.0  # Earth radius in feet
    dlat = np.radians(lat2 - lat1)
    dlon = np.radians(lon2 - lon1)
    a = (np.sin(dlat / 2) ** 2
         + np.cos(np.radians(lat1)) * np.cos(np.radians(lat2))
         * np.sin(dlon / 2) ** 2)
    return R * 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))


def latlon_to_rm(lat, lon, waypoints):
    """
    Convert a (lat, lon) position to an approximate river mile using the
    USACE waypoint table.  Finds the closest waypoint by Haversine distance,
    then refines by projecting onto the two adjacent RM segments.
    """
    dists = np.array([_haversine_ft(lat, lon, wp[1], wp[2]) for wp in waypoints])
    i_best = int(np.argmin(dists))
    best_rm = waypoints[i_best, 0]
    for i_start in [max(0, i_best - 1), i_best]:
        i_end = i_start + 1
        if i_end >= len(waypoints):
            continue
        ax_, ay_ = waypoints[i_start, 1], waypoints[i_start, 2]
        bx_, by_ = waypoints[i_end, 1],   waypoints[i_end, 2]
        abx, aby = bx_ - ax_, by_ - ay_
        apx, apy = lat - ax_, lon - ay_
        ab2 = abx * abx + aby * aby
        if ab2 < 1e-14:
            continue
        t = max(0.0, min(1.0, (apx * abx + apy * aby) / ab2))
        rm_proj = waypoints[i_start, 0] + t * (waypoints[i_end, 0] - waypoints[i_start, 0])
        plat = ax_ + t * abx
        plon = ay_ + t * aby
        d = _haversine_ft(lat, lon, plat, plon)
        if d < dists[i_best] + 100:
            best_rm = rm_proj
            break
    return best_rm


def _parse_ctd_full(filepath):
    """
    Parse a single CTD cast CSV and return the full depth/salinity profile.

    Returns
    -------
    dict  with keys:
        'lat', 'lon' : float
        'depths'     : np.ndarray – depth values (ft, positive down)
        'salinities' : np.ndarray – salinity values (ppt)
    or None if the file cannot be parsed.
    """
    import re as _re
    try:
        with open(str(filepath), "r", encoding="utf-8", errors="replace") as f:
            lines = f.readlines()
    except OSError:
        return None

    lat = lon = None
    for line in lines[:28]:
        if line.startswith("% Start latitude"):
            try:
                lat = float(line.split(",")[1].strip())
            except (IndexError, ValueError):
                pass
        elif line.startswith("% Start longitude"):
            try:
                lon = float(line.split(",")[1].strip())
            except (IndexError, ValueError):
                pass
    if lat is None or lon is None:
        return None

    # Column header at line 29 (0-indexed 28); data from line 30 onward.
    # Column 1 = Depth (ft), Column 5 = Salinity (PSS).
    depths, sals = [], []
    for line in lines[29:]:
        parts = line.strip().split(",")
        if len(parts) < 6:
            continue
        try:
            depths.append(float(parts[1]))
            sals.append(float(parts[5]))
        except (ValueError, IndexError):
            continue
    if not depths:
        return None
    return {"lat": lat, "lon": lon,
            "depths": np.array(depths), "salinities": np.array(sals)}


def load_field_observations(date_start, date_end, waypoints):
    """
    Scan the CTD field-data directory for measurement days that fall within
    [date_start, date_end].  For each day, return a list of parsed casts
    with depth-salinity profiles and computed river-mile positions.

    Returns
    -------
    dict  :  {datetime.date: list-of-cast-dicts}
        Each cast dict has keys:  'rm', 'depths', 'salinities'
    """
    import re as _re
    obs = {}  # date -> [cast, ...]
    if not FIELD_DATA_ROOT.exists():
        print(f"  WARNING: Field-data root not found: {FIELD_DATA_ROOT}")
        return obs

    for year_dir in sorted(FIELD_DATA_ROOT.iterdir()):
        if not year_dir.is_dir():
            continue
        try:
            year_val = int(year_dir.name)
        except ValueError:
            continue
        if year_val < date_start.year or year_val > date_end.year:
            continue

        for sww_dir in sorted(year_dir.iterdir()):
            if not sww_dir.is_dir():
                continue
            m = _re.match(r"SWW_(\d{8})", sww_dir.name)
            if not m:
                continue
            try:
                obs_date = dt.datetime.strptime(m.group(1), "%Y%m%d").date()
            except ValueError:
                continue
            if obs_date < date_start.date() or obs_date > date_end.date():
                continue

            casts = []
            for csv_file in sorted(sww_dir.glob("*.csv")):
                if csv_file.name.lower() == "summary.csv":
                    continue
                result = _parse_ctd_full(csv_file)
                if result is not None:
                    result["rm"] = latlon_to_rm(result["lat"], result["lon"], waypoints)
                    casts.append(result)
            if casts:
                obs[obs_date] = casts
    return obs


def make_calibration_animation(ref_date):
    """
    Produce an MP4 (or GIF fallback) animation identical to make_animation()
    but with observed CTD data overlaid.

    When the animation date matches a field-observation day, colour-filled
    boxes (1 RM wide × measured depth range) are superimposed on the
    salinity cross-section using the same colourmap.  The boxes persist
    until new observations appear or until 10 days have elapsed — whichever
    comes first.

    Saves to <OUTPUT_DIR>/MRSWAT_Calibration_<START>_<END>.mp4 (or .gif).
    """
    print("\n" + "=" * 70)
    print("  Generating Calibration Animation (Longitudinal + Observed CTD) ...")
    print("=" * 70)

    # Resolve ffmpeg
    ffmpeg_exe = _resolve_ffmpeg()
    if ffmpeg_exe:
        matplotlib.rcParams["animation.ffmpeg_path"] = ffmpeg_exe
        print(f"  ffmpeg located: {ffmpeg_exe}")
    else:
        print("  WARNING: ffmpeg not found -- GIF fallback will be used.")

    # -- Load modeled data -----------------------------------------------------
    disc_times, discharges = load_two_col_ts(str(DISCHARGE_FILE))
    stage_times, stages    = load_two_col_ts(str(STAGE_FILE))
    snapshots = read_mrswat_output()
    if not snapshots:
        print("  No data found -- skipping calibration animation.")
        return

    all_snap_times = [s["time"] for s in snapshots]
    t_min, t_max   = min(all_snap_times), max(all_snap_times)
    span_hw = max((t_max - t_min) * 0.008, 0.05)

    start_dt = ref_date + timedelta(days=int(np.floor(t_min)))
    end_dt   = ref_date + timedelta(days=int(np.floor(t_max)))
    out_base = f"MRSWAT_Calibration_{start_dt.strftime('%Y%m%d')}_{end_dt.strftime('%Y%m%d')}"

    # -- Load observed field data -----------------------------------------------
    waypoints = load_river_mile_waypoints()
    obs_by_date = load_field_observations(start_dt, end_dt, waypoints)
    obs_dates_sorted = sorted(obs_by_date.keys())
    print(f"  Observation days loaded: {len(obs_dates_sorted)}")

    # Convert observation dates to elapsed days since ref_date
    obs_elapsed = {d: (dt.datetime.combine(d, dt.time()) - ref_date).total_seconds() / 86400.0
                   for d in obs_dates_sorted}

    # -- Figure layout ---------------------------------------------------------
    fig = plt.figure(figsize=(14, 10))
    gs  = gridspec.GridSpec(2, 1, height_ratios=[2, 1], hspace=0.42)
    ax     = fig.add_subplot(gs[0])
    ax_ts  = fig.add_subplot(gs[1])
    ax_ts2 = ax_ts.twinx()

    # -- Vertical limits (global scan) -----------------------------------------
    all_bed, all_wsel = [], []
    for s in snapshots:
        d = s["data"]
        m = d[:, 1] <= MAX_RIVER_MILE
        if np.any(m):
            all_bed.extend(d[m, 4])
            all_wsel.extend(d[m, 2])
    min_elev = np.floor(np.min(all_bed))  - 5
    max_elev = np.ceil(np.max(all_wsel))  + 5
    y_grid   = np.linspace(min_elev, max_elev, 150)

    cmap, norm = build_salinity_cmap()

    # -- Initial plot elements -------------------------------------------------
    dummy_X, dummy_Y = np.meshgrid([0, 1], y_grid)
    mesh = ax.pcolormesh(dummy_X, dummy_Y, np.zeros_like(dummy_X),
                         cmap=cmap, norm=norm, shading="gouraud")
    fig.colorbar(mesh, ax=ax, label="Salinity (ppt)")

    line_bed,  = ax.plot([], [], color="saddlebrown", lw=2,
                         label="Bed", zorder=3)
    line_wsel, = ax.plot([], [], color="cyan", lw=2,
                         label="Water Surface", zorder=4)
    line_intf, = ax.plot([], [], color="white", lw=1.5, ls="--",
                         label="Interface", zorder=5)
    poly_bed   = ax.fill_between([], [], color="saddlebrown")

    ax.set_xlim(0, MAX_RIVER_MILE)
    ax.set_ylim(min_elev, max_elev)
    ax.set_xlabel("River Mile (from Head of Passes)")
    ax.set_ylabel("Elevation (ft, NAVD 88)")
    ax.legend(loc="upper right")
    ax.grid(True, alpha=0.3)
    title_text = ax.set_title("", fontsize=13, fontweight="bold")

    # -- Bottom panel (static) -------------------------------------------------
    ax_ts.plot(disc_times, discharges, color="steelblue", lw=1.5,
               label="Discharge (cfs)")
    ax_ts2.plot(stage_times, stages, color="darkorange", lw=1.0,
                label="Stage (ft)")
    ax_ts.set_xlim(t_min, t_max)
    ax_ts.set_xlabel("Elapsed Days Since Reference Date")
    ax_ts.set_ylabel("Discharge (cfs)", color="steelblue")
    ax_ts2.set_ylabel("Downstream Stage (ft)", color="darkorange")
    ax_ts.tick_params(axis="y",  labelcolor="steelblue")
    ax_ts2.tick_params(axis="y", labelcolor="darkorange")
    ax_ts.grid(True, alpha=0.3)
    ax_ts.set_title("Open Boundary Conditions: River Discharge & Tidal Stage",
                    fontsize=11)
    l1, lb1 = ax_ts.get_legend_handles_labels()
    l2, lb2 = ax_ts2.get_legend_handles_labels()
    ax_ts.legend(l1 + l2, lb1 + lb2, loc="upper right", fontsize=9)

    vline = ax_ts.axvline(x=t_min, color="crimson", lw=1.5,
                          alpha=0.90, zorder=6)
    vspan = [ax_ts.axvspan(t_min, t_min + span_hw,
                           color="crimson", alpha=0.20, zorder=5)]

    prev_int_day = [None]
    cached_title = [""]

    # State for observed-data overlay boxes
    active_obs_patches = []      # list of matplotlib patches to remove later
    active_obs_day     = [None]  # the datetime.date currently displayed
    active_obs_elapsed = [None]  # elapsed day when current obs were shown

    OBS_PERSIST_DAYS = 10  # boxes disappear after this many days unless replaced

    def _clear_obs_patches():
        """Remove all currently displayed observation boxes."""
        for p in active_obs_patches:
            p.remove()
        active_obs_patches.clear()
        active_obs_day[0] = None
        active_obs_elapsed[0] = None

    def _draw_obs_patches(obs_date, cmap, norm):
        """Draw observed CTD cast boxes for a given observation day."""
        _clear_obs_patches()
        casts = obs_by_date[obs_date]
        half_width = 0.5  # RM -- total box width = 1 RM

        for cast in casts:
            rm = cast["rm"]
            if rm < 0 or rm > MAX_RIVER_MILE:
                continue
            depths = cast["depths"]        # ft, positive down from surface
            sals   = cast["salinities"]
            if len(depths) < 2:
                continue

            # Convert depth-below-surface to elevation.
            # We need an approximate water surface at this RM.
            # Use the current snapshot's water surface elevation.
            # (interpolated from the modeled profile)
            # -- this is passed in via the closure over the current frame --
            # For now, use a fixed approach: estimate WSEL from model data
            # of the current frame (will be set before calling this).

            # Sort by depth (ascending)
            order = np.argsort(depths)
            depths = depths[order]
            sals   = sals[order]

            # Build vertical colour strips (each depth interval gets one patch)
            rm_lo = rm - half_width
            rm_hi = rm + half_width

            for k in range(len(depths) - 1):
                d_top = depths[k]
                d_bot = depths[k + 1]
                sal_k = 0.5 * (sals[k] + sals[k + 1])  # average of interval

                # Elevation: depth is measured positive downward from surface.
                # We need the water surface elevation (WSEL) at this RM.
                # We'll store it in _current_wsel_interp (set per frame).
                elev_top = _current_wsel_interp[0] - d_top
                elev_bot = _current_wsel_interp[0] - d_bot

                color = cmap(norm(sal_k))
                rect = plt.Rectangle(
                    (rm_lo, elev_bot),
                    rm_hi - rm_lo,           # width
                    elev_top - elev_bot,      # height
                    facecolor=color, edgecolor="black", linewidth=0.3,
                    alpha=0.9, zorder=8
                )
                ax.add_patch(rect)
                active_obs_patches.append(rect)

        active_obs_day[0] = obs_date
        active_obs_elapsed[0] = _current_elapsed[0]

    # Shared mutable state for the update closure
    _current_wsel_interp = [0.0]
    _current_elapsed     = [0.0]

    # -- Per-frame callback ----------------------------------------------------
    def update(frame):
        nonlocal mesh, poly_bed

        s = snapshots[frame]
        d = s["data"]
        t = s["time"]
        _current_elapsed[0] = t

        int_day = int(np.floor(t))
        current_date = (ref_date + timedelta(days=int_day)).date()

        if int_day != prev_int_day[0]:
            cd = ref_date + timedelta(days=int_day)
            cached_title[0] = (f"MRSWAT Salinity Profile  \u00b7  "
                               f"{cd.strftime('%b %d %Y')}  \u00b7  "
                               f"Day {int_day} since reference")
            prev_int_day[0] = int_day
        title_text.set_text(cached_title[0])

        msk = d[:, 1] <= MAX_RIVER_MILE
        if np.any(msk):
            rm_i = d[msk, 1]
            idx  = np.argsort(rm_i)
            rm_i   = rm_i[idx]
            bed_i  = d[msk, 4][idx]
            wsel_i = d[msk, 2][idx]
            intf_i = d[msk, 3][idx]
            ss_i   = d[msk, 5][idx]
            sb_i   = d[msk, 6][idx]

            X_new, Y_new, Z_new = interpolate_profile(
                rm_i, bed_i, wsel_i, intf_i, ss_i, sb_i, y_grid)
            mesh.remove()
            mesh = ax.pcolormesh(X_new, Y_new, Z_new,
                                 cmap=cmap, norm=norm,
                                 shading="gouraud", zorder=1)
            line_bed.set_data(rm_i, bed_i)
            line_wsel.set_data(rm_i, wsel_i)
            line_intf.set_data(rm_i, intf_i)
            poly_bed.remove()
            poly_bed = ax.fill_between(rm_i, min_elev, bed_i,
                                       color="saddlebrown", alpha=1.0, zorder=2)

            # Provide WSEL interpolation for obs-overlay at any RM
            # (simple: use the average WSEL for the cast RM lookup)
            _current_wsel_interp[0] = float(np.interp(0, [0], wsel_i[:1]))
            # Actually, store the arrays so we can look up per-cast RM
            _current_rm_arr   = rm_i
            _current_wsel_arr = wsel_i

            # --- Observation overlay logic ---
            # Determine which obs day to show
            show_new_obs = False
            if current_date in obs_by_date:
                # New obs day arrived
                if active_obs_day[0] != current_date:
                    show_new_obs = True

            # Check expiration (10-day rule)
            if (active_obs_day[0] is not None
                    and active_obs_elapsed[0] is not None
                    and (t - active_obs_elapsed[0]) > OBS_PERSIST_DAYS
                    and not show_new_obs):
                _clear_obs_patches()

            if show_new_obs:
                # For each cast, compute WSEL at its RM from current model
                # We override _current_wsel_interp per-cast inside the drawer
                casts_today = obs_by_date[current_date]
                _clear_obs_patches()
                half_width = 0.5
                for cast in casts_today:
                    rm = cast["rm"]
                    if rm < 0 or rm > MAX_RIVER_MILE:
                        continue
                    depths_ = cast["depths"]
                    sals_   = cast["salinities"]
                    if len(depths_) < 2:
                        continue
                    # Water surface at this cast's RM
                    wsel_at_rm = float(np.interp(rm, _current_rm_arr, _current_wsel_arr))

                    order_ = np.argsort(depths_)
                    depths_ = depths_[order_]
                    sals_   = sals_[order_]

                    rm_lo = rm - half_width
                    rm_hi = rm + half_width
                    for k in range(len(depths_) - 1):
                        d_top = depths_[k]
                        d_bot = depths_[k + 1]
                        sal_k = 0.5 * (sals_[k] + sals_[k + 1])
                        elev_top = wsel_at_rm - d_top
                        elev_bot = wsel_at_rm - d_bot
                        color = cmap(norm(sal_k))
                        rect = plt.Rectangle(
                            (rm_lo, elev_bot),
                            rm_hi - rm_lo,
                            elev_top - elev_bot,
                            facecolor=color, edgecolor="black", linewidth=0.3,
                            alpha=0.9, zorder=8
                        )
                        ax.add_patch(rect)
                        active_obs_patches.append(rect)
                active_obs_day[0] = current_date
                active_obs_elapsed[0] = t

        vline.set_xdata([t, t])
        vspan[0].remove()
        vspan[0] = ax_ts.axvspan(max(t_min, t - span_hw),
                                 min(t_max, t + span_hw),
                                 color="crimson", alpha=0.20, zorder=5)
        return (line_bed, line_wsel, line_intf, mesh, poly_bed,
                title_text, vline)

    # -- Build and save --------------------------------------------------------
    print(f"  Generating calibration animation with {len(snapshots)} frames ...")
    ani = animation.FuncAnimation(fig, update, frames=len(snapshots),
                                  interval=100, blit=False)
    mp4_path = OUTPUT_DIR / f"{out_base}.mp4"
    try:
        ani.save(str(mp4_path), writer="ffmpeg", fps=FPS, dpi=150)
        print(f"  Animation saved: {mp4_path}")
    except Exception as e:
        print(f"  MP4 save failed ({e}); trying GIF fallback ...")
        gif_path = OUTPUT_DIR / f"{out_base}.gif"
        try:
            ani.save(str(gif_path), writer="pillow", fps=FPS)
            print(f"  Animation saved: {gif_path}")
        except Exception as e2:
            print(f"  GIF save also failed ({e2}); displaying interactively.")
            plt.show()
            return
    plt.close(fig)


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 8C — COMBINED 4-QUADRANT ANIMATION
# ═══════════════════════════════════════════════════════════════════════════════

def make_combined_animation(ref_date):
    """
    Combined four-quadrant animation on a single black-background figure:

      Left column (full height) – Animated satellite salinity map with
          coloured RM markers, threshold circles (9 ppt white solid,
          0.25 ppt red dashed) and infrastructure labels.
      Top-right quadrant  – Longitudinal salinity cross-section profile
          showing the two-layer structure through time.
      Bottom-right quadrant – River discharge and tidal stage time series
          with a moving crimson marker.

      A single centred super-title at the top updates every frame with the
      current date and elapsed days after the reference date.

    Saves to <OUTPUT_DIR>/MRSWAT_Combined_<START>_<END>.mp4 (or .gif).
    """
    print("\n" + "=" * 70)
    print("  Generating Combined 4-Quadrant Animation …")
    print("=" * 70)

    try:
        import contextily as ctx
    except ImportError:
        print("  WARNING: contextily not installed – satellite basemap will be blank.")
        ctx = None

    from matplotlib.patches import Circle as MplCircle

    # ── Resolve ffmpeg ────────────────────────────────────────────────
    ffmpeg_exe = _resolve_ffmpeg()
    if ffmpeg_exe:
        matplotlib.rcParams["animation.ffmpeg_path"] = ffmpeg_exe
        print(f"  ffmpeg located: {ffmpeg_exe}")
    else:
        print("  WARNING: ffmpeg not found – GIF fallback will be used.")

    # ── Load data ────────────────────────────────────────────────
    disc_times, discharges = load_two_col_ts(str(DISCHARGE_FILE))
    stage_times, stages    = load_two_col_ts(str(STAGE_FILE))
    snapshots = read_mrswat_output()
    if not snapshots:
        print("  No data found – skipping combined animation.")
        return

    all_snap_times = [s["time"] for s in snapshots]
    t_min, t_max   = min(all_snap_times), max(all_snap_times)
    span_hw = max((t_max - t_min) * 0.008, 0.05)

    start_dt = ref_date + dt.timedelta(days=int(np.floor(t_min)))
    end_dt   = ref_date + dt.timedelta(days=int(np.floor(t_max)))

    # ── Figure & GridSpec layout ──────────────────────────────────────
    # Left col  (both rows) : satellite map
    # Right col top         : salinity cross-section profile
    # Right col bottom      : discharge / stage time series
    fig = plt.figure(figsize=(22, 14), facecolor="black")
    gs  = gridspec.GridSpec(
        2, 2,
        height_ratios=[2, 1],
        width_ratios=[1, 1],
        hspace=0.30,
        wspace=0.18,
        figure=fig,
    )
    ax_map  = fig.add_subplot(gs[:, 0])   # left, full-height
    ax_prof = fig.add_subplot(gs[0, 1])   # top-right
    ax_ts   = fig.add_subplot(gs[1, 1])   # bottom-right
    ax_ts2  = ax_ts.twinx()

    # ── MAP — static elements ─────────────────────────────────────────
    MAP_RM_MIN, MAP_RM_MAX = -19.0, 100.0
    marker_rms  = np.arange(int(np.ceil(MAP_RM_MIN)), int(np.floor(MAP_RM_MAX)) + 1)
    waypoints   = load_river_mile_waypoints()
    marker_lat, marker_lon = rm_to_latlon(marker_rms.astype(float), waypoints)
    marker_x, marker_y     = latlon_to_web_mercator(marker_lat, marker_lon)
    circ_radius = (marker_x.max() - marker_x.min()) * 0.006

    d0 = snapshots[0]["data"]
    m0 = (d0[:, 1] >= MAP_RM_MIN) & (d0[:, 1] <= MAP_RM_MAX)
    cl_lat, cl_lon = rm_to_latlon(d0[m0, 1], waypoints)
    cl_x, cl_y     = latlon_to_web_mercator(cl_lat, cl_lon)

    cmap, norm = build_salinity_cmap()

    ax_map.plot(cl_x, cl_y, color="white", lw=1.0, alpha=0.40, zorder=1)
    sc = ax_map.scatter(marker_x, marker_y, c=np.zeros(len(marker_rms)),
                        cmap=cmap, norm=norm, s=30,
                        edgecolors="k", linewidths=0.3, zorder=3)
    cbar_map = fig.colorbar(sc, ax=ax_map,
                            label="Bottom-Layer Salinity (ppt)",
                            shrink=0.40, pad=0.02)
    cbar_map.ax.yaxis.label.set_color("white")
    cbar_map.ax.tick_params(colors="white")

    for i, rm_val in enumerate(marker_rms):
        if rm_val % 10 == 0:
            ax_map.annotate(
                f"RM {rm_val:g}",
                xy=(marker_x[i], marker_y[i]),
                xytext=(6, 4), textcoords="offset points",
                fontsize=7, fontweight="bold", color="white",
                path_effects=[matplotlib.patheffects.withStroke(
                    linewidth=2, foreground="black")],
                zorder=4)

    INFRA_LABELS = {
        19: "Boothville",
        49: "Port Sulphur/PALH",
        76: "Belle Chasse",
        81: "Dalcour",
        88: "St. Bernard",
        96: "Algiers",
    }
    for rm_val, label_text in INFRA_LABELS.items():
        idx_arr = np.where(marker_rms == rm_val)[0]
        if len(idx_arr) == 0:
            continue
        ax_map.annotate(
            label_text,
            xy=(marker_x[idx_arr[0]], marker_y[idx_arr[0]]),
            xytext=(8, -10), textcoords="offset points",
            fontsize=7.5, fontweight="bold", color="yellow",
            path_effects=[matplotlib.patheffects.withStroke(
                linewidth=2, foreground="black")],
            zorder=6)

    pad_x = (marker_x.max() - marker_x.min()) * 0.08
    pad_y = (marker_y.max() - marker_y.min()) * 0.05
    ax_map.set_xlim(marker_x.min() - pad_x, marker_x.max() + pad_x)
    ax_map.set_ylim(marker_y.min() - pad_y, marker_y.max() + pad_y)

    if ctx is not None:
        try:
            ctx.add_basemap(ax_map, source=ctx.providers.Esri.WorldImagery,
                            zoom="auto", attribution=False)
            print("  Satellite basemap tiles loaded.")
        except Exception as _e:
            print(f"  WARNING: Could not fetch satellite tiles ({_e}).")

    ax_map.set_axis_off()

    # ── Load main toe from output file ──────────────────────────────────
    toe_time_f, toe_loc_f, _ = load_toe_data()

    _map_artists = {"circ_white": None, "circ_red": None, "circ_alt": None,
                    "label_9": None, "label_025": None, "label_alt": None,
                    "pool_artists": []}

    # ── PROFILE — static elements ──────────────────────────────────────
    all_bed, all_wsel = [], []
    for s in snapshots:
        _d = s["data"]
        _m = _d[:, 1] <= MAX_RIVER_MILE
        if np.any(_m):
            all_bed.extend(_d[_m, 4])
            all_wsel.extend(_d[_m, 2])
    min_elev = np.floor(np.min(all_bed)) - 5
    max_elev = np.ceil(np.max(all_wsel)) + 5
    y_grid   = np.linspace(min_elev, max_elev, 150)

    _dummy_X, _dummy_Y = np.meshgrid([0, 1], y_grid)
    mesh_holder = [
        ax_prof.pcolormesh(_dummy_X, _dummy_Y, np.zeros_like(_dummy_X),
                           cmap=cmap, norm=norm, shading="gouraud")
    ]
    fig.colorbar(mesh_holder[0], ax=ax_prof, label="Salinity (ppt)", shrink=0.80)

    line_bed,  = ax_prof.plot([], [], color="saddlebrown", lw=2,
                              label="Bed", zorder=3)
    line_wsel, = ax_prof.plot([], [], color="cyan", lw=2,
                              label="Water Surface", zorder=4)
    line_intf, = ax_prof.plot([], [], color="white", lw=1.5, ls="--",
                              label="Interface", zorder=5)
    poly_bed_holder = [ax_prof.fill_between([], [], color="saddlebrown")]

    ax_prof.set_xlim(0, MAX_RIVER_MILE)
    ax_prof.set_ylim(min_elev, max_elev)
    ax_prof.set_xlabel("River Mile (from Head of Passes)", color="white")
    ax_prof.set_ylabel("Elevation (ft, NAVD 88)", color="white")
    ax_prof.tick_params(colors="white")
    for sp in ax_prof.spines.values():
        sp.set_color("white")
    ax_prof.legend(loc="upper right", fontsize=8,
                   facecolor="#222222", labelcolor="white")
    ax_prof.grid(True, alpha=0.3)
    ax_prof.set_facecolor("#111111")
    ax_prof.set_title("Salinity Cross-Section Profile", color="white", fontsize=10)

    # ── TIME SERIES — static elements ────────────────────────────────
    ax_ts.plot(disc_times, discharges, color="steelblue", lw=1.5,
               label="Discharge (cfs)")
    ax_ts2.plot(stage_times, stages, color="darkorange", lw=1.0,
                label="Stage (ft)")
    ax_ts.set_xlim(t_min, t_max)
    ax_ts.set_xlabel("Elapsed Days Since Reference Date", color="white")
    ax_ts.set_ylabel("Discharge (cfs)", color="steelblue")
    ax_ts2.set_ylabel("Downstream Stage (ft)", color="darkorange")
    ax_ts.tick_params(axis="x", colors="white")
    ax_ts.tick_params(axis="y", labelcolor="steelblue")
    ax_ts2.tick_params(axis="y", labelcolor="darkorange")
    for _sp_name, _color in [("bottom", "white"), ("top", "white"),
                              ("left", "steelblue"), ("right", "darkorange")]:
        ax_ts.spines[_sp_name].set_color(_color)
    ax_ts2.spines["right"].set_color("darkorange")
    ax_ts.grid(True, alpha=0.3)
    ax_ts.set_facecolor("#111111")
    ax_ts.set_title("River Discharge & Tidal Stage", color="white", fontsize=10)
    _l1, _lb1 = ax_ts.get_legend_handles_labels()
    _l2, _lb2 = ax_ts2.get_legend_handles_labels()
    ax_ts.legend(_l1 + _l2, _lb1 + _lb2, loc="upper right", fontsize=8,
                 facecolor="#222222", labelcolor="white")

    vline = ax_ts.axvline(x=t_min, color="crimson", lw=1.5,
                          alpha=0.90, zorder=6)
    vspan = [ax_ts.axvspan(t_min, t_min + span_hw,
                           color="crimson", alpha=0.20, zorder=5)]

    # ── Centred super-title (updated each frame) ──────────────────────
    sup = fig.suptitle("", fontsize=14, fontweight="bold",
                       color="white", y=0.998)

    # ── Per-frame callback ──────────────────────────────────────────
    def _update_combined(frame):
        snap  = snapshots[frame]
        d     = snap["data"]
        t_day = snap["time"]
        int_day = int(np.floor(t_day))

        # — Super-title: date + elapsed days ————————————————————
        snap_date = (ref_date + dt.timedelta(days=int_day)).strftime("%b %d %Y")
        sup.set_text(
            f"Mississippi River Salinity  ·  {snap_date}  ·  "
            f"Day {int_day} after reference date")

        # — Map: scatter colours + threshold circles ————————————————
        _mask = (d[:, 1] >= MAP_RM_MIN) & (d[:, 1] <= MAP_RM_MAX)
        _rm   = d[_mask, 1]
        _m_sal  = np.interp(marker_rms, _rm, d[_mask, 6])
        _m_surf = np.interp(marker_rms, _rm, d[_mask, 5])
        sc.set_array(_m_sal)

        for _key in ("circ_white", "circ_red", "circ_alt",
                    "label_9", "label_025", "label_alt"):
            if _map_artists[_key] is not None:
                _map_artists[_key].remove()
                _map_artists[_key] = None
        for _pa in _map_artists["pool_artists"]:
            _pa.remove()
        _map_artists["pool_artists"] = []

        # White solid ring: main 9-ppt toe from file
        _main_toe = float(np.interp(t_day, toe_time_f, toe_loc_f))
        if not np.isnan(_main_toe):
            _i9 = int(np.argmin(np.abs(marker_rms - _main_toe)))
            _cw = MplCircle(
                (marker_x[_i9], marker_y[_i9]),
                radius=circ_radius, transform=ax_map.transData,
                facecolor="none", edgecolor="white", linewidth=1.5, zorder=6)
            ax_map.add_patch(_cw)
            _map_artists["circ_white"] = _cw
            _map_artists["label_9"] = ax_map.annotate(
                "9 ppt",
                xy=(marker_x[_i9], marker_y[_i9]),
                xytext=(8, 4), textcoords="offset points",
                fontsize=7.5, fontweight="bold", color="white",
                path_effects=[matplotlib.patheffects.withStroke(
                    linewidth=2, foreground="black")], zorder=7)

        # Alt toe + pool circles (conditional)
        if ALT_TOE_ENABLED:
            _alt_rm, _pools = compute_connected_toe(d)
            # Orange dashed ring: alternative connected toe
            if not np.isnan(_alt_rm):
                _ia = int(np.argmin(np.abs(marker_rms - _alt_rm)))
                _ca = MplCircle(
                    (marker_x[_ia], marker_y[_ia]),
                    radius=circ_radius, transform=ax_map.transData,
                    facecolor="none", edgecolor="orange", linewidth=1.5,
                    linestyle="--", zorder=6)
                ax_map.add_patch(_ca)
                _map_artists["circ_alt"] = _ca
                _map_artists["label_alt"] = ax_map.annotate(
                    "alt toe",
                    xy=(marker_x[_ia], marker_y[_ia]),
                    xytext=(8, -12), textcoords="offset points",
                    fontsize=7, fontweight="bold", color="orange",
                    path_effects=[matplotlib.patheffects.withStroke(
                        linewidth=2, foreground="black")], zorder=7)

            # Yellow dashed rings: salinity pools
            for _pool in _pools:
                _p_center = (_pool['lo'] + _pool['hi']) / 2.0
                _p_idx = int(np.argmin(np.abs(marker_rms - _p_center)))
                _cp = MplCircle(
                    (marker_x[_p_idx], marker_y[_p_idx]),
                    radius=circ_radius, transform=ax_map.transData,
                    facecolor="none", edgecolor="yellow", linewidth=1.5,
                    linestyle="--", zorder=6)
                ax_map.add_patch(_cp)
                _map_artists["pool_artists"].append(_cp)

        # Red dashed ring: surface salinity > 0.25 ppt
        _mask_025 = _m_surf > 0.25
        if np.any(_mask_025):
            _i025 = int(np.max(np.where(_mask_025)[0]))
            _cr = MplCircle(
                (marker_x[_i025], marker_y[_i025]),
                radius=circ_radius, transform=ax_map.transData,
                facecolor="none", edgecolor="red", linewidth=2.0,
                linestyle="--", zorder=7)
            ax_map.add_patch(_cr)
            _map_artists["circ_red"] = _cr
            _map_artists["label_025"] = ax_map.annotate(
                "0.25 ppt",
                xy=(marker_x[_i025], marker_y[_i025]),
                xytext=(-65, 6), textcoords="offset points",
                fontsize=7.5, fontweight="bold", color="red",
                path_effects=[matplotlib.patheffects.withStroke(
                    linewidth=2, foreground="black")], zorder=8)

        # — Profile: update pcolormesh & boundary lines ————————————
        _msk = d[:, 1] <= MAX_RIVER_MILE
        if np.any(_msk):
            _rm_i = d[_msk, 1]
            _idx  = np.argsort(_rm_i)
            _rm_i   = _rm_i[_idx]
            _bed_i  = d[_msk, 4][_idx]
            _wsel_i = d[_msk, 2][_idx]
            _intf_i = d[_msk, 3][_idx]
            _ss_i   = d[_msk, 5][_idx]
            _sb_i   = d[_msk, 6][_idx]

            _Xn, _Yn, _Zn = interpolate_profile(
                _rm_i, _bed_i, _wsel_i, _intf_i, _ss_i, _sb_i, y_grid)
            mesh_holder[0].remove()
            mesh_holder[0] = ax_prof.pcolormesh(
                _Xn, _Yn, _Zn, cmap=cmap, norm=norm,
                shading="gouraud", zorder=1)
            line_bed.set_data(_rm_i, _bed_i)
            line_wsel.set_data(_rm_i, _wsel_i)
            line_intf.set_data(_rm_i, _intf_i)
            poly_bed_holder[0].remove()
            poly_bed_holder[0] = ax_prof.fill_between(
                _rm_i, min_elev, _bed_i,
                color="saddlebrown", alpha=1.0, zorder=2)

        # — Time series: move crimson marker ————————————————————
        vline.set_xdata([t_day, t_day])
        vspan[0].remove()
        vspan[0] = ax_ts.axvspan(
            max(t_min, t_day - span_hw),
            min(t_max, t_day + span_hw),
            color="crimson", alpha=0.20, zorder=5)

        return (sc, sup, line_bed, line_wsel, line_intf, vline)

    # ── Build and save animation ─────────────────────────────────────
    print(f"  Generating combined animation with {len(snapshots)} frames …")
    ani = animation.FuncAnimation(fig, _update_combined,
                                  frames=len(snapshots),
                                  interval=150, blit=False)

    out_base = (f"MRSWAT_Combined_{start_dt.strftime('%Y%m%d')}"
                f"_{end_dt.strftime('%Y%m%d')}")
    mp4_path = OUTPUT_DIR / f"{out_base}.mp4"
    try:
        ani.save(str(mp4_path), writer="ffmpeg", fps=FPS, dpi=150,
                 savefig_kwargs={"facecolor": "black", "edgecolor": "none"})
        print(f"  ✓ Combined animation saved: {mp4_path}")
    except Exception as e:
        print(f"  MP4 save failed ({e}); trying GIF fallback …")
        gif_path = OUTPUT_DIR / f"{out_base}.gif"
        try:
            ani.save(str(gif_path), writer="pillow", fps=FPS,
                     savefig_kwargs={"facecolor": "black", "edgecolor": "none"})
            print(f"  ✓ Combined animation saved: {gif_path}")
        except Exception as e2:
            print(f"  GIF save also failed ({e2}).")
    plt.close(fig)


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 8B — EXCEL DATA RECORDING
# ═══════════════════════════════════════════════════════════════════════════════

# ── Minimal XLSX writer (stdlib only – no openpyxl/xlsxwriter needed) ────────
#    XLSX is a ZIP archive of XML files.  We build the XML by hand so the
#    script works on any Python installation without third-party packages.

import zipfile
import xml.etree.ElementTree as ET


def _col_letter(idx):
    """Convert a 0-based column index to an Excel column letter (A, B, …)."""
    result = ""
    while True:
        result = chr(idx % 26 + ord("A")) + result
        idx = idx // 26 - 1
        if idx < 0:
            break
    return result


def _build_sheet_xml(name, rows):
    """
    Build the XML content for a single worksheet.

    Parameters
    ----------
    name : str  – sheet name (not used inside the XML, but kept for clarity)
    rows : list of list  – each inner list is one row of cell values
                           (str or numeric)

    Returns
    -------
    bytes  – UTF-8 encoded XML
    """
    WS_NS = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
    ws = ET.Element("worksheet", xmlns=WS_NS)
    sd = ET.SubElement(ws, "sheetData")
    for r_idx, row in enumerate(rows, start=1):
        row_el = ET.SubElement(sd, "row", r=str(r_idx))
        for c_idx, val in enumerate(row):
            ref = f"{_col_letter(c_idx)}{r_idx}"
            if isinstance(val, (int, float)):
                cell = ET.SubElement(row_el, "c", r=ref)
                ET.SubElement(cell, "v").text = str(val)
            else:
                cell = ET.SubElement(row_el, "c", r=ref, t="inlineStr")
                is_el = ET.SubElement(cell, "is")
                ET.SubElement(is_el, "t").text = str(val)
    return ET.tostring(ws, encoding="UTF-8", xml_declaration=True)


def _write_xlsx(filepath, sheets):
    """
    Write a minimal .xlsx workbook.

    Parameters
    ----------
    filepath : str or Path
    sheets   : list of (name, rows)  – name=sheet tab name,
               rows=list-of-lists (header row first, then data rows)
    """
    import time as _time
    filepath = str(filepath)
    CT_NS  = "http://schemas.openxmlformats.org/package/2006/content-types"
    REL_NS = "http://schemas.openxmlformats.org/package/2006/relationships"
    WB_NS  = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
    R_NS   = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"

    # Retry loop — the file may be locked if open in Excel.
    # Retries every 10 seconds for up to 5 minutes, then raises.
    _max_attempts = 30
    _retry_interval = 10   # seconds
    for _attempt in range(_max_attempts):
        try:
            _zf_handle = zipfile.ZipFile(filepath, "w", zipfile.ZIP_DEFLATED)
            break
        except PermissionError:
            if _attempt == 0:
                print(f"  ⚠ Cannot write to {filepath}")
                print(f"    File is locked (open in Excel?). "
                      f"Retrying every {_retry_interval}s "
                      f"for up to {_max_attempts * _retry_interval // 60} minutes.")
                print(f"    Please close MRSWAT_Data.xlsx to allow writing.")
            else:
                print(f"    Still locked … attempt {_attempt + 1}/{_max_attempts}")
            if _attempt == _max_attempts - 1:
                raise PermissionError(
                    f"Could not write {filepath} after "
                    f"{_max_attempts} attempts. Close the file and rerun."
                )
            _time.sleep(_retry_interval)

    with _zf_handle as zf:
        # ---- [Content_Types].xml ----
        ct = ET.Element("Types", xmlns=CT_NS)
        ET.SubElement(ct, "Default", Extension="rels",
                      ContentType="application/vnd.openxmlformats-package.relationships+xml")
        ET.SubElement(ct, "Default", Extension="xml",
                      ContentType="application/xml")
        ET.SubElement(ct, "Override", PartName="/xl/workbook.xml",
                      ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml")
        for i in range(len(sheets)):
            ET.SubElement(ct, "Override", PartName=f"/xl/worksheets/sheet{i+1}.xml",
                          ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml")
        zf.writestr("[Content_Types].xml",
                    ET.tostring(ct, encoding="unicode", xml_declaration=True))

        # ---- _rels/.rels ----
        rels = ET.Element("Relationships", xmlns=REL_NS)
        ET.SubElement(rels, "Relationship", Id="rId1",
                      Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument",
                      Target="xl/workbook.xml")
        zf.writestr("_rels/.rels",
                    ET.tostring(rels, encoding="unicode", xml_declaration=True))

        # ---- xl/_rels/workbook.xml.rels ----
        wb_rels = ET.Element("Relationships", xmlns=REL_NS)
        for i in range(len(sheets)):
            ET.SubElement(wb_rels, "Relationship", Id=f"rId{i+1}",
                          Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet",
                          Target=f"worksheets/sheet{i+1}.xml")
        zf.writestr("xl/_rels/workbook.xml.rels",
                    ET.tostring(wb_rels, encoding="unicode", xml_declaration=True))

        # ---- xl/workbook.xml ----
        wb = ET.Element("workbook", xmlns=WB_NS)
        wb.set("xmlns:r", R_NS)
        wb_sheets = ET.SubElement(wb, "sheets")
        for i, (sname, _) in enumerate(sheets):
            ET.SubElement(wb_sheets, "sheet", name=sname,
                          sheetId=str(i + 1), **{"r:id": f"rId{i+1}"})
        zf.writestr("xl/workbook.xml",
                    ET.tostring(wb, encoding="unicode", xml_declaration=True))

        # ---- xl/worksheets/sheetN.xml ----
        for i, (sname, rows) in enumerate(sheets):
            zf.writestr(f"xl/worksheets/sheet{i+1}.xml",
                        _build_sheet_xml(sname, rows))


def _excel_serial_to_datestr(serial_str):
    """
    Convert an Excel serial-date number to an MM/DD/YYYY string.

    Excel stores dates as floating-point day counts from a 1900-01-00 epoch
    (with the Lotus 1-2-3 leap-year bug for day 60).  The minimal XLSX
    writer used by this script cannot round-trip Excel date *styles*, so
    we convert detected serial dates back to readable text.

    Returns the date string, or None if the value does not look like a
    plausible serial date (covers roughly 1950-01-01 to 2199-12-31).
    """
    try:
        serial = float(serial_str)
    except (ValueError, TypeError):
        return None
    if serial < 18264 or serial > 109574:          # outside 1950–2199
        return None
    from datetime import datetime as _dt, timedelta as _td
    # Excel epoch: day 1 == 1900-01-01, but due to Lotus bug the
    # arithmetic base that works for dates after 1900-02-28 is 1899-12-30.
    base = _dt(1899, 12, 30)
    d = base + _td(days=int(serial))
    return d.strftime("%m/%d/%Y")


def _read_xlsx_sheet(filepath, sheet_index=0, sheet_name=None):
    """
    Read a single sheet from an existing .xlsx file (stdlib only).

    Handles both the minimal inline-string format produced by this script
    **and** the shared-string-table (SST) format that Excel writes when the
    user opens and re-saves the workbook.

    Parameters
    ----------
    filepath    : str or Path
    sheet_index : int  – positional index (0-based), used only when
                         *sheet_name* is None.
    sheet_name  : str or None – if provided, look up the sheet by its
                  tab name instead of by position.  This is safer when
                  the user may reorder sheets in Excel.

    Returns
    -------
    list of list  – each inner list is one row of cell values (str).
                   Returns [] if the sheet is not found.
    """
    filepath = str(filepath)
    rows_out = []
    NS = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"

    with zipfile.ZipFile(filepath, "r") as zf:
        # ── Load the Shared String Table if present ───────────────────
        sst = []                            # list of string values by index
        if "xl/sharedStrings.xml" in zf.namelist():
            sst_xml = zf.read("xl/sharedStrings.xml")
            sst_root = ET.fromstring(sst_xml)
            for si_el in sst_root.iter(f"{{{NS}}}si"):
                # Each <si> may contain a single <t> or multiple <r> runs.
                t_el = si_el.find(f"{{{NS}}}t")
                if t_el is not None and t_el.text is not None:
                    sst.append(t_el.text)
                else:
                    # Rich-text: concatenate all <r><t>…</t></r> fragments.
                    parts = []
                    for r_el in si_el.iter(f"{{{NS}}}t"):
                        if r_el.text:
                            parts.append(r_el.text)
                    sst.append("".join(parts))

        # ── Build an {rId: Target} mapping from workbook.xml.rels ─────
        rels_xml = zf.read("xl/_rels/workbook.xml.rels")
        rels_root = ET.fromstring(rels_xml)
        rid_to_target = {}
        for rel in rels_root:
            tag = rel.tag.split("}")[-1] if "}" in rel.tag else rel.tag
            if tag == "Relationship" and "worksheet" in rel.get("Type", ""):
                rid_to_target[rel.get("Id")] = rel.get("Target")

        # ── Determine which sheet XML file to read ────────────────────
        sheet_path = None
        if sheet_name is not None:
            # Parse xl/workbook.xml to match name → rId → file path
            wb_xml = zf.read("xl/workbook.xml")
            wb_root = ET.fromstring(wb_xml)
            WB_NS = NS
            R_NS  = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
            for s_el in wb_root.iter(f"{{{WB_NS}}}sheet"):
                if s_el.get("name") == sheet_name:
                    rid = s_el.get(f"{{{R_NS}}}id")
                    if rid and rid in rid_to_target:
                        sheet_path = "xl/" + rid_to_target[rid]
                    break
        else:
            # Fall back to positional index
            ordered_targets = list(rid_to_target.values())
            if sheet_index < len(ordered_targets):
                sheet_path = "xl/" + ordered_targets[sheet_index]

        if sheet_path is None:
            return rows_out

        sheet_xml = zf.read(sheet_path)
        root = ET.fromstring(sheet_xml)

        for row_el in root.iter(f"{{{NS}}}row"):
            cells = []
            for c_el in row_el:
                tag = c_el.tag.split("}")[-1] if "}" in c_el.tag else c_el.tag
                if tag != "c":
                    continue
                t = c_el.get("t", "")
                v_el = c_el.find(f"{{{NS}}}v")
                is_el = c_el.find(f"{{{NS}}}is")

                if t == "inlineStr" and is_el is not None:
                    # Inline string (produced by this script's own writer)
                    t_el = is_el.find(f"{{{NS}}}t")
                    cells.append(t_el.text if t_el is not None and t_el.text else "")
                elif t == "s" and v_el is not None and v_el.text is not None:
                    # Shared-string reference (produced by Excel on re-save)
                    idx = int(v_el.text)
                    cells.append(sst[idx] if idx < len(sst) else "")
                elif v_el is not None and v_el.text:
                    cells.append(v_el.text)
                else:
                    cells.append("")
            rows_out.append(cells)
    return rows_out


def _determine_sim_type():
    """
    Read main.inp and return 'HS' if restart flag == 1, else 'CS'.
    """
    label = "Enter a 0 for no restart conditions, enter a 1 for restart conditions"
    with open(str(MASTER_DIR / "main.inp"), "r") as fh:
        lines = fh.readlines()
    for i, line in enumerate(lines):
        if line.strip() == label and i + 1 < len(lines):
            val = lines[i + 1].strip()
            return "HS" if val == "1" else "CS"
    return "CS"  # default to cold start if not found


def record_to_excel(ref_date, sim_date):
    """
    Check for / create MRSWAT_Data.xlsx in the parent of the simulation folder
    and append a results row for this simulation.

    Sheet 'Observations':  Date | 9ppt | 0.25ppt
    Sheet 'Simulations' :  Simulation | Type | 10d 9ppt | 10d 0.25ppt |
                           28d 9ppt | 28d 0.25ppt
    """
    print("\n" + "=" * 70)
    print("  Recording forecast data to MRSWAT_Data.xlsx …")
    print("=" * 70)

    # ── Paths ────────────────────────────────────────────────────────────
    parent_dir = MASTER_DIR.parent          # master folder holding all sims
    xlsx_path  = parent_dir / "MRSWAT_Data.xlsx"
    sim_name   = MASTER_DIR.name            # e.g. MRSWAT_20260224_02
    sim_type   = _determine_sim_type()      # 'HS' or 'CS'

    # ── Load toe-location data ───────────────────────────────────────────
    if not TOE_FILE.exists():
        print(f"  ⚠ Toe output file not found: {TOE_FILE}")
        print("    (has the simulation been run yet?)")
        print("    → Writing row with N/A values.")
        iso_10d_9 = iso_10d_025 = iso_28d_9 = iso_28d_025 = None
    else:
        toe_time, toe_loc, tmdl_loc = load_toe_data()
        # toe_loc  = 9-ppt bottom isohaline (River Mile)
        # tmdl_loc = 0.25-ppt surface isohaline (River Mile)

        now = dt.datetime.now()

        def _isohaline_at_day(days_ahead):
            """Interpolate isohaline positions at `days_ahead` from now."""
            target = now + dt.timedelta(days=days_ahead)
            td = (target - ref_date).total_seconds() / 86400.0
            if td < toe_time[0] or td > toe_time[-1]:
                return None, None
            val_9   = float(np.interp(td, toe_time, toe_loc))
            val_025 = float(np.interp(td, toe_time, tmdl_loc))
            return round(val_9, 2), round(val_025, 2)

        iso_10d_9, iso_10d_025   = _isohaline_at_day(10)
        iso_28d_9, iso_28d_025   = _isohaline_at_day(28)

    print(f"  Simulation name : {sim_name}")
    print(f"  Simulation type : {sim_type}")
    print(f"  10-day 9ppt     : {iso_10d_9}")
    print(f"  10-day 0.25ppt  : {iso_10d_025}")
    print(f"  28-day 9ppt     : {iso_28d_9}")
    print(f"  28-day 0.25ppt  : {iso_28d_025}")

    # ── Sheet headers ────────────────────────────────────────────────────
    obs_header = ["Date", "9ppt", "0.25ppt"]
    sim_header = ["Simulation", "Type", "10d 9ppt", "10d 0.25ppt",
                  "28d 9ppt", "28d 0.25ppt"]

    # ── Load existing workbook or start fresh ────────────────────────────
    #    Sheets are looked up by NAME so user-reordering in Excel is safe.
    #    The Observations sheet is read back as-is and never modified;
    #    only the Simulations sheet is updated.
    if xlsx_path.exists():
        obs_rows = _read_xlsx_sheet(xlsx_path, sheet_name="Observations")
        sim_rows = _read_xlsx_sheet(xlsx_path, sheet_name="Simulations")
        # Ensure at least the header row exists
        if not obs_rows:
            obs_rows = [obs_header]
        if not sim_rows:
            sim_rows = [sim_header]
    else:
        obs_rows = [obs_header]
        sim_rows = [sim_header]

    # ── Preserve Observations sheet (user-managed) ───────────────────────
    #    Column A contains dates.  Excel stores them as serial numbers with
    #    a format style that our minimal writer cannot reproduce, so we
    #    convert any detected serial-date value back to an MM/DD/YYYY text
    #    string.  Numeric columns (B, C, …) are converted to float so
    #    Excel treats them as numbers.
    def _typed(val):
        try:
            return float(val)
        except (ValueError, TypeError):
            return val

    typed_obs_rows = []
    for r_idx, row in enumerate(obs_rows):
        if r_idx == 0:
            typed_obs_rows.append(obs_header)
            continue
        typed_row = []
        for c_idx, val in enumerate(row):
            if c_idx == 0:
                # Column A = date: convert serial dates to readable text
                ds = _excel_serial_to_datestr(val)
                typed_row.append(ds if ds is not None else val)
            else:
                typed_row.append(_typed(val))
        typed_obs_rows.append(typed_row)

    # ── Convert Simulations rows to proper types for numeric cells ───────
    #    Columns A (name) and B (type) are always strings; C–F are numeric.
    typed_sim_rows = []
    for r_idx, row in enumerate(sim_rows):
        if r_idx == 0:
            typed_sim_rows.append(sim_header)
            continue
        typed_row = []
        for c_idx, val in enumerate(row):
            if c_idx < 2:
                # Simulation name and type — always keep as text
                typed_row.append(str(val) if val else "")
            else:
                typed_row.append(_typed(val))
        typed_sim_rows.append(typed_row)

    # ── Write results row to Simulations sheet ───────────────────────────
    new_row = [
        sim_name,
        sim_type,
        iso_10d_9   if iso_10d_9   is not None else "N/A",
        iso_10d_025 if iso_10d_025 is not None else "N/A",
        iso_28d_9   if iso_28d_9   is not None else "N/A",
        iso_28d_025 if iso_28d_025 is not None else "N/A",
    ]

    # Find top-most empty row after the header (row 0); "empty" means all
    # cells are blank strings.  If none found, append at the end.
    empty_idx = None
    for i, row in enumerate(typed_sim_rows):
        if i == 0:                          # skip header
            continue
        if all(str(v).strip() == "" for v in row):
            empty_idx = i
            break

    if empty_idx is not None:
        typed_sim_rows[empty_idx] = new_row
    else:
        typed_sim_rows.append(new_row)

    # ── Write workbook ───────────────────────────────────────────────────
    #    Observations is written first (preserved as read), then Simulations.
    _write_xlsx(xlsx_path, [
        ("Observations", typed_obs_rows),
        ("Simulations",  typed_sim_rows),
    ])
    print(f"  ✓ MRSWAT_Data.xlsx updated (Simulations sheet): {xlsx_path}")


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 9 — MAIN DRIVER
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    print("\n" + "=" * 70)
    print("  MRSWAT Visualization Suite")
    print(f"  Master folder : {MASTER_DIR}")
    print(f"  Output folder : {OUTPUT_DIR}")
    print(f"  Switches:")
    print(f"    Summary Figure     : {'ON' if MAKE_SUMMARY_FIGURE     else 'OFF'}")
    print(f"    Salinity Map       : {'ON' if MAKE_SALINITY_MAP       else 'OFF'}")
    print(f"    Map Animation      : {'ON' if MAKE_MAP_ANIMATION      else 'OFF'}")
    print(f"    Animation          : {'ON' if MAKE_ANIMATION          else 'OFF'}")
    print(f"    Calibration Anim.  : {'ON' if MAKE_CALIBRATION_ANIMATION else 'OFF'}")
    print(f"    Combined Animation : {'ON' if MAKE_COMBINED_ANIMATION else 'OFF'}")
    print(f"    Record to Excel    : {'ON' if RECORD_TO_EXCEL        else 'OFF'}")
    print(f"    Alt Toe Enabled    : {'ON' if ALT_TOE_ENABLED       else 'OFF'}")
    if ALT_TOE_ENABLED:
        print(f"      Salinity         : {ALT_TOE_SALINITY} ppt")
        print(f"      Min Thickness    : {ALT_TOE_MIN_THICKNESS} ft")
        print(f"      Max Gap          : {ALT_TOE_MAX_GAP} RM")
        print(f"      Min Pool Depth   : {ALT_TOE_MIN_POOL_DEPTH} ft")
    print("=" * 70)

    ref_date, sim_date = load_reference_and_sim_dates()
    print(f"  Reference date : {ref_date.strftime('%Y-%m-%d')}")
    print(f"  Simulation date: {sim_date}")

    if MAKE_SUMMARY_FIGURE:
        make_summary_figure(ref_date, sim_date)

    if MAKE_SALINITY_MAP:
        make_salinity_map(ref_date)

    if MAKE_MAP_ANIMATION:
        make_salinity_map_animation(ref_date)

    if MAKE_ANIMATION:
        make_animation(ref_date)

    if MAKE_CALIBRATION_ANIMATION:
        make_calibration_animation(ref_date)

    if MAKE_COMBINED_ANIMATION:
        make_combined_animation(ref_date)

    if RECORD_TO_EXCEL:
        record_to_excel(ref_date, sim_date)

    print("\n" + "=" * 70)
    print("  All requested outputs complete.")
    print("=" * 70)


if __name__ == "__main__":
    main()
