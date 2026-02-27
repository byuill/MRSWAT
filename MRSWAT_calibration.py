"""
MRSWAT Calibration Script
=========================
Compare modeled salt-wedge toe location against vessel-based CTD field
observations and produce a two-panel calibration figure:

  Top panel    – Time series of modeled toe (red line) with observed toe
                 locations shown as open circles with vertical error bars.
  Bottom panel – Quantile–quantile (Q-Q) plot of observed vs. modeled toe.

Data sources
------------
  Modeled : ../mr-toe-location-output.txt   (elapsed days → River Mile)
  Observed: \\\\mvd\\mvn\\Data_EDR\\02 Stream Gaging Section\\...\\Field Data
            CTD cast CSV files organised by year → SWW_YYYYMMDD folders.

River-mile conversion uses the USACE waypoint table in river_miles.csv.
"""

# ═══════════════════════════════════════════════════════════════════════════════
# IMPORTS
# ═══════════════════════════════════════════════════════════════════════════════
from pathlib import Path
import datetime as dt
import csv
import re
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import dates as mdates

# ═══════════════════════════════════════════════════════════════════════════════
# PATHS
# ═══════════════════════════════════════════════════════════════════════════════
SCRIPT_DIR   = Path(__file__).resolve().parent           # …/Code Storage/
SIM_DIR      = SCRIPT_DIR.parent                         # …/MRSWAT_YYYYMMDD_NN/
TOE_FILE     = SIM_DIR / "mr-toe-location-output.txt"
DATE_FILE    = SIM_DIR / "Date_Information.txt"
RM_CSV       = SCRIPT_DIR / "river_miles.csv"
OUTPUT_DIR   = SIM_DIR / "Output"
OUTPUT_DIR.mkdir(exist_ok=True)

FIELD_DATA_ROOT = Path(
    r"\\mvd\mvn\Data_EDR\02 Stream Gaging Section"
    r"\Gaging Stations\Gages_RandomProjects\Saltwater Wedge\Field Data"
)

# Salinity threshold for toe identification (ppt)
TOE_SALINITY = 9.0


# ═══════════════════════════════════════════════════════════════════════════════
# UTILITY — RIVER-MILE WAYPOINTS & LAT/LON → RM CONVERSION
# ═══════════════════════════════════════════════════════════════════════════════

def load_river_mile_waypoints():
    """
    Load USACE river-mile marker coordinates from the companion CSV.

    Returns
    -------
    np.ndarray, shape (K, 3) – columns: model_RM, latitude, longitude
    """
    mi_rows, sw_rows = [], []
    with open(str(RM_CSV), "r") as cf:
        reader = csv.DictReader(cf)
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


def _haversine_ft(lat1, lon1, lat2, lon2):
    """Great-circle distance between two points in feet."""
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
    then refines by projecting onto the two adjacent river-mile segments.

    Parameters
    ----------
    lat, lon : float
    waypoints : np.ndarray, shape (K, 3)

    Returns
    -------
    float – estimated river mile
    """
    dists = np.array([
        _haversine_ft(lat, lon, wp[1], wp[2]) for wp in waypoints
    ])
    i_best = int(np.argmin(dists))

    # Try projecting onto the segment before and after the nearest waypoint
    best_rm = waypoints[i_best, 0]
    for i_start in [max(0, i_best - 1), i_best]:
        i_end = i_start + 1
        if i_end >= len(waypoints):
            continue
        ax, ay = waypoints[i_start, 1], waypoints[i_start, 2]
        bx, by = waypoints[i_end, 1],   waypoints[i_end, 2]
        # Project (lat, lon) onto line segment (a -> b) in lat/lon space
        abx, aby = bx - ax, by - ay
        apx, apy = lat - ax, lon - ay
        ab2 = abx * abx + aby * aby
        if ab2 < 1e-14:
            continue
        t = max(0.0, min(1.0, (apx * abx + apy * aby) / ab2))
        rm_a = waypoints[i_start, 0]
        rm_b = waypoints[i_end, 0]
        rm_proj = rm_a + t * (rm_b - rm_a)
        # Distance from point to projected location on segment
        plat = ax + t * abx
        plon = ay + t * aby
        d = _haversine_ft(lat, lon, plat, plon)
        if d < dists[i_best] + 100:  # allow small tolerance
            best_rm = rm_proj
            break

    return best_rm


# ═══════════════════════════════════════════════════════════════════════════════
# LOAD MODELED DATA
# ═══════════════════════════════════════════════════════════════════════════════

def load_reference_date():
    """Read reference date from Date_Information.txt."""
    with open(str(DATE_FILE), "r") as f:
        lines = f.readlines()
    header = lines[0].strip().split()
    values = lines[1].strip().split()
    ref_str = values[header.index("Ref_Date")]
    return dt.datetime.strptime(ref_str, "%Y%m%d")


def load_modeled_toe():
    """
    Read mr-toe-location-output.txt.

    Returns
    -------
    dates : list of datetime
    toe_rm : np.ndarray
    """
    ref_date = load_reference_date()
    times, toes = [], []
    with open(str(TOE_FILE), "r") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 2:
                continue
            try:
                t = float(parts[0])
                r = float(parts[1])
            except ValueError:
                continue
            times.append(t)
            toes.append(r)
    times = np.array(times)
    toes  = np.array(toes)
    dates = [ref_date + dt.timedelta(days=float(t)) for t in times]
    return dates, times, toes, ref_date


def load_two_col_ts(path):
    """Load a simple two-column time series file (TIME value) into arrays.

    Expects a header row with names and two whitespace- or tab-separated
    columns. Returns (times, values) where `times` are floats (days).
    """
    times, vals = [], []
    p = Path(path)
    if not p.exists():
        return np.array([]), np.array([])
    with open(str(p), "r") as f:
        # skip header
        header = f.readline()
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = re.split(r"\s+|,", line)
            try:
                t = float(parts[0])
                v = float(parts[1])
            except Exception:
                continue
            times.append(t)
            vals.append(v)
    return np.array(times), np.array(vals)


# ═══════════════════════════════════════════════════════════════════════════════
# LOAD OBSERVED DATA
# ═══════════════════════════════════════════════════════════════════════════════

def _parse_ctd_csv(filepath):
    """
    Parse a single CTD cast CSV file.

    Returns
    -------
    dict with keys:
        'lat'      : float – start latitude
        'lon'      : float – start longitude
        'max_sal'  : float – maximum salinity observed in the cast (ppt)
    or None if the file cannot be parsed.
    """
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

    # Data matrix starts at row 30 (0-indexed line 29).
    # Row 29 (0-indexed 28) is the column header.
    # Salinity is column index 5.
    max_sal = 0.0
    for line in lines[29:]:
        parts = line.strip().split(",")
        if len(parts) < 6:
            continue
        try:
            sal = float(parts[5])
            if sal > max_sal:
                max_sal = sal
        except (ValueError, IndexError):
            continue

    return {"lat": lat, "lon": lon, "max_sal": max_sal}


def load_observed_toe(date_start, date_end, waypoints):
    """
    Scan the field-data directory tree for CTD measurement days that fall
    within [date_start, date_end].  For each measurement day, determine
    the observed toe position.

    The observed toe is defined as the upstream-most (highest RM) cast
    location where at least one depth sample has salinity ≥ 9 ppt.

    Error bars capture the spatial discretisation uncertainty:
      • lower bound = RM of the toe cast itself (highest saline cast)
      • upper bound = RM of the next cast upstream (first fresh cast),
                      or the toe RM if no fresh cast exists upstream.

    Returns
    -------
    obs_dates      : list of datetime.date
    obs_toe_rm     : np.ndarray – estimated toe RM for each day
    obs_toe_lo     : np.ndarray – lower error-bar extent  (RM)
    obs_toe_hi     : np.ndarray – upper error-bar extent  (RM)
    """
    obs_dates  = []
    obs_toe_rm = []
    obs_toe_lo = []
    obs_toe_hi = []

    if not FIELD_DATA_ROOT.exists():
        print(f"  WARNING: Field data root not found: {FIELD_DATA_ROOT}")
        return obs_dates, np.array([]), np.array([]), np.array([])

    # Iterate year folders
    for year_dir in sorted(FIELD_DATA_ROOT.iterdir()):
        if not year_dir.is_dir():
            continue
        # Only consider years that could overlap with the simulation period
        try:
            year_val = int(year_dir.name)
        except ValueError:
            continue
        if year_val < date_start.year or year_val > date_end.year:
            continue

        # Iterate SWW_YYYYMMDD folders
        for sww_dir in sorted(year_dir.iterdir()):
            if not sww_dir.is_dir():
                continue
            m = re.match(r"SWW_(\d{8})", sww_dir.name)
            if not m:
                continue
            try:
                obs_date = dt.datetime.strptime(m.group(1), "%Y%m%d").date()
            except ValueError:
                continue

            if obs_date < date_start.date() or obs_date > date_end.date():
                continue

            # Parse all individual cast CSVs in this folder
            casts = []
            for csv_file in sorted(sww_dir.glob("*.csv")):
                if csv_file.name.lower() == "summary.csv":
                    continue
                result = _parse_ctd_csv(csv_file)
                if result is not None:
                    result["rm"] = latlon_to_rm(
                        result["lat"], result["lon"], waypoints)
                    casts.append(result)

            if not casts:
                continue

            # Sort casts by river mile (downstream → upstream)
            casts.sort(key=lambda c: c["rm"])

            # Identify the upstream-most cast with max salinity ≥ 9 ppt
            saline_casts = [c for c in casts if c["max_sal"] >= TOE_SALINITY]
            if not saline_casts:
                continue  # no salinity ≥ 9 ppt observed this day

            # Toe = upstream-most saline cast
            toe_cast = max(saline_casts, key=lambda c: c["rm"])
            toe_rm   = toe_cast["rm"]

            # Error bounds:
            #   lower = toe cast RM (can't be farther downstream than this)
            #   upper = next upstream cast that is FRESH, giving the bracket
            fresh_upstream = [c for c in casts
                             if c["rm"] > toe_rm and c["max_sal"] < TOE_SALINITY]
            if fresh_upstream:
                # The true toe is between the last saline cast and the first
                # fresh cast upstream
                first_fresh = min(fresh_upstream, key=lambda c: c["rm"])
                upper_rm = first_fresh["rm"]
            else:
                # All casts upstream of the toe are also saline, or there
                # are no casts upstream — can't bound from above
                upper_rm = toe_rm

            obs_dates.append(obs_date)
            obs_toe_rm.append(toe_rm)
            obs_toe_lo.append(toe_rm)
            obs_toe_hi.append(upper_rm)

            n_casts = len(casts)
            n_saline = len(saline_casts)
            print(f"    {obs_date}  –  {n_casts} casts, {n_saline} saline, "
                  f"toe ≈ RM {toe_rm:.1f}  "
                  f"[{toe_rm:.1f} – {upper_rm:.1f}]")

    return (obs_dates,
            np.array(obs_toe_rm),
            np.array(obs_toe_lo),
            np.array(obs_toe_hi))


# ═══════════════════════════════════════════════════════════════════════════════
# PLOTTING
# ═══════════════════════════════════════════════════════════════════════════════

def make_calibration_figure(mod_dates, mod_toe, ref_date,
                            obs_dates, obs_toe, obs_lo, obs_hi):
    """
    Three-panel calibration figure:

    Top:    time series — modeled toe (red line) + observed toe (open
            circles with vertical error bars).
    Mid:    discharge & stage time series from Discharges.txt and Stages.txt.
    Bottom: Q-Q plot — observed (x) vs. modeled (y).  Legend and stats
            box placed to the right to avoid overlap.
    """
    fig = plt.figure(figsize=(13, 10))
    gs = fig.add_gridspec(3, 1, height_ratios=[2, 1, 1], hspace=0.28)
    ax1 = fig.add_subplot(gs[0, 0])
    ax_mid = fig.add_subplot(gs[1, 0])
    ax3 = fig.add_subplot(gs[2, 0])

    # ── Top panel: time-series comparison ────────────────────────────────
    ax1.plot(mod_dates, mod_toe, color="red", lw=1.2, label="Modeled toe (9 ppt)")

    if len(obs_dates) > 0:
        obs_datetimes = [dt.datetime.combine(d, dt.time()) for d in obs_dates]
        yerr_lo = obs_toe - obs_lo
        yerr_hi = obs_hi - obs_toe
        ax1.errorbar(obs_datetimes, obs_toe,
                     yerr=[yerr_lo, yerr_hi],
                     fmt="o", mfc="none", mec="blue", ecolor="blue",
                     capsize=4, ms=6, lw=1.0,
                     label="Observed toe (CTD casts)")

    # Top panel: don't show calendar x-axis labels (middle panel will show them)
    ax1.set_xlabel("")
    ax1.xaxis.set_major_formatter(matplotlib.ticker.NullFormatter())
    ax1.set_ylabel("Toe Location (River Mile)")
    ax1.set_title("MRSWAT Calibration — Modeled vs. Observed Salt-Wedge Toe",
                  fontsize=13, fontweight="bold")
    ax1.legend(loc="best")
    ax1.grid(True, alpha=0.3)

    # ── Middle panel: Discharge & Stage time series ─────────────────────
    disc_t, disc_v = load_two_col_ts(SIM_DIR / "Discharges.txt")
    stage_t, stage_v = load_two_col_ts(SIM_DIR / "Stages.txt")
    if disc_t.size > 0:
        disc_dates = [ref_date + dt.timedelta(days=float(t)) for t in disc_t]
        ax_mid.plot(disc_dates, disc_v, color="steelblue", lw=1.2, label="Discharge (cfs)")
        ax_mid.set_ylabel("Discharge (cfs)", color="steelblue")
        ax_mid.tick_params(axis="y", labelcolor="steelblue")
    if stage_t.size > 0:
        stage_dates = [ref_date + dt.timedelta(days=float(t)) for t in stage_t]
        ax_stage = ax_mid.twinx()
        ax_stage.plot(stage_dates, stage_v, color="darkorange", lw=1.0, label="Stage (ft)")
        ax_stage.set_ylabel("Stage (ft)", color="darkorange")
        ax_stage.tick_params(axis="y", labelcolor="darkorange")

    # Middle panel: show calendar date ticks
    ax_mid.set_xlabel("Date")
    locator = mdates.AutoDateLocator()
    formatter = mdates.ConciseDateFormatter(locator)
    ax_mid.xaxis.set_major_locator(locator)
    ax_mid.xaxis.set_major_formatter(formatter)
    ax_mid.tick_params(axis="x", rotation=30, labelsize=8)
    ax_mid.set_title("Hydrodynamic Forcing — Discharge & Stage", fontsize=11, fontweight="bold")
    ax_mid.grid(True, alpha=0.2)

    # ── Bottom panel: Q-Q plot ───────────────────────────────────────────
    if len(obs_dates) > 0:
        mod_times_num = np.array([(d - ref_date).total_seconds() / 86400.0 for d in mod_dates])
        obs_times_num = np.array([(dt.datetime.combine(d, dt.time()) - ref_date).total_seconds() / 86400.0 for d in obs_dates])
        mod_at_obs = np.interp(obs_times_num, mod_times_num, mod_toe)

        sc = ax3.scatter(obs_toe, mod_at_obs, edgecolors="blue", facecolors="none",
                         s=50, zorder=3, label="Observation days")

        all_vals = np.concatenate([obs_toe, mod_at_obs])
        lo, hi = np.floor(all_vals.min()) - 2, np.ceil(all_vals.max()) + 2
        ax3.plot([lo, hi], [lo, hi], "k--", lw=0.8, label="1:1 line")
        ax3.set_xlim(lo, hi)
        ax3.set_ylim(lo, hi)
        ax3.set_aspect("equal", adjustable="box")

        # Stats annotation -> place to the right of axes
        residuals = mod_at_obs - obs_toe
        rmse = np.sqrt(np.mean(residuals ** 2))
        bias = np.mean(residuals)
        r = np.corrcoef(obs_toe, mod_at_obs)[0, 1] if len(obs_toe) > 1 else np.nan
        stats_text = (f"n = {len(obs_toe)}\n"
                      f"RMSE = {rmse:.2f} RM\n"
                      f"Bias = {bias:+.2f} RM\n"
                      f"r = {r:.3f}")
        ax3.text(1.02, 0.95, stats_text, transform=ax3.transAxes,
                 fontsize=9, verticalalignment="top", family="monospace",
                 bbox=dict(boxstyle="round,pad=0.4", facecolor="#f0f0f0",
                           edgecolor="gray", lw=0.8), ha="left")

        # Keep the stats box to the right; remove the Q-Q legend to reduce clutter
        pass
    else:
        ax3.text(0.5, 0.5, "No matching observations found",
                 transform=ax3.transAxes, ha="center", va="center",
                 fontsize=12, color="gray")

    ax3.set_xlabel("Observed Toe (River Mile)")
    ax3.set_ylabel("Modeled Toe (River Mile)")
    ax3.set_title("Quantile–Quantile Plot", fontsize=11, fontweight="bold")
    ax3.grid(True, alpha=0.3)

    # ── Save ─────────────────────────────────────────────────────────────
    fig_path = OUTPUT_DIR / "MRSWAT_Calibration.jpg"
    plt.savefig(str(fig_path), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"\n  ✓ Calibration figure saved: {fig_path}")


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    print("\n" + "=" * 70)
    print("  MRSWAT Calibration — Modeled vs. Observed Toe Location")
    print("=" * 70)

    # ── Load modeled data ────────────────────────────────────────────────
    mod_dates, mod_times, mod_toe, ref_date = load_modeled_toe()
    print(f"  Reference date   : {ref_date.strftime('%Y-%m-%d')}")
    print(f"  Modeled period   : {mod_dates[0].strftime('%Y-%m-%d')} → "
          f"{mod_dates[-1].strftime('%Y-%m-%d')}")
    print(f"  Modeled snapshots: {len(mod_dates)}")

    # ── Load river-mile waypoints ────────────────────────────────────────
    waypoints = load_river_mile_waypoints()
    print(f"  Waypoints loaded : {len(waypoints)} river-mile markers")

    # ── Scan for matching field observations ─────────────────────────────
    print(f"\n  Scanning field data: {FIELD_DATA_ROOT}")
    obs_dates, obs_toe, obs_lo, obs_hi = load_observed_toe(
        mod_dates[0], mod_dates[-1], waypoints)

    if len(obs_dates) == 0:
        print("\n  ⚠ No field observations found within the modeled period.")
        print("    The figure will show the modeled time series only.")
    else:
        print(f"\n  Observation days found: {len(obs_dates)}")

    # ── Generate calibration figure ──────────────────────────────────────
    make_calibration_figure(mod_dates, mod_toe, ref_date,
                            obs_dates, obs_toe, obs_lo, obs_hi)

    print("\n" + "=" * 70)
    print("  Calibration complete.")
    print("=" * 70)


if __name__ == "__main__":
    main()
