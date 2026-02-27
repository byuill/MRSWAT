#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
=============================================================================
 MRSWAT2000.py  –  Mississippi River Salt-Water Intrusion Pre-Processor
=============================================================================
 Purpose
 -------
   Automatically retrieves, processes, and exports the boundary-condition
   time-series data that drive the MRSWAT2000 salt-water intrusion model.

   The script performs five main tasks:
     1. Download observed **water levels** (tide gauge) from the NOAA CO-OPS
        API and predicted tidal water levels from the noaa_coops library.
     2. Download observed **instantaneous discharge** from the USGS NWIS API.
     3. Optionally develop a **stage-discharge rating curve** using USGS daily
        data and polynomial regression.
     4. Download **forecasted stage/flow** from the National Water Prediction
        Service (NWPS) API, and convert forecast stages to discharge via the
        rating curve.
     5. Combine observed + forecast data, produce diagnostic plots (saved as
        PNG), and export tab-delimited boundary-condition files for the model.

 Data Sources
 ------------
   • NOAA CO-OPS Tides & Currents API  (water levels)
   • USGS NWIS Instantaneous-Value API (discharge)
   • USGS NWIS Daily-Value API          (stage & flow for rating curve)
   • NOAA NWPS (api.water.noaa.gov)     (forecast stage/flow)

 Outputs  (saved to a date-stamped simulation folder)
 -------
   • Water_Level_BC_<folder>.png   – combined hourly stage plot
   • Discharge_BC_<folder>.png     – combined daily discharge plot
   • Stages.txt                    – elapsed-time ↔ stage table  (tab-delimited)
   • Discharges.txt                – elapsed-time ↔ discharge table (tab-delimited)

 Author / History
 ----------------
   Assembled from multiple standalone scripts and reorganised into a single
   pre-processing pipeline.
=============================================================================
"""

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 0 — IMPORTS
# ═══════════════════════════════════════════════════════════════════════════════
# Standard library
from datetime import datetime, timedelta, timezone
import io
import time

# Third-party – data & numerics
import numpy as np
import pandas as pd
import requests

# Third-party – plotting
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.ticker as ticker

# Third-party – NOAA helper
from noaa_coops import Station

# File-system
from pathlib import Path
import shutil


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 1 — USER-CONFIGURABLE PARAMETERS
# ═══════════════════════════════════════════════════════════════════════════════
# These constants control which stations are queried, the look-back window,
# datum / unit conventions, and whether a rating curve is regenerated.

# -- Flags ----------------------------------------------------------------
generate_RC = 1                # Set to 1 to build a rating curve from USGS
                               # daily data; 0 to skip that step entirely.

# -- Station identifiers --------------------------------------------------
STATION_ID       = "8760922"   # NOAA CO-OPS: Pilots Station East, SW Pass, LA
USGS_STATION_ID  = "07374000"  # USGS NWIS:   Mississippi River at Baton Rouge
NWPS_STATION_ID  = "BTRL1"    # NWPS gauge:  Baton Rouge (Stage)

# -- Datum / units / time zone --------------------------------------------
DATUM     = "MLLW"             # Options: MLLW, NAVD, MSL
UNITS     = "english"          # Options: english, metric
TIME_ZONE = "lst_ldt"          # Options: gmt, lst_ldt (local standard/daylight)

# -- Temporal parameters ---------------------------------------------------
Days_Back   = 180              # Number of past days to retrieve for observed data
MLLW2NAVD88 = -0.68           # Additive offset (ft) to convert MLLW → NAVD88

# -- Forecast Calibration --------------------------------------------------
Discharge_Forecast_Calibration = 0.9
    # Multiplicative scaling factor for forecast discharge (valid range: 0.5–2.0).
    # The scaling strength ramps with temporal distance from the current day:
    #   • Day +1  →  10 % of the factor applied
    #   • Day +7+ → 100 % of the factor applied  (linear ramp in between)
    # A value of 1.0 means no scaling.

Stage_Forecast_Calibration = "STM"
    # Controls how the forecasted stage time-series is adjusted:
    #   True    → shift the forecast so its temporal mean matches the
    #             observed-stage temporal mean (bias correction).
    #   False   → no adjustment applied to forecast stage.
    #   "STM"  → Short Term Match: compute the median of the last 24 h of
    #             observed stage and the median of the first 24 h of forecast
    #             stage; apply the difference as an additive offset to the
    #             entire forecast so that those two medians align.
    #   <float> → additive offset (ft) applied to the entire forecast
    #             stage time-series (e.g. 0.25 adds 0.25 ft).


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 2 — BOOKKEEPING  (resolve output folder)
# ═══════════════════════════════════════════════════════════════════════════════
# Outputs are written directly into the folder one level above the directory
# that contains this script (i.e. the parent of "Code Storage/").

today_str         = datetime.now().strftime('%Y%m%d')
reference_date    = datetime.now() - timedelta(days=Days_Back)  # used later for elapsed-time export
forecast_end_date = datetime.now() + timedelta(days=30)         # target end of all forecast data
final_date_str    = forecast_end_date.strftime('%Y%m%d')

# Output directory: one level above the directory containing this script
folder = Path(__file__).resolve().parent.parent

print("\n" + "="*80)
print(f"  MRSWAT2000 Pre-Processor – {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("="*80)
print(f"  ✓ Output folder: {folder}\n")

# -----------------------------------------------------------------------------
# Write Date_Information.txt into the simulation folder
# Columns: Ref_Date<TAB>Sim_Date<TAB>Final_Date
# Row 2   : dates as YYYYMMDD
# Row 3   : dates as elapsed days after reference date
# -----------------------------------------------------------------------------
_ref_elapsed   = 0
_sim_elapsed   = Days_Back
_final_elapsed = Days_Back + 30

try:
    date_info_path = folder / 'Date_Information.txt'
    with open(date_info_path, 'w') as f:
        f.write('Ref_Date\tSim_Date\tFinal_Date\n')
        f.write(f"{reference_date.strftime('%Y%m%d')}\t{today_str}\t{final_date_str}\n")
        f.write(f"{_ref_elapsed}\t{_sim_elapsed}\t{_final_elapsed}\n")
    print(f"  ✓ Wrote Date_Information.txt: {date_info_path.name}")
    print(f"    Ref_Date={reference_date.strftime('%Y%m%d')}  "
          f"Sim_Date={today_str}  Final_Date={final_date_str}")
    print(f"    Elapsed days from ref: Ref=0  Sim={_sim_elapsed}  Final={_final_elapsed}\n")
except Exception as e:
    print(f"  ⚠ Failed to write Date_Information.txt: {e}\n")


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 3 — FUNCTION DEFINITIONS
# ═══════════════════════════════════════════════════════════════════════════════
# All data-retrieval logic is encapsulated in functions so that each step can
# be tested or reused independently.
# ---------------------------------------------------------------------------


def get_noaa_data(days_back=Days_Back):
    """
    Fetch **observed** water levels from the NOAA CO-OPS API.

    The API imposes a 31-day maximum per request, so this function splits the
    requested window into ≤ 30-day chunks, concatenates the results, and
    returns a single DataFrame indexed by 'Date Time'.

    Parameters
    ----------
    days_back : int
        Number of past days to retrieve (default: ``Days_Back``).

    Returns
    -------
    pd.DataFrame or None
        DataFrame with a ``'Water Level'`` column (plus quality flags).
        Returns ``None`` if no data could be retrieved.
    """
    base_url = "https://api.tidesandcurrents.noaa.gov/api/prod/datagetter"

    # --- define the full time window (UTC) --------------------------------
    end_date   = datetime.now(timezone.utc)
    start_date = end_date - timedelta(days=days_back)

    all_data      = []
    current_start = start_date

    # --- loop in ≤ 30-day chunks ------------------------------------------
    while current_start < end_date:
        current_end = min(current_start + timedelta(days=30), end_date)

        begin_str = current_start.strftime("%Y%m%d")
        end_str   = current_end.strftime("%Y%m%d")
        print(f"Fetching data: {begin_str} to {end_str}...")

        params = {
            "begin_date":  begin_str,
            "end_date":    end_str,
            "station":     STATION_ID,
            "product":     "water_level",
            "datum":       DATUM,
            "units":       UNITS,
            "time_zone":   TIME_ZONE,
            "format":      "csv",
            "application": "Python_Script",
        }

        try:
            response = requests.get(base_url, params=params)
            if response.status_code == 200 and "Error" not in response.text:
                chunk_df = pd.read_csv(io.StringIO(response.text))
                chunk_df.columns = chunk_df.columns.str.strip()   # prevent KeyErrors
                all_data.append(chunk_df)
            else:
                print(f"  Warning: No data or error for {begin_str}-{end_str}")
        except Exception as e:
            print(f"  Error: {e}")

        time.sleep(0.5)                                          # polite rate-limit
        current_start = current_end + timedelta(days=1)

    # --- concatenate all chunks -------------------------------------------
    if all_data:
        final_df = pd.concat(all_data, ignore_index=True)
        final_df['Date Time'] = pd.to_datetime(final_df['Date Time'])
        final_df.set_index('Date Time', inplace=True)
        return final_df
    else:
        return None


# ---------------------------------------------------------------------------

def get_usgs_instantaneous_discharge(site_id, days_back):
    """
    Fetch **instantaneous** discharge (parameter 00060) from the USGS NWIS
    Instantaneous-Value ('iv') service.

    Parameters
    ----------
    site_id   : str   – USGS station number (e.g. ``"07374000"``).
    days_back : int   – Number of past days to retrieve.

    Returns
    -------
    pd.DataFrame
        DataFrame with column ``'Discharge_cfs'`` indexed by ``'Date'``
        (timezone-naive).  Empty DataFrame on failure.
    """
    url = "https://waterservices.usgs.gov/nwis/iv/?format=json"

    params = {
        "sites":       site_id,
        "period":      f"P{days_back}D",   # ISO-8601 duration
        "parameterCd": "00060",            # Discharge
        "siteStatus":  "all",
    }

    try:
        print(f"Fetching USGS Instantaneous data for site {site_id}...")
        response = requests.get(url, params=params)
        response.raise_for_status()
        data = response.json()

        ts_data = []
        if ('value' in data
                and 'timeSeries' in data['value']
                and len(data['value']['timeSeries']) > 0):
            values = data['value']['timeSeries'][0]['values'][0]['value']
            for item in values:
                ts_data.append({
                    'Date':          pd.to_datetime(item['dateTime'], utc=True),
                    'Discharge_cfs': float(item['value']),
                })

            df = pd.DataFrame(ts_data)
            if not df.empty:
                df.set_index('Date', inplace=True)
                df.index = df.index.tz_convert(None)   # drop TZ for simpler downstream use
                return df
        else:
            print("USGS Response contained no timeSeries data.")

        return pd.DataFrame()

    except Exception as e:
        print(f"Error fetching USGS data: {e}")
        return pd.DataFrame()


# ---------------------------------------------------------------------------

def get_usgs_daily_data(site_id, days_back, param_cd):
    """
    Fetch **daily-mean** data for a single parameter from the USGS NWIS
    Daily-Value ('dv') service.

    Used primarily for building the rating curve (stage vs. discharge).

    Parameters
    ----------
    site_id   : str   – USGS station number.
    days_back : int   – Number of past days.
    param_cd  : str   – USGS parameter code:
                          ``'00060'`` = Discharge,
                          ``'00065'`` = Gage Height.

    Returns
    -------
    pd.DataFrame
        Single-column ``'Value'`` indexed by ``'Date'``.
        Empty DataFrame on failure.
    """
    url = "https://waterservices.usgs.gov/nwis/dv/?format=json"

    end_date   = datetime.now()
    start_date = end_date - timedelta(days=days_back)

    params = {
        "sites":       site_id,
        "startDT":     start_date.strftime("%Y-%m-%d"),
        "endDT":       end_date.strftime("%Y-%m-%d"),
        "parameterCd": param_cd,
        "statCd":      "00003",          # Mean value
    }

    try:
        print(f"Fetching parameter {param_cd} for site {site_id}...")
        response = requests.get(url, params=params)
        response.raise_for_status()
        data = response.json()

        ts_data = []
        if 'value' in data and 'timeSeries' in data['value']:
            for ts in data['value']['timeSeries']:
                for item in ts['values'][0]['value']:
                    ts_data.append({
                        'Date':  pd.to_datetime(item['dateTime']),
                        'Value': float(item['value']),
                    })

        df = pd.DataFrame(ts_data)
        if not df.empty:
            df.set_index('Date', inplace=True)
        return df

    except Exception as e:
        print(f"Error fetching param {param_cd}: {e}")
        return pd.DataFrame()


# ---------------------------------------------------------------------------

def qaqc_timeseries(series, threshold_multiplier=2.5, label="series"):
    """
    Quality-assure a time series by detecting and replacing erroneous spikes.

    Algorithm
    ---------
    1. Compute the absolute value of the step-to-step differences.
    2. Calculate the mean of those absolute differences (mean_abs_delta).
    3. Flag any value whose absolute change from the previous valid value
       exceeds ``threshold_multiplier × mean_abs_delta`` as erroneous.
    4. Set flagged (and any pre-existing NaN) values to NaN.
    5. Linearly interpolate across the gaps using the surrounding valid values.

    Parameters
    ----------
    series              : pd.Series  – time-indexed series to clean.
    threshold_multiplier: float      – spike-detection multiplier (default 2.5).
    label               : str        – descriptive name used in printed messages.

    Returns
    -------
    tuple(pd.Series, pd.Series)
        ``(cleaned, replaced_mask)`` where ``cleaned`` is the QA'd series and
        ``replaced_mask`` is a boolean Series that is ``True`` at every index
        whose value was either a pre-existing NaN or a detected spike (i.e.
        every position that was filled by linear interpolation).
    """
    s = series.copy().astype(float)

    # Record which positions were already NaN before any processing
    original_nan_mask = s.isna()
    n_original_nan    = original_nan_mask.sum()

    # Step 1 & 2: mean absolute step-to-step change
    abs_deltas     = s.diff().abs()
    mean_abs_delta = abs_deltas.mean()

    # Accumulate a mask of all positions that get set to NaN (spikes)
    spike_mask = pd.Series(False, index=s.index)

    if mean_abs_delta == 0 or pd.isna(mean_abs_delta):
        print(f"  • QAQC ({label}): mean Δ = 0 or NaN – no spike removal performed")
        s.interpolate(method='time', inplace=True)
        s.ffill(inplace=True)
        s.bfill(inplace=True)
        replaced_mask = original_nan_mask
        return s, replaced_mask

    threshold = threshold_multiplier * mean_abs_delta

    # Step 3 & 4: iteratively flag spikes
    # We iterate because removing a spike changes the neighbours' deltas.
    flagged_total = 0
    max_passes    = 10
    for _pass in range(max_passes):
        abs_deltas = s.diff().abs()
        spikes     = abs_deltas > threshold
        n_spikes   = spikes.sum()
        if n_spikes == 0:
            break
        spike_mask |= spikes
        s[spikes]   = np.nan
        flagged_total += n_spikes

    # Step 5: linear interpolation over NaN gaps
    s.interpolate(method='time', inplace=True)
    s.ffill(inplace=True)   # fill any leading NaNs
    s.bfill(inplace=True)   # fill any trailing NaNs

    # Combined mask: every index that received an interpolated value
    replaced_mask = original_nan_mask | spike_mask

    total_replaced = int(replaced_mask.sum())
    print(f"  • QAQC ({label}): mean |Δ| = {mean_abs_delta:.4g}, "
          f"threshold = {threshold:.4g} "
          f"→ {flagged_total} spike(s) flagged, "
          f"{total_replaced} value(s) interpolated (incl. {n_original_nan} pre-existing NaN)")
    return s, replaced_mask


# ---------------------------------------------------------------------------

def get_nwps_forecast(station_id="BRLL1"):
    """
    Retrieve observed + forecasted stage/flow from the **National Water
    Prediction Service (NWPS)** API  (api.water.noaa.gov, 2025+ endpoint).

    The API returns separate ``"observation"`` and ``"forecast"`` blocks; this
    function merges them into a single chronological DataFrame.

    Parameters
    ----------
    station_id : str – NWPS gauge identifier (e.g. ``"BTRL1"``).

    Returns
    -------
    pd.DataFrame or None
        Columns: ``Stage_ft``, ``Flow_kcfs``, ``Type``, ``Flow_cfs``.
        Indexed by ``'Timestamp'``.  Returns ``None`` on failure.

    Notes
    -----
    • A ``User-Agent`` header is **required** by the NWPS API; omitting it
      results in HTTP 403.
    • Flow is reported in *thousands of cfs* (kcfs) by the API; this function
      adds a ``Flow_cfs`` column converted to cfs for consistency.
    """
    url = f"https://api.water.noaa.gov/nwps/v1/gauges/{station_id}/stageflow"

    # User-Agent is STRICTLY REQUIRED by the new API
    headers = {
        "User-Agent": "my-weather-script (contact@example.com)",
        "Accept":     "application/json",
    }

    print(f"Connecting to NWPS API for {station_id}...")

    try:
        response = requests.get(url, headers=headers, timeout=10)
        response.raise_for_status()
        data = response.json()

        # ---- 1. Parse Forecast Data --------------------------------------
        forecast_points = []
        if 'forecast' in data and 'data' in data['forecast']:
            for point in data['forecast']['data']:
                forecast_points.append({
                    'Timestamp': point.get('validTime'),
                    'Stage_ft':  point.get('primary'),    # primary = stage
                    'Flow_kcfs': point.get('secondary'),  # secondary = flow (kcfs)
                    'Type':      'Forecast',
                })

        # ---- 2. Parse Observed Data (past) --------------------------------
        observed_points = []
        if 'observation' in data and 'data' in data['observation']:
            for point in data['observation']['data']:
                observed_points.append({
                    'Timestamp': point.get('validTime'),
                    'Stage_ft':  point.get('primary'),
                    'Flow_kcfs': point.get('secondary'),
                    'Type':      'Observed',
                })

        # ---- 3. Combine & build DataFrame ---------------------------------
        all_data = observed_points + forecast_points
        if not all_data:
            print("No data found in JSON response.")
            return None

        df = pd.DataFrame(all_data)
        df['Timestamp'] = pd.to_datetime(df['Timestamp'])
        df.set_index('Timestamp', inplace=True)
        df.sort_index(inplace=True)

        # Convert kcfs → cfs (handle None / NaN gracefully)
        if 'Flow_kcfs' in df.columns:
            df['Flow_cfs'] = df['Flow_kcfs'].apply(
                lambda x: x * 1000.0 if pd.notnull(x) else None
            )

        return df

    except requests.exceptions.HTTPError as e:
        # 404 → wrong Station ID;  403 → missing User-Agent
        print(f"HTTP Error: {e}")
    except Exception as e:
        print(f"General Error: {e}")

    return None


# ---------------------------------------------------------------------------

def get_noaa_predictions_rest(station_id, days_forward=31,
                              datum="MLLW", units="english",
                              time_zone="lst_ldt"):
    """
    Fetch **predicted** tidal water levels directly from the NOAA CO-OPS REST
    API, as a fallback when the ``noaa_coops`` library fails.

    The request window runs from tomorrow through *days_forward* days ahead,
    split into ≤ 30-day chunks to respect the API limit.

    Parameters
    ----------
    station_id   : str  – NOAA CO-OPS station ID (e.g. ``"8760922"``).
    days_forward : int  – Total number of future days to retrieve (default 31).
    datum        : str  – Tidal datum (default ``"MLLW"``).
    units        : str  – ``"english"`` or ``"metric"`` (default ``"english"``).
    time_zone    : str  – ``"gmt"`` or ``"lst_ldt"`` (default ``"lst_ldt"``).

    Returns
    -------
    pd.DataFrame or None
        DataFrame indexed by ``'Date Time'`` with a ``'Predicted'`` column.
        Returns ``None`` if no data could be retrieved.
    """
    base_url   = "https://api.tidesandcurrents.noaa.gov/api/prod/datagetter"
    now        = datetime.now()
    start_date = now + timedelta(days=1)
    end_date   = now + timedelta(days=days_forward)

    all_chunks    = []
    current_start = start_date

    while current_start < end_date:
        current_end = min(current_start + timedelta(days=30), end_date)
        params = {
            "begin_date":  current_start.strftime("%Y%m%d"),
            "end_date":    current_end.strftime("%Y%m%d"),
            "station":     station_id,
            "product":     "predictions",
            "datum":       datum,
            "units":       units,
            "time_zone":   time_zone,
            "interval":    "h",          # hourly predictions
            "format":      "csv",
            "application": "Python_Script",
        }
        try:
            resp = requests.get(base_url, params=params, timeout=15)
            if resp.status_code == 200 and "Error" not in resp.text and resp.text.strip():
                chunk = pd.read_csv(io.StringIO(resp.text))
                chunk.columns = chunk.columns.str.strip()
                all_chunks.append(chunk)
            else:
                print(f"    REST fallback: no data for "
                      f"{current_start.strftime('%Y%m%d')}–"
                      f"{current_end.strftime('%Y%m%d')}")
        except Exception as exc:
            print(f"    REST fallback error: {exc}")
        time.sleep(0.5)
        current_start = current_end + timedelta(days=1)

    if not all_chunks:
        return None

    df = pd.concat(all_chunks, ignore_index=True)
    df['Date Time'] = pd.to_datetime(df['Date Time'])
    df.set_index('Date Time', inplace=True)
    return df


# ---------------------------------------------------------------------------

def add_elapsed_days_xaxis(ax, ref_date, offset=-0.18):
    """
    Attach a secondary x-axis **below** *ax* that labels time as elapsed days
    after *ref_date*.  Matplotlib date numbers (days since 0001-01-01) are
    shifted by the reference-date number so that the reference date maps to 0.

    The secondary axis is placed at *offset* in axes-fraction coordinates
    (default -0.18), which sits below the primary axis's own tick labels and
    xlabel so that nothing overlaps.  ``tight_layout()`` at the call site
    expands the figure's bottom margin to fit.

    Parameters
    ----------
    ax       : matplotlib.axes.Axes – the primary (date) axis.
    ref_date : datetime             – the reference date (elapsed-day origin).
    offset   : float               – vertical position in axes coordinates;
                                     negative values place the axis below the
                                     plot area (default -0.18).

    Returns
    -------
    matplotlib.axes.Axes  – the secondary axis object.
    """
    ref_num = mdates.date2num(ref_date)
    ax2 = ax.secondary_xaxis(
        offset,
        functions=(lambda x: x - ref_num, lambda x: x + ref_num),
    )
    ax2.set_xlabel("Elapsed Days After Reference Date", labelpad=8)
    ax2.tick_params(axis='x', labelsize=8, pad=4)
    return ax2


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 4 — MAIN EXECUTION
# ═══════════════════════════════════════════════════════════════════════════════
# Everything below runs top-to-bottom when the script is executed directly.
# The work is broken into numbered steps that mirror the logical pipeline.
# ═══════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":

    # ╔═════════════════════════════════════════════════════════════════════╗
    # ║  STEP 4-A:  RETRIEVE OBSERVED WATER LEVELS  (NOAA CO-OPS)        ║
    # ╠═════════════════════════════════════════════════════════════════════╣
    # ║  Source : NOAA Tides & Currents REST API                         ║
    # ║  Station: Pilots Station East, SW Pass  (8760922)                ║
    # ║  Product: observed water_level                                   ║
    # ║  Output : wl_observed – DataFrame of 6-minute water levels       ║
    # ╚═════════════════════════════════════════════════════════════════════╝

    print("\n[STEP 1/8] Retrieving Observed Water Levels from NOAA CO-OPS...")
    print(f"  Station: {STATION_ID} (Pilots Station East, SW Pass, LA)")
    print(f"  Look-back period: {Days_Back} days")
    df = get_noaa_data(days_back=Days_Back)
    wl_observed = df
    if wl_observed is not None:
        print(f"  ✓ Retrieved {len(wl_observed):,} observed water level records\n")
    else:
        print("  ✗ Failed to retrieve observed water level data\n")

    # ╔═════════════════════════════════════════════════════════════════════╗
    # ║  STEP 4-B:  RETRIEVE PREDICTED WATER LEVELS  (noaa_coops)        ║
    # ╠═════════════════════════════════════════════════════════════════════╣
    # ║  Source : NOAA CO-OPS predictions via the `noaa_coops` library   ║
    # ║  Window : tomorrow → +31 days                                    ║
    # ║  Output : df_predictions – DataFrame of predicted tidal levels   ║
    # ╚═════════════════════════════════════════════════════════════════════╝

    print("[STEP 2/8] Retrieving Predicted Water Levels (Tides)...")
    df_predictions = None   # always defined; set to valid DataFrame on success
    try:
        now = datetime.now()
        prediction_begin_str = (now + timedelta(days=1)).strftime('%Y%m%d')
        prediction_end_str   = (now + timedelta(days=31)).strftime('%Y%m%d')

        station = Station(id=STATION_ID)

        df_predictions = station.get_data(
            begin_date=prediction_begin_str,
            end_date=prediction_end_str,
            product="predictions",
            datum=DATUM,
            units=UNITS,
            time_zone=TIME_ZONE,
        )

        print(f"  ✓ Retrieved {len(df_predictions):,} tidal prediction records")
        print(f"  Forecast window: {prediction_begin_str} to {prediction_end_str}\n")

    except Exception as e:
        print(f"  ⚠ noaa_coops library failed for station {STATION_ID}: {e}")
        print("  → Attempting direct NOAA CO-OPS REST API fallback...")
        try:
            df_predictions = get_noaa_predictions_rest(
                station_id=STATION_ID,
                days_forward=31,
                datum=DATUM,
                units=UNITS,
                time_zone=TIME_ZONE,
            )
            if df_predictions is not None and not df_predictions.empty:
                print(f"  ✓ REST fallback: retrieved {len(df_predictions):,} "
                      f"tidal prediction records\n")
            else:
                print("  ✗ REST fallback also returned no data. "
                      "Stage forecast will be omitted.\n")
                df_predictions = None
        except Exception as e2:
            print(f"  ✗ REST fallback failed: {e2}\n")
            df_predictions = None

    # ╔═════════════════════════════════════════════════════════════════════╗
    # ║  STEP 4-C:  COMBINE OBSERVED + PREDICTED STAGE                   ║
    # ╠═════════════════════════════════════════════════════════════════════╣
    # ║  Both datasets are resampled to 1-hour frequency and converted   ║
    # ║  from MLLW to NAVD88.  A combined plot is saved as a PNG.        ║
    # ║  Output : combined_df – hourly water-level DataFrame (NAVD88)    ║
    # ╚═════════════════════════════════════════════════════════════════════╝

    # Identify the observed-data value column and resample to hourly
    obs_col    = 'Water Level' if 'Water Level' in wl_observed.columns else wl_observed.columns[0]
    obs_hourly = wl_observed[[obs_col]].resample('1h').mean()

    # Resample predicted data to hourly means (rename column to match observed)
    if df_predictions is not None and not df_predictions.empty:
        # Determine the prediction column: noaa_coops uses 'v'; REST API uses 'Predicted'
        _pred_candidates = [c for c in ('v', 'Predicted', 'prediction')
                            if c in df_predictions.columns]
        pred_col = _pred_candidates[0] if _pred_candidates else df_predictions.columns[0]
        df_predictions.index = pd.to_datetime(df_predictions.index)
        pred_hourly = (
            df_predictions[[pred_col]]
            .rename(columns={pred_col: obs_col})
            .resample('1h')
            .mean()
        )
    else:
        print("  ⚠ No tidal prediction data available – stage boundary condition "
              "will use observed data only")
        pred_hourly = pd.DataFrame(columns=[obs_col])

    # ---- Apply Stage Forecast Calibration --------------------------------
    if isinstance(Stage_Forecast_Calibration, bool):
        if Stage_Forecast_Calibration:
            obs_mean  = obs_hourly[obs_col].mean()
            pred_mean = pred_hourly[obs_col].mean()
            stage_offset = obs_mean - pred_mean
            pred_hourly[obs_col] = pred_hourly[obs_col] + stage_offset
            print(f"  ✓ Stage Forecast Calibration (auto): shifted forecast by "
                  f"{stage_offset:+.3f} ft to match observed mean ({obs_mean:.3f} ft)")
        else:
            print("  • Stage Forecast Calibration: disabled")
    elif isinstance(Stage_Forecast_Calibration, str) and Stage_Forecast_Calibration.strip().upper() == 'STM':
        # Short Term Match: align the median of the last 24 h of observed stage
        # with the median of the first 24 h of forecast stage.
        _obs_last24  = obs_hourly[obs_col].last('24h')
        _pred_first24 = pred_hourly[obs_col].first('24h')
        if _obs_last24.empty or _pred_first24.empty:
            print("  ⚠ Stage Forecast Calibration (STM): insufficient data for "
                  "short-term match – calibration skipped")
        else:
            _obs_median  = _obs_last24.median()
            _pred_median = _pred_first24.median()
            stage_offset = _obs_median - _pred_median
            pred_hourly[obs_col] = pred_hourly[obs_col] + stage_offset
            print(f"  ✓ Stage Forecast Calibration (STM): shifted forecast by "
                  f"{stage_offset:+.3f} ft "
                  f"(obs last-24 h median={_obs_median:.3f} ft, "
                  f"fcst first-24 h median={_pred_median:.3f} ft)")
    elif isinstance(Stage_Forecast_Calibration, (int, float)):
        pred_hourly[obs_col] = pred_hourly[obs_col] + Stage_Forecast_Calibration
        print(f"  ✓ Stage Forecast Calibration (manual): applied additive offset "
              f"of {Stage_Forecast_Calibration:+.3f} ft to forecast stage")

    # Concatenate and apply the MLLW → NAVD88 datum shift
    print("[STEP 3/8] Combining & Processing Water Level Data...")
    combined_df = pd.concat([obs_hourly, pred_hourly]).sort_index() + MLLW2NAVD88
    UNITS = "ft NAVD88"   # update label to reflect converted datum
    print(f"  ✓ Combined {len(combined_df):,} hourly records (observed + forecast)")
    print(f"  ✓ Converted datum: MLLW → NAVD88 ({MLLW2NAVD88:+.2f} ft offset)")

    # ---- Plot & save combined water level --------------------------------
    plt.figure(figsize=(14, 7))
    plt.plot(combined_df.index, combined_df[obs_col],
             label='Hourly Water Level', color='navy')

    # Shade the forecast period
    if not pred_hourly.empty:
        plt.axvspan(pred_hourly.index.min(), pred_hourly.index.max(),
                    color='lightblue', alpha=0.5, label='Forecast Period')

    plt.title(f"Combined Hourly River Stage - Station {STATION_ID}")
    plt.ylabel(f"Water Level ({UNITS} {DATUM})")
    plt.xlabel("Date Time")
    plt.legend()
    plt.grid(True)

    # Save instead of showing interactively
    sim_folder = folder
    png_path   = sim_folder / f"Water_Level_BC_{sim_folder.name}.png"
    add_elapsed_days_xaxis(plt.gca(), reference_date)
    plt.tight_layout()
    plt.savefig(png_path, dpi=150)
    plt.close()
    print(f"  ✓ Saved water level plot: {png_path.name}\n")

    # ╔═════════════════════════════════════════════════════════════════════╗
    # ║  STEP 4-D:  RETRIEVE OBSERVED INSTANTANEOUS DISCHARGE  (USGS)    ║
    # ╠═════════════════════════════════════════════════════════════════════╣
    # ║  Source : USGS NWIS Instantaneous-Value API                      ║
    # ║  Station: Mississippi River at Baton Rouge (07374000)             ║
    # ║  Output : usgs_df – high-frequency (15-min) discharge            ║
    # ╚═════════════════════════════════════════════════════════════════════╝

    print("[STEP 4/8] Retrieving Observed Discharge from USGS NWIS...")
    print(f"  Station: {USGS_STATION_ID} (Mississippi River at Baton Rouge, LA)")
    usgs_df = get_usgs_instantaneous_discharge(USGS_STATION_ID, Days_Back)

    if not usgs_df.empty:
        print(f"  ✓ Downloaded {len(usgs_df):,} instantaneous discharge records\n")

    # ╔═════════════════════════════════════════════════════════════════════╗
    # ║  STEP 4-E:  DEVELOP STAGE-DISCHARGE RATING CURVE  (optional)     ║
    # ╠═════════════════════════════════════════════════════════════════════╣
    # ║  Controlled by the `generate_RC` flag (Section 1).               ║
    # ║  Uses 5 years of daily stage + flow from USGS to fit a 2nd-order ║
    # ║  polynomial.  A separate "low-water" fit (< 450 000 cfs) is also ║
    # ║  produced.                                                        ║
    # ║  Outputs: coeffs      – global curve [a, b, c]                   ║
    # ║           coeffs_low  – low-water curve [a, b, c]                ║
    # ╚═════════════════════════════════════════════════════════════════════╝

    if generate_RC == 1:
        print("[STEP 5/8] Developing Stage-Discharge Rating Curve...")
        print(f"  Using {5} years of USGS daily data for polynomial regression")
        RC_STATION   = USGS_STATION_ID       # same gauge for rating curve
        RC_YEARS     = 5
        RC_DAYS_BACK = RC_YEARS * 365

        # ---- Acquire stage & flow daily data ----------------------------
        df_flow  = get_usgs_daily_data(RC_STATION, RC_DAYS_BACK, '00060')
        df_flow.rename(columns={'Value': 'Flow_cfs'}, inplace=True)

        df_stage = get_usgs_daily_data(RC_STATION, RC_DAYS_BACK, '00065')
        df_stage.rename(columns={'Value': 'Stage_ft'}, inplace=True)

        if not df_flow.empty and not df_stage.empty:
            # Merge on date index (inner join keeps only overlapping days)
            df_merged = pd.merge(df_stage, df_flow,
                                 left_index=True, right_index=True, how='inner')
            print(f"  ✓ Merged {len(df_merged):,} daily stage-discharge pairs")

            df_clean = df_merged.dropna()
            x = df_clean['Stage_ft']
            y = df_clean['Flow_cfs']

            # ---- Global polynomial regression (degree 2) ----------------
            #   Flow = a * Stage^2 + b * Stage + c
            coeffs   = np.polyfit(x, y, 2)
            poly_eqn = np.poly1d(coeffs)

            print("\n  Global Rating Curve (2nd-order polynomial):")
            print(f"    Flow = {coeffs[0]:.2f}·Stage² + {coeffs[1]:.2f}·Stage + {coeffs[2]:.2f}")

            x_line = np.linspace(x.min(), x.max(), 100)
            y_line = poly_eqn(x_line)

            # ---- Low-water regression (< 450 000 cfs) -------------------
            df_low     = df_clean[df_clean['Flow_cfs'] < 450000]
            y_line_low = None

            if not df_low.empty:
                x_low = df_low['Stage_ft']
                y_low = df_low['Flow_cfs']

                coeffs_low   = np.polyfit(x_low, y_low, 2)
                poly_eqn_low = np.poly1d(coeffs_low)

                print("\n  Low-Water Rating Curve (< 450,000 cfs):")
                print(f"    Flow = {coeffs_low[0]:.2f}·Stage² + {coeffs_low[1]:.2f}·Stage + {coeffs_low[2]:.2f}")
                print(f"  ✓ Rating curves developed successfully\n")

                x_line_low = np.linspace(x_low.min(), x_low.max(), 100)
                y_line_low = poly_eqn_low(x_line_low)
            else:
                print("\nNot enough low water data for specific regression.")
        else:
            print("Error: Could not retrieve necessary data to build rating curve.")

    # ╔═════════════════════════════════════════════════════════════════════╗
    # ║  STEP 4-F:  RETRIEVE NWPS FORECAST  (stage + flow)               ║
    # ╠═════════════════════════════════════════════════════════════════════╣
    # ║  Source : NOAA NWPS API                                          ║
    # ║  Station: BTRL1 (Baton Rouge)                                    ║
    # ║  Output : df_forecast – forecast-only rows from the NWPS data    ║
    # ║           calc_flow   – discharge computed from forecast stage    ║
    # ║                         using the rating-curve coefficients       ║
    # ╚═════════════════════════════════════════════════════════════════════╝

    print("[STEP 6/8] Retrieving Forecast Stage/Flow from NWPS...")
    STATION_ID = NWPS_STATION_ID           # switch context to NWPS gauge
    df = get_nwps_forecast(STATION_ID)

    if df is not None:
        obs_count = len(df[df['Type'] == 'Observed'])
        fcst_count = len(df[df['Type'] == 'Forecast'])
        print(f"  ✓ Retrieved {obs_count:,} observed + {fcst_count:,} forecast records")

        # Save raw NWPS data to CSV
        # filename = f"NWPS_{STATION_ID}_Combined.csv"
        # df.to_csv(filename)
        # print(f"  ✓ Saved NWPS data: {filename}")

        # Isolate forecast rows
        df_forecast = df[df['Type'] == 'Forecast'].copy()

        # ---- Compute calc_flow from forecast stage using rating curve ----
        #   Prefer low-water coefficients; fall back to global fit.
        try:
            if 'Stage_ft' in df_forecast.columns and not df_forecast['Stage_ft'].isna().all():
                stages = df_forecast['Stage_ft'].astype(float)

                if 'coeffs_low' in globals():
                    calc_flow = (coeffs_low[0] * stages**2
                                 + coeffs_low[1] * stages
                                 + coeffs_low[2])
                    print("  ✓ Converted forecast stage to discharge (low-water curve)\n")
                elif 'coeffs' in globals():
                    calc_flow = (coeffs[0] * stages**2
                                 + coeffs[1] * stages
                                 + coeffs[2])
                    print("  ✓ Converted forecast stage to discharge (global curve)\n")
                else:
                    calc_flow = pd.Series(index=stages.index, dtype=float)
                    print("Warning: No rating curve coefficients found; "
                          "'calc_flow' will be empty.")
            else:
                calc_flow = pd.Series(dtype=float)
                print("Warning: 'df_forecast' has no Stage_ft values; "
                      "'calc_flow' is empty.")
        except Exception as e:
            calc_flow = pd.Series(dtype=float)
            print(f"Error calculating 'calc_flow': {e}")
    else:
        print("Failed to retrieve NWPS data.")

    # ╔═════════════════════════════════════════════════════════════════════╗
    # ║  STEP 4-G:  COMBINE OBSERVED + FORECAST DISCHARGE & PLOT         ║
    # ╠═════════════════════════════════════════════════════════════════════╣
    # ║  Merges USGS instantaneous discharge (resampled daily) with the  ║
    # ║  rating-curve-derived forecast discharge.  `combine_first` keeps ║
    # ║  observed values where the two overlap.                          ║
    # ║  Output : combined_daily – daily discharge Series                ║
    # ╚═════════════════════════════════════════════════════════════════════╝

    if 'usgs_df' in globals() and 'calc_flow' in globals():
        print("[STEP 7/8] Combining Observed & Forecast Discharge...")

        # ---- Observed → daily mean --------------------------------------
        obs_daily = usgs_df['Discharge_cfs'].resample('D').mean()

        # ---- Forecast → daily mean (strip timezone first) ----------------
        forecast_series = calc_flow.copy()
        if forecast_series.index.tz is not None:
            forecast_series.index = forecast_series.index.tz_convert(None)
        forecast_daily = forecast_series.resample('D').mean()

        # ---- Apply Discharge Forecast Calibration ------------------------
        #   Scaling ramps linearly from 10 % at day +1 to 100 % at day +7.
        #   effective_scale = 1 + weight(d) × (factor − 1)
        if Discharge_Forecast_Calibration != 1.0:
            _cal = np.clip(Discharge_Forecast_Calibration, 0.5, 2.0)
            now_date = pd.Timestamp.now().normalize()
            days_ahead = np.maximum((forecast_daily.index - now_date).days, 1)
            weight = np.minimum(1.0, 0.10 + (days_ahead - 1) * 0.90 / 6.0)
            effective_scale = 1.0 + weight * (_cal - 1.0)
            forecast_daily = forecast_daily * effective_scale
            print(f"  ✓ Discharge Forecast Calibration: factor={_cal:.2f} "
                  f"(10 % weight at day +1 → 100 % at day +7)")
        else:
            print("  • Discharge Forecast Calibration: factor=1.00 (no scaling)")

        # ---- Merge: observed takes priority on overlapping dates ----------
        combined_daily = obs_daily.combine_first(forecast_daily).sort_index()

        # ---- Extend to 30 days from today via linear extrapolation -------
        # The NWPS forecast typically covers only ~28 days; fill any gap up
        # to forecast_end_date (today + 30 days) using a linear trend fitted
        # to the last 7 available daily values.
        _d30_end = pd.Timestamp(forecast_end_date.date())
        if combined_daily.index.max() < _d30_end:
            _tail = combined_daily.dropna().tail(7)
            if len(_tail) >= 2:
                _xi     = np.arange(len(_tail), dtype=float)
                _yi     = _tail.values.astype(float)
                _slope, _intercept = np.polyfit(_xi, _yi, 1)
                _missing_dates = pd.date_range(
                    start=combined_daily.index.max() + timedelta(days=1),
                    end=_d30_end,
                    freq='D',
                )
                _n_extra    = len(_missing_dates)
                _extra_xi   = np.arange(len(_tail), len(_tail) + _n_extra, dtype=float)
                _extra_vals = _slope * _extra_xi + _intercept
                _extra_series = pd.Series(_extra_vals, index=_missing_dates)
                combined_daily = pd.concat([combined_daily, _extra_series]).sort_index()
                # Extend forecast_daily mask so the shading covers extrapolated days
                forecast_daily = pd.concat([forecast_daily, _extra_series]).sort_index()
                print(f"  ✓ Linearly extrapolated discharge for {_n_extra} day(s) "
                      f"beyond NWPS data → {_d30_end.strftime('%Y-%m-%d')} "
                      f"(slope={_slope:+.0f} cfs/day)")
            else:
                print("  ⚠ Not enough data points for linear extrapolation of discharge")

        # ---- Plot & save combined discharge ------------------------------
        plt.figure(figsize=(14, 7))
        plt.plot(combined_daily.index, combined_daily,
                 color='navy', linewidth=2, label='Combined Discharge')

        # Shade forecast window
        if not forecast_daily.empty:
            forecast_start = forecast_daily.index.min()
            plt.axvspan(forecast_start, combined_daily.index.max(),
                        color='lightblue', alpha=0.5, label='Forecast Period')

        plt.title(f"Combined Daily River Discharge (USGS Observed + NWPS Forecast)")
        plt.ylabel("Discharge (cfs)")
        plt.xlabel("Date")
        plt.legend()
        plt.grid(True, linestyle=':', alpha=0.6)

        # Format Y-axis labels with commas (e.g. 250,000)
        plt.gca().get_yaxis().set_major_formatter(
            ticker.FuncFormatter(lambda x, p: format(int(x), ','))
        )

        # Save instead of showing interactively
        sim_folder = folder
        png_path   = sim_folder / f"Discharge_BC_{sim_folder.name}.png"
        add_elapsed_days_xaxis(plt.gca(), reference_date)
        plt.tight_layout()
        plt.savefig(png_path, dpi=150)
        plt.close()
        print(f"  ✓ Combined {len(combined_daily):,} daily records (observed + forecast)")
        print(f"  ✓ Saved discharge plot: {png_path.name}\n")

    else:
        print("Error: Necessary data ('usgs_df' or 'calc_flow') is missing. "
              "Please ensure previous steps executed successfully.")

    # ╔═════════════════════════════════════════════════════════════════════╗
    # ║  STEP 4-H:  EXPORT BOUNDARY-CONDITION FILES                      ║
    # ╠═════════════════════════════════════════════════════════════════════╣
    # ║  Two tab-delimited text files are written to the simulation       ║
    # ║  folder for direct ingestion by the MRSWAT2000 Fortran model:    ║
    # ║                                                                   ║
    # ║    Stages.txt      – columns: TIME (elapsed days)  |  STAGE      ║
    # ║    Discharges.txt  – columns: TIME (elapsed days)  |  DISCHARGE  ║
    # ║                                                                   ║
    # ║  TIME is measured in fractional days from the earliest timestamp  ║
    # ║  in each respective dataset.                                      ║
    # ╚═════════════════════════════════════════════════════════════════════╝

    print("[STEP 8/8] Exporting Model Boundary-Condition Files...")
    try:
        sim_folder = folder

        # ╔═══════════════════════════════════════════════════════════════╗
        # ║  QAQC — spike detection & linear interpolation              ║
        # ║  Applied to both the stage and discharge series before      ║
        # ║  writing to file.  Any value whose step-to-step change      ║
        # ║  exceeds 2.5 × the mean absolute change is flagged as       ║
        # ║  erroneous, set to NaN, and replaced by linear              ║
        # ║  interpolation from the surrounding valid values.           ║
        # ╚═══════════════════════════════════════════════════════════════╝
        print("\n  [QAQC] Running quality control on time series...")

        stage_replaced_mask    = None
        discharge_replaced_mask = None

        if 'combined_df' in globals() and combined_df is not None and not combined_df.empty:
            _col = combined_df.columns[0]
            combined_df[_col], stage_replaced_mask = qaqc_timeseries(
                combined_df[_col],
                threshold_multiplier=2.5,
                label="Stage (ft NAVD88)",
            )

        if 'combined_daily' in globals() and combined_daily is not None and not combined_daily.empty:
            combined_daily, discharge_replaced_mask = qaqc_timeseries(
                combined_daily,
                threshold_multiplier=2.5,
                label="Discharge (cfs)",
            )

        print("  [QAQC] Quality control complete\n")

        # ---- Regenerate plots with replaced values highlighted -----------
        # Both PNGs are overwritten so the saved files always reflect the
        # QA-corrected data.  Interpolated / replaced points are shown as
        # filled orange circles so they stand out clearly from the main line.

        # -- Stage plot ---------------------------------------------------
        if (combined_df is not None and not combined_df.empty
                and stage_replaced_mask is not None):
            _col      = combined_df.columns[0]
            _obs_end  = pred_hourly.index.min() if 'pred_hourly' in globals() and not pred_hourly.empty else None

            fig, ax = plt.subplots(figsize=(14, 7))
            ax.plot(combined_df.index, combined_df[_col],
                    color='navy', linewidth=1.0, label='Hourly Water Level', zorder=2)

            if _obs_end is not None:
                ax.axvspan(pred_hourly.index.min(), combined_df.index.max(),
                           color='lightblue', alpha=0.5, label='Forecast Period', zorder=1)

            # Overlay interpolated/replaced points
            _rep = stage_replaced_mask[stage_replaced_mask]
            if not _rep.empty:
                ax.scatter(_rep.index, combined_df.loc[_rep.index, _col],
                           color='orangered', s=20, zorder=5,
                           label=f'Interpolated / replaced ({len(_rep):,})')

            ax.set_title(f"Combined Hourly River Stage – Station {NWPS_STATION_ID} [QAQC Applied]")
            ax.set_ylabel(f"Water Level (ft NAVD88)")
            ax.set_xlabel("Date Time")
            ax.legend()
            ax.grid(True)
            add_elapsed_days_xaxis(ax, reference_date)
            fig.tight_layout()
            png_path = sim_folder / f"Water_Level_BC_{sim_folder.name}.png"
            fig.savefig(png_path, dpi=150)
            plt.close(fig)
            print(f"  ✓ Regenerated stage plot with QAQC highlights: {png_path.name}")

        # -- Discharge plot -----------------------------------------------
        if (combined_daily is not None and not combined_daily.empty
                and discharge_replaced_mask is not None):
            _fcst_start = forecast_daily.index.min() if 'forecast_daily' in globals() and not forecast_daily.empty else None

            fig, ax = plt.subplots(figsize=(14, 7))
            ax.plot(combined_daily.index, combined_daily,
                    color='navy', linewidth=2, label='Combined Discharge', zorder=2)

            if _fcst_start is not None:
                ax.axvspan(_fcst_start, combined_daily.index.max(),
                           color='lightblue', alpha=0.5, label='Forecast Period', zorder=1)

            # Overlay interpolated/replaced points
            _rep = discharge_replaced_mask[discharge_replaced_mask]
            if not _rep.empty:
                ax.scatter(_rep.index, combined_daily.loc[_rep.index],
                           color='orangered', s=60, zorder=5,
                           label=f'Interpolated / replaced ({len(_rep):,})')

            ax.set_title("Combined Daily River Discharge – USGS Observed + NWPS Forecast [QAQC Applied]")
            ax.set_ylabel("Discharge (cfs)")
            ax.set_xlabel("Date")
            ax.legend()
            ax.grid(True, linestyle=':', alpha=0.6)
            ax.get_yaxis().set_major_formatter(
                ticker.FuncFormatter(lambda x, p: format(int(x), ','))
            )
            add_elapsed_days_xaxis(ax, reference_date)
            fig.tight_layout()
            png_path = sim_folder / f"Discharge_BC_{sim_folder.name}.png"
            fig.savefig(png_path, dpi=150)
            plt.close(fig)
            print(f"  ✓ Regenerated discharge plot with QAQC highlights: {png_path.name}\n")

        # ---- Stages.txt --------------------------------------------------
        if 'combined_df' in globals() and combined_df is not None and not combined_df.empty:
            stages = combined_df.copy()

            # Ensure a naive (timezone-unaware) datetime index
            if getattr(stages.index, 'tz', None) is not None:
                stages.index = stages.index.tz_convert(None)

            ref_time     = stages.index.min()
            elapsed_days = (stages.index - ref_time).total_seconds() / 86400.0
            stage_values = stages.iloc[:, 0].values

            out_stages = pd.DataFrame({'TIME': elapsed_days, 'STAGE': stage_values})
            stages_path = sim_folder / 'Stages.txt'
            out_stages.to_csv(stages_path, sep='\t', index=False, float_format='%.6f')
            print(f"  ✓ Exported Stages.txt ({len(out_stages):,} hourly records)")
        else:
            print("Skipping Stages.txt (combined_df missing or empty)")

        # ---- Discharges.txt ----------------------------------------------
        if 'combined_daily' in globals() and combined_daily is not None and not combined_daily.empty:
            dis = combined_daily.copy()

            if getattr(dis.index, 'tz', None) is not None:
                dis.index = dis.index.tz_convert(None)

            ref_time_d     = dis.index.min()
            elapsed_days_d = (dis.index - ref_time_d).total_seconds() / 86400.0
            discharge_vals = dis.values

            out_dis = pd.DataFrame({
                'TIME':      elapsed_days_d,
                'DISCHARGE': discharge_vals.flatten(),
            })
            dis_path = sim_folder / 'Discharges.txt'
            out_dis.to_csv(dis_path, sep='\t', index=False, float_format='%.6f')
            print(f"  ✓ Exported Discharges.txt ({len(out_dis):,} daily records)")
        else:
            print("Skipping Discharges.txt (combined_daily missing or empty)")

    except Exception as e:
        print(f"Error writing time series files: {e}")

