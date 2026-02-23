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

Stage_Forecast_Calibration = True
    # Controls how the forecasted stage time-series is adjusted:
    #   True   → shift the forecast so its temporal mean matches the
    #            observed-stage temporal mean (bias correction).
    #   False  → no adjustment applied to forecast stage.
    #   <float> → additive offset (ft) applied to the entire forecast
    #            stage time-series (e.g. 0.25 adds 0.25 ft).


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 2 — BOOKKEEPING  (create date-stamped simulation folder)
# ═══════════════════════════════════════════════════════════════════════════════
# A new subfolder is created for every run so that outputs never overwrite a
# previous simulation.  The folder name includes the date and a zero-padded
# sequence number (e.g. MRSWAT_20260212_00, MRSWAT_20260212_01, …).

today_str      = datetime.now().strftime('%Y%m%d')
reference_date = datetime.now() - timedelta(days=Days_Back)  # used later for elapsed-time export

i = 0
while True:
    folder = Path(f"MRSWAT_{today_str}_{i:02d}")
    if not folder.exists():
        folder.mkdir()
        print("\n" + "="*80)
        print(f"  MRSWAT2000 Pre-Processor – {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print("="*80)
        print(f"  ✓ Created simulation folder: {folder}\n")
        break
    i += 1

# -----------------------------------------------------------------------------
# Write Date_Information.txt into the simulation folder
# Columns: Ref_Date<TAB>Sim_Date
# Row 2 : reference_date (YYYYMMDD) <TAB> current date (YYYYMMDD)
# -----------------------------------------------------------------------------
try:
    date_info_path = folder / 'Date_Information.txt'
    with open(date_info_path, 'w') as f:
        f.write('Ref_Date\tSim_Date\n')
        f.write(f"{reference_date.strftime('%Y%m%d')}\t{today_str}\n")
    print(f"  ✓ Wrote Date_Information.txt: {date_info_path.name}\n")
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
        print(f"Could not retrieve or plot predictions for station {STATION_ID}.")
        print(f"Error detail: {e}")

    # ╔═════════════════════════════════════════════════════════════════════╗
    # ║  STEP 4-C:  COMBINE OBSERVED + PREDICTED STAGE                   ║
    # ╠═════════════════════════════════════════════════════════════════════╣
    # ║  Both datasets are resampled to 1-hour frequency and converted   ║
    # ║  from MLLW to NAVD88.  A combined plot is saved as a PNG.        ║
    # ║  Output : combined_df – hourly water-level DataFrame (NAVD88)    ║
    # ╚═════════════════════════════════════════════════════════════════════╝

    # Identify the value column in each DataFrame
    obs_col  = 'Water Level' if 'Water Level' in wl_observed.columns else wl_observed.columns[0]
    pred_col = 'v' if 'v' in df_predictions.columns else df_predictions.columns[0]

    # Resample observed data to hourly means
    obs_hourly = wl_observed[[obs_col]].resample('1h').mean()

    # Resample predicted data to hourly means (rename column to match observed)
    df_predictions.index = pd.to_datetime(df_predictions.index)
    pred_hourly = (
        df_predictions[[pred_col]]
        .rename(columns={pred_col: obs_col})
        .resample('1h')
        .mean()
    )

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
        filename = f"NWPS_{STATION_ID}_Combined.csv"
        df.to_csv(filename)
        print(f"  ✓ Saved NWPS data: {filename}")

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


# -----------------------------------------------------------------------------
# Copy required template files into the newly created simulation folder
# The files are copied from the project's `Template` directory into the
# simulation folder so the MRSWAT Fortran executable and auxiliary files
# are present for downstream workflow steps.
# -----------------------------------------------------------------------------
print("\nCopying template files to simulation folder...")
try:
    template_dir = Path.cwd() / 'Template'
    files_to_copy = [
        'MRSWAT-111825.exe',
        'river-miles-for-selected-output.txt',
        'IdealizedClassicSill.txt',
        'mr-diversions.txt',
        'input.txt',
    ]

    copied_count = 0
    for fname in files_to_copy:
        src = template_dir / fname
        dst = folder / fname
        if src.exists():
            shutil.copy2(src, dst)
            copied_count += 1
        else:
            print(f"  ⚠ Template file not found (skipping): {src}")
    print(f"  ✓ Copied {copied_count} template files")
    # -----------------------------------------------------------------
    # Copy `main.inp`: prefer the latest previous simulation folder's
    # `main.inp` if any exist; otherwise fall back to Template/main.inp
    # -----------------------------------------------------------------
    cwd = Path.cwd()
    sim_dirs = [p for p in cwd.iterdir() if p.is_dir() and p.name.startswith('MRSWAT_') and p.name != folder.name]
    src_main = None
    if sim_dirs:
        sim_dirs_sorted = sorted(sim_dirs)
        latest_prev = sim_dirs_sorted[-1]
        candidate = latest_prev / 'main.inp'
        if candidate.exists():
            src_main = candidate
            print(f"Found previous simulation folder: {latest_prev.name}; using its main.inp")
        else:
            src_main = template_dir / 'main.inp'
            print(f"Previous simulation folder found but no main.inp there; using Template/main.inp")
    else:
        src_main = template_dir / 'main.inp'
        print("No previous simulation folders found; using Template/main.inp")

    dst_main = folder / 'main.inp'
    if src_main and src_main.exists():
        shutil.copy2(src_main, dst_main)
        print(f"  ✓ Copied main.inp from: {src_main.parent.name}/{src_main.name}")
    else:
        print(f"  ✗ main.inp source not found: {src_main}")

    # -------------------------------------------------------------------------
    # Extract hotstart conditions from previous simulation if available
    # -------------------------------------------------------------------------
    print("\nExtracting hotstart conditions from previous simulation...")
    if sim_dirs and latest_prev:
        try:
            # Read Date_Information.txt to get reference date
            date_info_path = latest_prev / 'Date_Information.txt'
            output_path = latest_prev / 'mr-full-output.txt'
            
            if date_info_path.exists() and output_path.exists():
                # Parse Date_Information.txt
                date_df = pd.read_csv(date_info_path, sep='\t')
                ref_date_str = str(date_df['Ref_Date'].iloc[0])
                ref_date = pd.to_datetime(ref_date_str, format='%Y%m%d')
                sim_date_str = str(date_df['Sim_Date'].iloc[0])
                sim_date = pd.to_datetime(sim_date_str, format='%Y%m%d')
                print(f"  ✓ Found reference date: {ref_date.strftime('%Y-%m-%d')}")
                print(f"  ✓ Found simulation date: {sim_date.strftime('%Y-%m-%d')}")
                
                # Read the entire output file
                with open(output_path, 'r') as f:
                    lines = f.readlines()
                
                # Parse matrices and their times
                matrices = []
                current_time = None
                current_matrix = []
                in_matrix = False
                header_lines = 4
                
                for i, line in enumerate(lines):
                    line_stripped = line.strip()
                    
                    # Check for time marker
                    if line_stripped.startswith('Time ='):
                        # Save previous matrix if exists
                        if current_time is not None and current_matrix:
                            matrices.append({
                                'time': current_time,
                                'matrix': current_matrix.copy()
                            })
                            current_matrix = []
                        
                        # Extract new time value
                        time_str = line_stripped.split('=')[1].strip()
                        current_time = float(time_str)
                        in_matrix = False
                        continue
                    
                    # Skip header lines after time marker
                    if current_time is not None and not in_matrix:
                        if i > 0 and lines[i-1].strip().startswith('Time ='):
                            # Count down header lines
                            header_count = 0
                            for j in range(i, min(i+10, len(lines))):
                                if lines[j].strip() and not lines[j].strip().replace('.','').replace('-','').replace('+','').replace('E','').replace('e','').replace(' ','').replace('\t','').isdigit():
                                    header_count += 1
                                else:
                                    break
                            if header_count >= 4:
                                continue
                    
                    # Check if this is a data line (contains numbers)
                    if line_stripped and current_time is not None:
                        # Try to parse as numerical data
                        parts = line_stripped.split()
                        if len(parts) > 0:
                            try:
                                # Check if first element is a number
                                float(parts[0])
                                current_matrix.append(line_stripped)
                                in_matrix = True
                            except ValueError:
                                # Not a data line (probably header)
                                pass
                
                # Save last matrix
                if current_time is not None and current_matrix:
                    matrices.append({
                        'time': current_time,
                        'matrix': current_matrix.copy()
                    })
                
                print(f"  ✓ Parsed {len(matrices)} output matrices from previous simulation")
                
                # Convert times to dates and find matrix matching the simulation date
                sim_date_target = sim_date.date()
                target_matrix = None
                closest_diff = None
                
                for mat in matrices:
                    mat_date = (ref_date + pd.Timedelta(days=mat['time'])).date()
                    mat['date'] = mat_date
                    
                    # Find closest matrix to the simulation date
                    diff = abs((mat_date - sim_date_target).days)
                    if closest_diff is None or diff < closest_diff:
                        closest_diff = diff
                        target_matrix = mat
                
                if target_matrix:
                    hotstart_path = folder / 'hotstart.txt'
                    with open(hotstart_path, 'w') as f:
                        for line in target_matrix['matrix']:
                            f.write(line + '\n')
                    
                    print(f"  ✓ Extracted hotstart conditions from: {target_matrix['date'].strftime('%Y-%m-%d')}")
                    print(f"    (Time = {target_matrix['time']:.4f} days, matching Sim_Date from previous run)")
                    print(f"  ✓ Saved hotstart.txt ({len(target_matrix['matrix'])} rows)")
                else:
                    print("  ⚠ No suitable matrix found for hotstart")
            else:
                if not date_info_path.exists():
                    print(f"  ⚠ Date_Information.txt not found in {latest_prev.name}")
                if not output_path.exists():
                    print(f"  ⚠ mr-full-output.txt not found in {latest_prev.name}")
        except Exception as e:
            print(f"  ✗ Error extracting hotstart: {e}")
    else:
        print("  ⚠ No previous simulation folder found; skipping hotstart extraction")

    print("\n" + "="*80)
    print(f"  ✓ PREPROCESSING COMPLETE – All outputs saved to: {folder}")
    print("="*80 + "\n")
except Exception as e:
    print(f"  ✗ Error copying template files: {e}")
