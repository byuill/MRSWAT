"""
Plot saltwater intrusion model results.

Loads model output (toe location), river discharge boundary condition,
and downstream stage boundary condition. Produces an 8.5 x 11 inch
figure with two panels:
    1. Time series of intrusion distance, TMDL threshold, and discharge.
         (x-axis shown as calendar dates computed from a user-provided
         reference date)
    2. Scatter of discharge vs intrusion distance coloured by downstream
         stage relative to its temporal mean.

This version asks for a reference date (YYYY-MM-DD) via `--ref-date`
or interactive input, computes calendar dates for the provided time
series (days since reference), identifies future times relative to
the current system date as forecast, shades the forecast portion in
the top plot, and draws a dashed sill line at intrusion distance = 64
on both plots.
"""

import datetime as dt
import numpy as np
import matplotlib
import os
if os.environ.get("MPLBACKEND") is None:
        # prefer non-interactive when running in batch; allow override
        matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.dates as mdates

# ── working directory paths ──────────────────────────────────────────
toe_file   = "mr-toe-location-output.txt"
flow_file  = "Discharges.txt"
stage_file = "Stages.txt"

# ── load data ────────────────────────────────────────────────────────
# Toe-location output (skip header row)
toe_data = np.loadtxt(toe_file, skiprows=1)
toe_time = toe_data[:, 0]          # days since reference
toe_loc  = toe_data[:, 1]          # river mile  (Toe-Location)
tmdl_loc = toe_data[:, 2]          # river mile  (0.25-ppt-at-surface)

# River discharge (skip header row)
flow_data = np.loadtxt(flow_file, skiprows=1)
flow_time = flow_data[:, 0]        # days since reference
flow_q    = flow_data[:, 1]        # cfs

# Downstream stage (skip header row; ignore incomplete rows)
with open(stage_file, "r") as _f:
    _stage_rows = []
    for _line in _f:
        _parts = _line.split()
        if len(_parts) == 2:
            try:
                _stage_rows.append([float(_parts[0]), float(_parts[1])])
            except ValueError:
                pass   # skip header / non-numeric lines
stage_data = np.array(_stage_rows)
stage_time = stage_data[:, 0]      # days since reference
stage_wsel = stage_data[:, 1]      # ft

# ── compute temporally-averaged downstream stage ─────────────────────
stage_mean = np.mean(stage_wsel)

# ---------------- parse reference and simulation dates -----------------
date_info_file = "Date_Information.txt"
try:
    date_info = np.loadtxt(date_info_file, dtype=str)
except Exception as e:
    raise SystemExit(f"Cannot read {date_info_file}: {e}")

if date_info.shape[0] < 2 or date_info.shape[1] < 2:
    raise SystemExit(f"{date_info_file} must have at least two rows and two columns.")

ref_str = date_info[1, 0].strip()
sim_str = date_info[1, 1].strip()

try:
    ref_date = dt.datetime.strptime(ref_str, "%Y%m%d")
    sim_date = dt.datetime.strptime(sim_str, "%Y%m%d").date()
except Exception as e:
    raise SystemExit(f"Cannot parse dates in {date_info_file}: {e}")

today = dt.datetime.now().date()
if sim_date != today:
    offset_days = (today - sim_date).days
    print(
        f"WARNING: simulation is {abs(offset_days)} days before the current date "
        f"and {offset_days} is the offset in days."
    )

# convert times (days) -> datetimes
toe_dates = np.array([ref_date + dt.timedelta(days=float(d)) for d in toe_time])
flow_dates = np.array([ref_date + dt.timedelta(days=float(d)) for d in flow_time])
stage_dates = np.array([ref_date + dt.timedelta(days=float(d)) for d in stage_time])

# determine 'now' and mark forecast points as dates strictly after now
now = dt.datetime.now()
forecast_mask = toe_dates > now
forecast_start = None
if np.any(forecast_mask):
    forecast_start = toe_dates[np.argmax(forecast_mask)]

# ── build figure (8.5 x 11 in, portrait) ────────────────────────────
# Matplotlib figure
fig = plt.figure(figsize=(8.5, 11))

# Layout constants (inches → figure-fraction)
page_w, page_h = 8.5, 11.0
margin = 0.25                       # inch margin on all sides
plot_h = 2.0                        # each plot is 2 inches tall
left   = 1.00 / page_w             # left margin for y-label room
right  = 0.90 / page_w             # extra right margin for 2nd y-axis
gap    = 0.60                       # vertical gap between plots (inches)

ax_w = 1.0 - left - right          # axes width in figure fraction
ax_h = plot_h / page_h             # axes height in figure fraction

# --- first (top) axes: starts at top with 0.25" margin ---------------
top1 = 1.0 - margin / page_h - ax_h
ax1  = fig.add_axes([left, top1, ax_w, ax_h])

# --- second axes: below the first with a gap --------------------------
top2 = top1 - (gap + 0.5) / page_h - ax_h
ax2  = fig.add_axes([left, top2, ax_w, ax_h])

# =====================================================================
# PLOT 1 – Time series
# =====================================================================
color_toe  = "#1f77b4"
color_tmdl = "#d62728"
color_flow = "#2ca02c"

ln1 = ax1.plot(toe_dates, toe_loc,  color=color_toe,  linewidth=1.2,
               label="intrusion distance")
ln2 = ax1.plot(toe_dates, tmdl_loc, color=color_tmdl, linewidth=1.2,
               label="TMDL threshold")

ax1.set_xlabel("Date", fontsize=10)
ax1.set_ylabel("River Mile", fontsize=10, color="black")
ax1.tick_params(axis="y", labelcolor="black")

# Secondary y-axis for discharge
ax1r = ax1.twinx()
ln3 = ax1r.plot(flow_dates, flow_q / 1000.0, color=color_flow,
                linewidth=1.0, linestyle="--", alpha=0.8,
                label="Discharge")
ax1r.set_ylabel("Discharge (×1 000 cfs)", fontsize=10, color=color_flow)
ax1r.tick_params(axis="y", labelcolor=color_flow)

lns = ln1 + ln2 + ln3
labs = [l.get_label() for l in lns]
ax1.legend(lns, labs, loc="lower right", fontsize=8, framealpha=0.9)
ax1.set_title("Saltwater Intrusion – Time Series", fontsize=12,
              fontweight="bold")
ax1.grid(True, linewidth=0.3, alpha=0.6)

# format x-axis as dates
ax1.xaxis_date()
ax1.xaxis.set_major_locator(mdates.AutoDateLocator())
ax1.xaxis.set_major_formatter(mdates.DateFormatter('%b-%d'))
fig.autofmt_xdate(rotation=30)

# shade forecast region (light blue) on the top plot only
from matplotlib.patches import Patch
if forecast_start is not None:
    ax1.axvspan(mdates.date2num(forecast_start),
                mdates.date2num(toe_dates[-1]),
                color='lightblue', alpha=0.40, zorder=0,
                label='_nolegend_')
    # add 'Forecast' label to the existing legend
    forecast_patch = Patch(facecolor='lightblue', alpha=0.60, label='Forecast')
    lns_with_fc = lns + [forecast_patch]
    labs_with_fc = [l.get_label() for l in lns] + ['Forecast']
    ax1.legend(lns_with_fc, labs_with_fc, loc='lower right', fontsize=8, framealpha=0.9)

# draw dashed sill line at intrusion distance = 64 on top plot
ax1.axhline(64.0, color='k', linestyle='--', linewidth=1.0)

# =====================================================================
# PLOT 2 – Discharge vs Intrusion Distance (scatter)
# =====================================================================
# For each toe-location time step find the closest-in-time discharge
# and downstream stage values. Use numeric day offsets for interpolation.
q_matched     = np.interp(toe_time, flow_time, flow_q)
stage_matched = np.interp(toe_time, stage_time, stage_wsel)

# Colour: red if stage > mean, blue if stage <= mean
colors = np.where(stage_matched > stage_mean, "#FF0000", "#0000FF")

base_s = 18
ax2.scatter(q_matched / 1000.0, toe_loc, c=colors, s=base_s, alpha=0.8,
            edgecolors="none")

# highlight current date point (open black circle) and last date (open blue)
# find index nearest to 'now' and last index
if len(toe_dates) > 0:
    # compute seconds difference safely
    diffs = np.array([abs((d - now).total_seconds()) for d in toe_dates])
    idx_now = int(np.argmin(diffs))
    idx_last = len(toe_dates) - 1

    x_now = q_matched[idx_now] / 1000.0
    y_now = toe_loc[idx_now]
    x_last = q_matched[idx_last] / 1000.0
    y_last = toe_loc[idx_last]

    highlight_s = base_s * 2
    ax2.scatter([x_now], [y_now], facecolors='none', edgecolors='k', s=highlight_s,
                linewidths=1.2, zorder=5, label='current date')
    ax2.scatter([x_last], [y_last], facecolors='none', edgecolors='b', s=highlight_s,
                linewidths=1.2, zorder=5, label='last date')

ax2.set_xlabel("Discharge (×1 000 cfs)", fontsize=10)
ax2.set_ylabel("Intrusion Distance (River Mile)", fontsize=10)
ax2.set_title("Discharge vs Intrusion Distance", fontsize=12,
              fontweight="bold")
# draw dashed sill line at intrusion distance = 64 on scatter plot
ax2.axhline(64.0, color='k', linestyle='--', linewidth=1.0)
ax2.grid(True, linewidth=0.3, alpha=0.6)

# Legend for colour meaning
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], marker="o", color="w", markerfacecolor="#FF0000",
           markersize=8, label=f"Stage > mean ({stage_mean:+.2f} ft)"),
    Line2D([0], [0], marker="o", color="w", markerfacecolor="#0000FF",
           markersize=8, label=f"Stage ≤ mean ({stage_mean:+.2f} ft)"),
]
ax2.legend(handles=legend_elements, loc="upper right", fontsize=8,
           framealpha=0.9)

# ── 10-day and 28-day forecast values ───────────────────────────────
def _forecast_values(days_ahead):
    """Return (toe_loc, tmdl_loc) at now + days_ahead, or (None, None)."""
    target = now + dt.timedelta(days=days_ahead)
    # check that target falls within the result time series
    if target < toe_dates[0] or target > toe_dates[-1]:
        return None, None
    # interpolate using numeric day offsets
    target_day = (target - ref_date).total_seconds() / 86400.0
    t_interp = np.interp(target_day, toe_time, toe_loc)
    tmdl_interp = np.interp(target_day, toe_time, tmdl_loc)
    return t_interp, tmdl_interp

toe_10,  tmdl_10  = _forecast_values(10)
toe_28,  tmdl_28  = _forecast_values(28)

# Build text lines
lines_10 = (
    f"  10-day forecast  – Toe Location: {toe_10:.2f} RM  |  0.25-ppt Location: {tmdl_10:.2f} RM"
    if toe_10 is not None
    else "  10-day forecast  – not available in the current result time series."
)
lines_28 = (
    f"  28-day forecast  – Toe Location: {toe_28:.2f} RM  |  0.25-ppt Location: {tmdl_28:.2f} RM"
    if toe_28 is not None
    else "  28-day forecast  – not available in the current result time series."
)

forecast_text = "Forecast summary:\n" + lines_10 + "\n" + lines_28

# Print to terminal
print(forecast_text)

# Add as text annotation below ax2 in the figure
fig.text(
    0.5,
    top2 - 0.55 / page_h,          # just below the scatter plot
    forecast_text,
    ha='center', va='top',
    fontsize=9, family='monospace',
    bbox=dict(boxstyle='round,pad=0.4', facecolor='#f0f0f0', edgecolor='gray', linewidth=0.8)
)

# ── save & show ──────────────────────────────────────────────────────
plt.savefig("model_results_figure.png", dpi=300, bbox_inches='tight')
plt.show()
print("Figure saved to model_results_figure.png")
