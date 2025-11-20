import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import os
import matplotlib as mpl
from area_ini import *
mpl.use("Agg")  # For non-interactive backend

# === Parameters ===
sea_gp_num = 145023      # Total sea grid points

# === Input directory and file ===
indir = work_dir
nc_file = os.path.join(indir, "basin_modes_pow_med.nc")
ds = xr.open_dataset(nc_file)

# === Extract mode period variables (m?_T) ===
modes_vars = [var for var in ds.data_vars if var.startswith("m") and "_T" in var]
print(f"Found {len(modes_vars)} mode period variables.")

# === Stack values, convert to hours ===
mode_data = np.stack([
    ds[var].data.astype("timedelta64[s]").astype("float") for var in modes_vars
], axis=-1)
periods_in_hours = mode_data.flatten() / 3600

# === Clean values ===
periods_series = pd.Series(periods_in_hours, name="Period")
periods_series = periods_series.dropna()
periods_series = periods_series[(periods_series > 0) & (periods_series <= th_filter)]

if periods_series.empty:
    print("No valid period values found between 0 and 40 hours.")
    exit()

# === Round and save all periods ===
rounded_periods = periods_series.round(2)
rounded_periods = rounded_periods[rounded_periods > 0] 
df_all = rounded_periods.value_counts().reset_index()
df_all.columns = ["Period", "Count"]
df_all["%"] = (df_all["Count"] / sea_gp_num * 100).round(2)
#df_all = df_all.sort_values("Count", ascending=False).reset_index(drop=True)
df_all = df_all.sort_values("Period", ascending=False).reset_index(drop=True)
print ('Order by period..')
df_all.to_csv(os.path.join(indir, "periods_all_pow.csv"), index=False)
print("Saved: periods_all_pow.csv")

# === Plot histogram and top 10 table ===
df_top10 = df_all.head(10)
table_data = list(zip(
    df_top10["Period"].round(2),
    df_top10["Count"],
    df_top10["%"].round(1)
))
column_labels = ["Period (h)", "Frequency (grid points)", "Percentage (%)"]

fig, (ax1, ax2) = plt.subplots(nrows=2, figsize=(10, 8), gridspec_kw={'height_ratios': [2.5, 1]})
ax1.bar(df_all["Period"], df_all["Count"], width=0.06, color="tab:orange", edgecolor="black")
ax1.set_xlabel("Period (hours)")
ax1.set_ylabel("Frequency (grid points)")
ax1.set_title("Frequency of Mode Periods in the Mediterranean Sea")
ax1.grid(axis="y", linestyle="--", alpha=0.6)
ax1.tick_params(axis='x', rotation=90)
ax1.set_ylim(0, 150000)

ax2.axis('off')
table = ax2.table(cellText=table_data,
                  colLabels=column_labels,
                  cellLoc='center',
                  loc='center')
table.auto_set_font_size(False)
table.set_fontsize(10)
table.scale(1.2, 1.3)
for (row, col), cell in table.get_celld().items():
    if row == 0:
        cell.set_facecolor('#D3D3D3')
        cell.set_text_props(weight='bold')

plt.tight_layout()
plt.savefig(os.path.join(indir, "hist_all_pow_table10_singletable.png"), dpi=300)
print("Saved: histogram with top 10 table (all values)")

# === Grouping algorithm ===
if flag_var_unc == 0:
    # Fixed tolerance grouping
    tolerance = fixed_uncertainty
    remaining = rounded_periods.copy()
    greedy_groups = []

    while not remaining.empty:
        mode = remaining.mode()[0]
        group = remaining[np.abs(remaining - mode) <= tolerance]
        greedy_groups.append((round(group.mean(), 2), len(group)))
        remaining = remaining.drop(group.index)

    df_greedy = pd.DataFrame(greedy_groups, columns=["Grouped_Period", "Count"])
    df_greedy["%"] = (df_greedy["Count"] / sea_gp_num * 100).round(2)
    #df_greedy = df_greedy.sort_values("Count", ascending=False).reset_index(drop=True)
    df_greedy = df_greedy.sort_values("Grouped_Period", ascending=False).reset_index(drop=True)
    print ('Order by period..')
    df_greedy.to_csv(os.path.join(indir, "periods_grouped_pow.csv"), index=False)
    print("Saved: periods_grouped_pow.csv (fixed tolerance)")

    # Plot
    df_top10 = df_greedy.head(10)
    table_data = list(zip(
        df_top10["Grouped_Period"].round(2),
        df_top10["Count"],
        df_top10["%"].round(1)
    ))
    column_labels = ["Grouped Period (h)", "Frequency (grid points)", "Percentage (%)"]

    fig, (ax1, ax2) = plt.subplots(nrows=2, figsize=(10, 8), gridspec_kw={'height_ratios': [2.5, 1]})
    ax1.bar(df_greedy["Grouped_Period"], df_greedy["Count"], width=fixed_uncertainty, color="tab:green", edgecolor="black")
    ax1.set_xlabel(f"Grouped Period (hours ±{tolerance}h)")
    ax1.set_ylabel("Frequency (grid points)")
    ax1.set_title("Grouped Mode Periods (Fixed Tolerance)")
    ax1.grid(axis="y", linestyle="--", alpha=0.6)
    ax1.tick_params(axis='x', rotation=45)
    ax1.set_ylim(0, 150000)

    ax2.axis('off')
    table = ax2.table(cellText=table_data,
                      colLabels=column_labels,
                      cellLoc='center',
                      loc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.2, 1.3)
    for (row, col), cell in table.get_celld().items():
        if row == 0:
            cell.set_facecolor('#D3D3D3')
            cell.set_text_props(weight='bold')

    plt.tight_layout()
    plt.savefig(os.path.join(indir, "hist_grouped_pow_table10_singletable.png"), dpi=300)
    print("Saved: grouped histogram and top 10 table")

elif flag_var_unc == 1:
    # Variable tolerance based on spectral resolution
    segment_len_hours = segment_len_days * 24
    delta_f = 1 / segment_len_hours  # Spectral resolution in cph
    remaining = rounded_periods.copy()
    greedy_groups = []
    tolerances_list = []

    def round_to_1_sigfig(x):
        if x == 0:
            return 0
        return round(x, -int(np.floor(np.log10(abs(x)))))

    while not remaining.empty:
        mode = remaining.mode()[0]
        tolerance = (mode ** 2) * delta_f  
        tolerance = tolerance + extra_unc * tolerance
        tolerance = round_to_1_sigfig(tolerance)
        if tolerance < min_unc :
           tolerance = min_unc
        print ('Tolerance:',tolerance)
        group = remaining[np.abs(remaining - mode) <= tolerance]
        if len(group) == 0:
            remaining = remaining.drop(remaining[remaining == mode].index)
            continue
        group_mean = round(group.mean(), -int(np.floor(np.log10(abs(tolerance)))))
        greedy_groups.append((group_mean, len(group)))
        tolerances_list.append(tolerance)
        remaining = remaining.drop(group.index)

    df_greedy = pd.DataFrame(greedy_groups, columns=["Grouped_Period", "Count"])
    df_greedy["Tolerance"] = tolerances_list
    df_greedy["Percentage"] = df_greedy["Count"] / sea_gp_num * 100
    #df_greedy = df_greedy.sort_values("Count", ascending=False).reset_index(drop=True)
    df_greedy = df_greedy.sort_values("Grouped_Period", ascending=False).reset_index(drop=True)
    print ('Order by period..')
    df_greedy.to_csv(os.path.join(indir, "periods_grouped_pow.csv"), index=False)
    print("Saved: periods_grouped_pow.csv (variable tolerance)")

    df_top10 = df_greedy.head(10)
    def format_group_period(row):
        prec = -int(np.floor(np.log10(abs(row["Tolerance"]))))
        return f"{row['Grouped_Period']:.{prec}f} ±{row['Tolerance']:.{prec}f} h"

    table_data = list(zip(
        df_top10.apply(format_group_period, axis=1),
        df_top10["Count"],
        df_top10["Percentage"].round(1)
    ))
    column_labels = ["Grouped Period (h ± res)", "Frequency (grid points)", "Percentage (%)"]

    fig, (ax1, ax2) = plt.subplots(nrows=2, figsize=(10, 8), gridspec_kw={'height_ratios': [2.5, 1]})
    ax1.bar(df_greedy["Grouped_Period"], df_greedy["Count"],
            width=df_greedy["Tolerance"], color="tab:blue", edgecolor="black")
    ax1.set_xlabel("Grouped Period (hours ± resolution)")
    ax1.set_ylabel("Frequency (grid points)")
    ax1.set_title("Grouped Mode Periods (Variable Tolerance)")
    ax1.grid(axis="y", linestyle="--", alpha=0.6)
    ax1.tick_params(axis='x', rotation=45)
    ax1.set_ylim(0, 150000)

    ax2.axis('off')
    table = ax2.table(cellText=table_data,
                      colLabels=column_labels,
                      cellLoc='center',
                      loc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.2, 1.3)
    for (row, col), cell in table.get_celld().items():
        if row == 0:
            cell.set_facecolor('#D3D3D3')
            cell.set_text_props(weight='bold')

    plt.tight_layout()
    plt.savefig(os.path.join(indir, "hist_grouped_pow_table10_singletable.png"), dpi=300)
    print("Saved: grouped histogram and top 10 table (variable tolerance)")
