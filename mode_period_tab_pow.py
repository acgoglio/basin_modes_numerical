import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import os
import matplotlib as mpl
mpl.use("Agg")  # For non-interactive backend

# Tot sea grid points in the Med 
sea_gp_num=int(145023)

# Load dataset
indir = "/work/cmcc/ag15419/basin_modes_20not/"
ds = xr.open_dataset(os.path.join(indir, "basin_modes_pow_med.nc"))

# Extract m?_T variables
modes_vars = [var for var in ds.data_vars if var.startswith("m") and "_T" in var]
mode_data = np.stack([ds[var].values for var in modes_vars], axis=-1)
flattened_data = mode_data.flatten()

# Convert to hours
periods_in_hours = pd.to_timedelta(flattened_data).total_seconds() / 3600

# Drop NaN and invalid values
periods_in_hours = pd.to_numeric(periods_in_hours, errors="coerce")
periods_in_hours = periods_in_hours.dropna()
periods_in_hours = periods_in_hours[periods_in_hours > 0]
periods_in_hours = periods_in_hours[periods_in_hours <= 40]

if periods_in_hours.empty:
    print("Nessun valore valido di periodo tra 0 e 40 ore trovato.")
    exit()

# Round periods to second digit
rounded_periods = pd.Series(np.round(periods_in_hours, 2), name="Period")

# Save full list (ordered by frequency)
df_all = rounded_periods.value_counts().reset_index()
df_all.columns = ["Period", "Count"]
df_all["%"] = (df_all["Count"] / sea_gp_num * 100).round(2)
df_all = df_all.sort_values("Count", ascending=False).reset_index(drop=True)
df_all.to_csv(os.path.join(indir, "periods_all_pow.csv"), index=False)

print("Totale valori prima del filtro:", len(flattened_data))
print("Dopo rimozione NaN:", periods_in_hours.shape[0])
print("Valori min/max:", periods_in_hours.min(), periods_in_hours.max())

# --- Histogram of all periods with table (Top 10) ---
# --- Top 10 tab ---
df_top10 = df_all.head(10)

# --- Tab ---
table_data = list(zip(
    df_top10["Period"].round(2),
    df_top10["Count"],
    df_top10["%"].round(1)
))
column_labels = ["Period (h)", "Frequency (grid points)", "Percentage (%)"]

# --- Fig con 2 subplots ---
fig, (ax1, ax2) = plt.subplots(nrows=2, figsize=(10, 8), gridspec_kw={'height_ratios': [2.5, 1]})

# --- Subplot 1: Histogram ---
ax1.bar(df_all["Period"], df_all["Count"],
        width=0.06, color="tab:orange", edgecolor="black")

ax1.set_xlabel("Period (hours)")
ax1.set_ylabel("Frequency (grid points)")
ax1.set_title("Frequency of Mode Periods in the Mediterranean Sea")
ax1.grid(axis="y", linestyle="--", alpha=0.6)
ax1.tick_params(axis='x', rotation=90)
ax1.set_ylim(0,150000)

# --- Subplot 2: Table ---
ax2.axis('off')
#ax2.text(0.5, 1.0, "Top 10 Main Modes", fontsize=12, ha='center', transform=ax2.transAxes)
table = ax2.table(cellText=table_data,
                  colLabels=column_labels,
                  cellLoc='center',
                  loc='center')

table.auto_set_font_size(False)
table.set_fontsize(10)
table.scale(1.2, 1.3)

# Tab header grigio e bold
for (row, col), cell in table.get_celld().items():
    if row == 0:
        cell.set_facecolor('#D3D3D3')
        cell.set_text_props(weight='bold')

# --- Save ---
plt.tight_layout()
plt.savefig(os.path.join(indir, "hist_all_pow_table10_singletable.png"), dpi=300)

############################
# --- Grouping algorithm ---
tolerance = 0.4
remaining = rounded_periods.copy()
greedy_groups = []

while not remaining.empty:
    mode = remaining.mode()[0]
    group = remaining[np.abs(remaining - mode) <= tolerance]
    greedy_groups.append((round(group.mean(), 2), len(group)))
    remaining = remaining.drop(group.index)

df_greedy = pd.DataFrame(greedy_groups, columns=["Grouped_Period", "Count"])
df_greedy["%"] = (df_greedy["Count"] / sea_gp_num * 100).round(2)
df_greedy = df_greedy.sort_values("Count", ascending=False).reset_index(drop=True)
df_greedy.to_csv(os.path.join(indir, "periods_grouped_pow.csv"), index=False)
df_greedy["Percentage"] = df_greedy["Count"] / sea_gp_num * 100

# --- Histogram and table of grouped periods ---
# --- Top 10 tab ---
df_top10 = df_greedy.nlargest(10, "Count")

# --- Tab ---
table_data = list(zip(
    df_top10["Grouped_Period"].round(2),
    df_top10["Count"],
    df_top10["Percentage"].round(1)
))
column_labels = ["Grouped Period (h)", "Frequency (grid points)", "Percentage (%)"]

# --- Fig 2 subplots ---
fig, (ax1, ax2) = plt.subplots(nrows=2, figsize=(10, 8), gridspec_kw={'height_ratios': [2.5, 1]})

# --- Subplot 1: Histogram ---
ax1.bar(df_greedy["Grouped_Period"], df_greedy["Count"],
        width=0.4, color="tab:green", edgecolor="black")

ax1.set_xlabel(f"Grouped Period (hours ±{tolerance}h)")
ax1.set_ylabel("Frequency (grid points)")
ax1.set_title("Frequency of Mode Periods in the Mediterranean Sea (Grouped)")
ax1.grid(axis="y", linestyle="--", alpha=0.6)
ax1.tick_params(axis='x', rotation=45)
ax1.set_ylim(0,150000)

# --- Subplot 2: Table ---
ax2.axis('off') 
#ax2.text(0.5, 1.0, "Top 10 Main Modes (Grouped)", fontsize=12, ha='center', transform=ax2.transAxes)
table = ax2.table(cellText=table_data,
                  colLabels=column_labels,
                  cellLoc='center',
                  loc='center')

table.auto_set_font_size(False)
table.set_fontsize(10)
table.scale(1.2, 1.3)  

# Tab header
for (row, col), cell in table.get_celld().items():
    if row == 0:
        cell.set_facecolor('#D3D3D3')  # gray
        cell.set_text_props(weight='bold')
# --- Save ---
plt.tight_layout()
plt.savefig(os.path.join(indir, "hist_grouped_pow_table10_singletable.png"), dpi=300)
