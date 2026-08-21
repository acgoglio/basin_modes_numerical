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

# === Stack values, convert to hours if needed ===
mode_data_list = []
for var in modes_vars:
    arr = ds[var].data
    # Se di tipo timedelta converti in ore
    if np.issubdtype(arr.dtype, np.timedelta64):
        arr = arr.astype('timedelta64[s]').astype(float) / 3600
    else:
        arr = arr.astype(float)
    mode_data_list.append(arr)
mode_data = np.stack(mode_data_list, axis=-1)   # shape (y, x, K)
periods_in_hours = mode_data.flatten()

# ----------------------------------------------------------------------
# OPZIONE C: etichetta di provenienza di ogni valore.
# Il flatten() perde l'informazione su QUALE grid point ha prodotto ogni
# periodo. Se un punto contribuisce due picchi allo stesso gruppo (cosa
# possibile, perche' la tolleranza puo' contenere piu' bin dello spettro),
# quel punto veniva contato due volte e la colonna "%" poteva superare il
# 100%. Qui si conserva l'indice del punto per poter contare i punti
# DISTINTI in fase di raggruppamento.
# Il flatten() e' in ordine C su (y, x, K), quindi il valore in posizione j
# proviene dal grid point j // K.
# ----------------------------------------------------------------------
n_modes_stacked = mode_data.shape[-1]
point_id_all = np.arange(periods_in_hours.size) // n_modes_stacked

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
        # OPZIONE C: punti distinti, non numero di valori
        n_points = len(np.unique(point_id_all[group.index]))
        greedy_groups.append((round(group.mean(), 2), n_points, len(group)))
        remaining = remaining.drop(group.index)

    df_greedy = pd.DataFrame(greedy_groups, columns=["Grouped_Period", "Count", "N_values"])
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
        # OPZIONE C: punti distinti, non numero di valori
        n_points = len(np.unique(point_id_all[group.index]))
        greedy_groups.append((group_mean, n_points, len(group)))
        tolerances_list.append(tolerance)
        remaining = remaining.drop(group.index)

    df_greedy = pd.DataFrame(greedy_groups, columns=["Grouped_Period", "Count", "N_values"])
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
            width=2*df_greedy["Tolerance"], color="tab:blue", edgecolor="black")
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

# ======================================================================
# FIGURA COMBINATA (Fig. 3.5 della tesi)
#   pannello (a): periodi rilevati punto per punto, prima del raggruppamento
#   pannello (b): gruppi risultanti dalla procedura di raggruppamento
# Sostituisce la composizione manuale di due PNG separati.
# Scelte:
#  - asse y in % dei punti mare (il numero assoluto non e' informativo);
#  - pannello (a) con vlines: le barre a width=0.06 venivano coperte dal
#    proprio bordo nero e apparivano nere invece che colorate;
#  - pannello (b) con width = 2*Tolerance, perche' il raggruppamento usa
#    |T - mode| <= tolerance, cioe' una finestra larga 2*tolerance;
#  - taglio a TMIN_CUT: sotto quel periodo i risultati non sono mostrati;
#  - banda grigia TMIN_CUT..TMIN_SHADE: modi discussi in appendice.
# ======================================================================
TMIN_CUT   = 3.0    # taglio duro dell'asse x [h]
TMIN_SHADE = 6.0    # limite superiore della banda ombreggiata [h]
ANNOT_PCT  = 80.0   # annota i gruppi che superano questa copertura [%]
LOG_X      = True   # asse x logaritmico (False per lineare)
LOG_Y      = False  # asse y logaritmico (attenzione: su barre la lunghezza
                    # non e' piu' proporzionale al valore e dipende da YMIN_LOG)
YMIN_LOG   = 0.1    # fondo dell'asse y quando LOG_Y=True [%]
BAR_FULL   = True   # True: larghezza = 2*tol (intervallo effettivo di
                    # raggruppamento). False: larghezza = tol (risoluzione
                    # spettrale, criterio di Rayleigh)

from matplotlib.patches import Patch
from matplotlib.lines import Line2D

# copertura in percentuale dei punti mare
pct_all    = 100.0 * df_all["Count"].values / sea_gp_num
pct_greedy = 100.0 * df_greedy["Count"].values / sea_gp_num

sel_a = df_all["Period"].values >= TMIN_CUT
sel_b = df_greedy["Grouped_Period"].values >= TMIN_CUT

fig, (axa, axb) = plt.subplots(nrows=2, figsize=(8.5, 7), sharex=True)

for ax in (axa, axb):
    ax.axvspan(TMIN_CUT, TMIN_SHADE, color="0.85", alpha=0.6, zorder=0)
    ax.grid(axis="y", linestyle="--", alpha=0.4, zorder=0)
    if LOG_Y:
        ax.set_yscale("log")
        ax.set_ylim(YMIN_LOG, 1500)
        ax.set_yticks([YMIN_LOG, 1, 10, 100])
        ax.yaxis.set_major_formatter(
            mpl.ticker.FuncFormatter(lambda v, pos: f"{v:g}"))
    else:
        ax.set_ylim(0, 145)
        ax.set_yticks(np.arange(0, 101, 20))
    ax.set_ylabel("Coverage (% of sea grid points)")

# --- (a) periodi grezzi ---
axa.vlines(df_all["Period"].values[sel_a],
           (YMIN_LOG if LOG_Y else 0), pct_all[sel_a],
           color="tab:red", linewidth=0.9, zorder=3)
axa.text(0.012, 0.94, "(a) before grouping", transform=axa.transAxes,
         fontsize=10, va="top")
axa.legend(handles=[
    Line2D([0], [0], color="tab:red", lw=1.4,
           label="period detected at individual grid points"),
    Patch(facecolor="0.85", edgecolor="none",
          label=f"{TMIN_CUT:.0f}-{TMIN_SHADE:.0f} h: see Appendix")],
    loc="upper right", fontsize=7.5, framealpha=0.95, ncol=1,
    borderpad=0.4, handlelength=1.6)

# --- (b) periodi raggruppati ---
axb.bar(df_greedy["Grouped_Period"].values[sel_b], pct_greedy[sel_b],
        width=(2 if BAR_FULL else 1) * df_greedy["Tolerance"].values[sel_b],
        bottom=(YMIN_LOG if LOG_Y else 0),
        color="tab:blue", edgecolor="black", linewidth=0.4, zorder=3)
axb.text(0.012, 0.94, "(b) after grouping", transform=axb.transAxes,
         fontsize=10, va="top")
axb.set_xlabel("Period (hours)")
axb.legend(handles=[
    Patch(facecolor="tab:blue", edgecolor="black", linewidth=0.4,
          label=("grouped mode (bar width = grouping window)" if BAR_FULL
                 else "grouped mode (bar width = spectral resolution)")),
    Patch(facecolor="0.85", edgecolor="none",
          label=f"{TMIN_CUT:.0f}-{TMIN_SHADE:.0f} h: see Appendix")],
    loc="upper right", fontsize=7.5, framealpha=0.95, ncol=1,
    borderpad=0.4, handlelength=1.6)

# etichette sui gruppi con copertura > ANNOT_PCT
for T, p, tol in zip(df_greedy["Grouped_Period"].values[sel_b],
                     pct_greedy[sel_b],
                     df_greedy["Tolerance"].values[sel_b]):
    if p > ANNOT_PCT:
        prec = max(0, -int(np.floor(np.log10(abs(tol)))))
        axb.annotate(f"{T:.{prec}f}", xy=(T, p), xytext=(0, 4),
                     textcoords="offset points", ha="center",
                     fontsize=7, rotation=90)

if LOG_X:
    axb.set_xscale("log")
    axb.set_xlim(TMIN_CUT, th_filter + 2)
    axb.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(
        lambda v, pos: f"{v:g}"))
    axb.set_xticks([3, 4, 5, 6, 8, 10, 15, 20, 30, 40])
else:
    axb.set_xlim(TMIN_CUT, th_filter + 1)
    axb.set_xticks(np.arange(5, th_filter + 1, 5))

plt.tight_layout()
plt.savefig(os.path.join(indir, "fig_mode_periods_pow.png"), dpi=300)
plt.close()
print("Saved: fig_mode_periods_pow.png (combined 2-panel figure)")

