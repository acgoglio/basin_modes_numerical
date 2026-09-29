"""
Plot extra (file separato): differenze tra gli eventi (pannello c) e i
riferimenti annuale/trimestrale (pannelli a/b).

Pannello 1: events[m,e] - annual[m, anno_dell_evento] (stesso modo,
            confrontato con l'anno in cui l'evento e' realmente
            accaduto, non con la media dei 4 anni)
Pannello 2: events[m,e] - data_mean[m, trimestre_corrispondente]
            (trimestre determinato dal mese di inizio dell'evento,
            preso da forced_windows_ini.EVENTS - es. Storm Gloria,
            17 gennaio -> JFM)

Stesso approccio "disegnato a mano" di tab_3month_extrev_v10.py (niente
BoundaryNorm/Colorbar automatici), ma con una scala DIVERGENTE
(le differenze possono essere negative) invece che a magnitudine.
"""

import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import Rectangle
mpl.use("Agg")

import forced_windows_ini as ini

RESULTS_FILE = "/work/cmcc/ag15419/basin_modes_new/basin_modes_EAS9NT/merged/forced_modes_results.npz"

npz = np.load(RESULTS_FILE, allow_pickle=True)

data = npz["data"]      # (10, 4, 4) modo x trimestre x anno
events = npz["events"]  # (10, 7) modo x evento
annual = npz["annual"]  # (10, 4) modo x anno

mode_areas = list(npz["modes"])
trimesters = list(npz["trimesters"])
events_names = list(npz["events_names"])
years = list(npz["years"])

modes = [f"Mode {i}\n({area})" for i, area in enumerate(mode_areas)]

data_mean = np.nanmean(data, axis=2)      # (10,4) modo x trimestre

# ---- Etichette eventi: nome + anno a capo (come nel plot principale) ----
event_year_map = {name: start_date[:4] for name, start_date, _ in ini.EVENTS}
events_labels = []
for name in events_names:
    match = re.search(r"\d{4}", name)
    if match:
        events_labels.append(f"{name[:match.start()].strip()}\n{match.group()}")
    else:
        events_labels.append(f"{name}\n{event_year_map[name]}")

year_index = {y: i for i, y in enumerate(years)}

# ---- Trimestre corrispondente a ciascun evento, dal mese di inizio ----
def month_to_trimester(month):
    if month in (1, 2, 3):
        return "JFM"
    elif month in (4, 5, 6):
        return "AMJ"
    elif month in (7, 8, 9):
        return "JAS"
    else:
        return "OND"

event_trimester_map = {name: month_to_trimester(int(start_date[4:6]))
                        for name, start_date, _ in ini.EVENTS}
trimester_index = {t: i for i, t in enumerate(trimesters)}

print("Trimestre assegnato a ciascun evento:")
for name in events_names:
    print(f"  {name}: {event_trimester_map[name]}")

# ---- Calcolo delle differenze ----
n_modes = events.shape[0]
n_events = events.shape[1]

diff_annual = np.full((n_modes, n_events), np.nan)
diff_trimester = np.full((n_modes, n_events), np.nan)

for e, name in enumerate(events_names):
    y_idx = year_index[int(event_year_map[name])]
    diff_annual[:, e] = events[:, e] - annual[:, y_idx]
    t_idx = trimester_index[event_trimester_map[name]]
    diff_trimester[:, e] = events[:, e] - data_mean[:, t_idx]

ref_annual = np.full((n_modes, n_events), np.nan)
ref_trimester = np.full((n_modes, n_events), np.nan)
for e, name in enumerate(events_names):
    y_idx = year_index[int(event_year_map[name])]
    ref_annual[:, e] = annual[:, y_idx]
    t_idx = trimester_index[event_trimester_map[name]]
    ref_trimester[:, e] = data_mean[:, t_idx]

plt.rcParams.update({
    "font.size": 20,
    "axes.titlesize": 22,
    "axes.labelsize": 20,
})

# =====================================================================
# Colori scelti a mano - scala DIVERGENTE, simmetrica attorno allo 0.
# Cambia VMAX_DIFF se le tue differenze reali sono molto piu' grandi o
# piccole di +/-2 cm (controlla il min/max stampato sotto).
# =====================================================================
VMAX_DIFF = 2.0
BIN_STEP = 0.5
boundaries = np.round(np.arange(-VMAX_DIFF, VMAX_DIFF + BIN_STEP/2, BIN_STEP), 2)
n_bins = len(boundaries) - 1

BAND_COLORS = [
    "#08306b",  # -2.0 - -1.5
    "#4292c6",  # -1.5 - -1.0
    "#c6dbef",  # -1.0 - -0.5
    "#f7f7f7",  # -0.5 -  0.0
    "#fee0d2",  #  0.0 -  0.5
    "#fb6a4a",  #  0.5 -  1.0
    "#cb181d",  #  1.0 -  1.5
    "#67000d",  #  1.5 -  2.0
]
UNDER_COLOR = "#00e5ff"  # ciano, per differenze < -VMAX_DIFF
OVER_COLOR = "#ff00ff"   # magenta, per differenze > +VMAX_DIFF
assert len(BAND_COLORS) == n_bins, f"Servono esattamente {n_bins} colori, trovati {len(BAND_COLORS)}"

print(f"\nRange reale differenze annuale: {np.nanmin(diff_annual):.2f} - {np.nanmax(diff_annual):.2f} cm")
print(f"Range reale differenze trimestre: {np.nanmin(diff_trimester):.2f} - {np.nanmax(diff_trimester):.2f} cm")
print(f"(VMAX_DIFF attuale: +/-{VMAX_DIFF} cm - alza/abbassa se molti valori finiscono in ciano/magenta)")


def value_to_color(v):
    if np.isnan(v):
        return "#d0d0d0"
    if v < -VMAX_DIFF:
        return UNDER_COLOR
    if v >= VMAX_DIFF:
        return OVER_COLOR
    idx = np.searchsorted(boundaries, v, side="right") - 1
    idx = max(0, min(idx, n_bins - 1))
    return BAND_COLORS[idx]


def draw_manual_panel(ax, matrix, ref_matrix, xticklabels, title, xlabel, rotate_x=0, ha=None):
    n_rows, n_cols = matrix.shape
    for r in range(n_rows):
        for c in range(n_cols):
            v = matrix[r, c]
            color = value_to_color(v)
            ax.add_patch(Rectangle((c, n_rows - 1 - r), 1, 1, facecolor=color,
                                    edgecolor="white", linewidth=0.5))
            if not np.isnan(v):
                rgb = mpl.colors.to_rgb(color)
                luminance = 0.299 * rgb[0] + 0.587 * rgb[1] + 0.114 * rgb[2]
                text_color = "black" if luminance > 0.6 else "white"
                sign = "+" if v > 0 else ""
                ref = ref_matrix[r, c]
                if not np.isnan(ref) and ref != 0:
                    pct = v / ref * 100
                    psign = "+" if pct > 0 else ""
                    label = f"{sign}{v:.1f}\n({psign}{pct:.0f}%)"
                else:
                    label = f"{sign}{v:.1f}\n(n/a)"
                ax.text(c + 0.5, n_rows - 1 - r + 0.5, label,
                        ha="center", va="center", fontsize=15, color=text_color)

    ax.set_xlim(0, n_cols)
    ax.set_ylim(0, n_rows)
    ax.set_xticks(np.arange(n_cols) + 0.5)
    ax.set_xticklabels(xticklabels, rotation=rotate_x, ha=ha if ha else "center", fontsize=19)
    ax.set_yticks(np.arange(n_rows) + 0.5)
    ax.set_yticklabels(modes[::-1], fontsize=17)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Modes")
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)


def draw_manual_colorbar(fig, cbar_ax_rect, label):
    cax = fig.add_axes(cbar_ax_rect)
    # Estensione sotto (ciano)
    cax.add_patch(plt.Polygon([[0, 0], [1, 0], [0.5, -0.6]],
                               facecolor=UNDER_COLOR, edgecolor="black", linewidth=0.5))
    for i, color in enumerate(BAND_COLORS):
        cax.add_patch(Rectangle((0, i), 1, 1, facecolor=color, edgecolor="black", linewidth=0.5))
    # Estensione sopra (magenta)
    cax.add_patch(plt.Polygon([[0, n_bins], [1, n_bins], [0.5, n_bins + 0.6]],
                               facecolor=OVER_COLOR, edgecolor="black", linewidth=0.5))
    cax.set_xlim(0, 1)
    cax.set_ylim(-0.6, n_bins + 0.6)
    cax.set_xticks([])
    cax.set_yticks(range(n_bins + 1))
    cax.set_yticklabels([f"{b:+g}" if b != 0 else "0" for b in boundaries], fontsize=17)
    cax.yaxis.tick_right()
    cax.set_ylabel(label, fontsize=20)
    cax.yaxis.set_label_position("right")
    for spine in cax.spines.values():
        spine.set_visible(False)


# ---- Figura: differenze eventi vs annuale/trimestre ----
fig = plt.figure(figsize=(18, 22))
gs = fig.add_gridspec(2, 1, right=0.85)
ax0 = fig.add_subplot(gs[0])
ax1 = fig.add_subplot(gs[1])

draw_manual_panel(ax0, diff_annual, ref_annual, events_labels,
                   "Events minus corresponding annual amplitude", "Events",
                   rotate_x=45, ha="right")
draw_manual_panel(ax1, diff_trimester, ref_trimester, events_labels,
                   "Events minus corresponding trimester amplitude", "Events",
                   rotate_x=45, ha="right")
draw_manual_colorbar(fig, [0.87, 0.3, 0.03, 0.4], "Difference (cm)")

plt.savefig("heatmap_events_vs_annual_trimester_diff.png", dpi=300, bbox_inches="tight")
print("\nSalvato: heatmap_events_vs_annual_trimester_diff.png")
plt.close(fig)
