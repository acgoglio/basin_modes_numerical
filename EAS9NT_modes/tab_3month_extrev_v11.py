"""
Versione di tab_3month_extrev.py che carica i dati da
forced_modes_results.npz.

v10: colorbar e celle disegnate COMPLETAMENTE A MANO (rettangoli
matplotlib.patches.Rectangle con colore calcolato esplicitamente in
numpy), niente BoundaryNorm/extend/Colorbar automatici - dopo tre
tentativi falliti con i meccanismi automatici di matplotlib (che su
matplotlib 3.3.3 si comportano diversamente dalle versioni recenti),
questo approccio non dipende da nessun comportamento automatico
version-specific: e' solo geometria e colori scelti da noi.
"""

import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import Rectangle
from matplotlib.colors import LinearSegmentedColormap
mpl.use("Agg")

import forced_windows_ini as ini

RESULTS_FILE = "/work/cmcc/ag15419/basin_modes_new/basin_modes_EAS9NT/merged/forced_modes_results.npz"

npz = np.load(RESULTS_FILE, allow_pickle=True)

data = npz["data"]
events = npz["events"]
annual = npz["annual"]

mode_areas = list(npz["modes"])
trimesters = list(npz["trimesters"])
events_names = list(npz["events_names"])
years = list(npz["years"])

if np.isnan(data).any():
    print(f"ATTENZIONE: {np.isnan(data).sum()} valori mancanti (NaN) in 'data'.")
if np.isnan(events).any():
    print(f"ATTENZIONE: {np.isnan(events).sum()} valori mancanti (NaN) in 'events'.")
if np.isnan(annual).any():
    print(f"ATTENZIONE: {np.isnan(annual).sum()} valori mancanti (NaN) in 'annual'.")

data_mean = np.nanmean(data, axis=2)

modes = [f"Mode {i}\n({area})" for i, area in enumerate(mode_areas)]

event_year_map = {name: start_date[:4] for name, start_date, _ in ini.EVENTS}
events_labels = []
for name in events_names:
    match = re.search(r"\d{4}", name)
    if match:
        events_labels.append(f"{name[:match.start()].strip()}\n{match.group()}")
    else:
        events_labels.append(f"{name}\n{event_year_map[name]}")

plt.rcParams.update({
    "font.size": 20,
    "axes.titlesize": 22,
    "axes.labelsize": 20,
})

# =====================================================================
# Colori scelti A MANO (nessun campionamento automatico da colormap
# continue, nessuna dipendenza da BoundaryNorm): 6 colori distinti per
# le bande 0-3cm + magenta per "oltre 3cm". Puoi cambiare questi colori
# semplicemente editando questa lista.
# =====================================================================
vmin, vmax = 0.0, 3.0
BIN_STEP = 0.5
boundaries = np.round(np.arange(vmin, vmax + BIN_STEP/2, BIN_STEP), 2)  # [0,0.5,1,...,3]
n_bins = len(boundaries) - 1

BAND_COLORS = [
    "#ffffff",  # 0.0 - 0.5
    "#d9d966",  # 0.5 - 1.0
    "#9e9e7d",  # 1.0 - 1.5
    "#6e6edc",  # 1.5 - 2.0
    "#153d7a",  # 2.0 - 2.5
    "#dc0c18",  # 2.5 - 3.0
]
OVER_COLOR = "#ff00ff"  # magenta, per valori oltre vmax
assert len(BAND_COLORS) == n_bins, f"Servono esattamente {n_bins} colori, trovati {len(BAND_COLORS)}"


def value_to_color(v):
    """Ritorna il colore (stringa hex) per un valore scalare, secondo
    le bande BAND_COLORS/boundaries definite sopra."""
    if np.isnan(v):
        return "#f0f0f0"  # grigio chiaro per i NaN, invece di un colore delle bande
    if v >= vmax:
        return OVER_COLOR
    idx = np.searchsorted(boundaries, v, side="right") - 1
    idx = max(0, min(idx, n_bins - 1))
    return BAND_COLORS[idx]


def draw_manual_panel(ax, matrix, xticklabels, title, xlabel, rotate_x=0, ha=None):
    """Disegna un pannello heatmap interamente a mano: un Rectangle
    colorato per cella + testo con il valore, niente pcolormesh."""
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
                ax.text(c + 0.5, n_rows - 1 - r + 0.5, f"{v:.1f}",
                        ha="center", va="center", fontsize=18, color=text_color)

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
    """Colorbar disegnata interamente a mano: un Rectangle per banda +
    un triangolo per l'estensione magenta, niente Colorbar/BoundaryNorm."""
    cax = fig.add_axes(cbar_ax_rect)
    for i, color in enumerate(BAND_COLORS):
        cax.add_patch(Rectangle((0, i), 1, 1, facecolor=color, edgecolor="black", linewidth=0.5))
    cax.add_patch(plt.Polygon([[0, n_bins], [1, n_bins], [0.5, n_bins + 0.6]],
                               facecolor=OVER_COLOR, edgecolor="black", linewidth=0.5))
    cax.set_xlim(0, 1)
    cax.set_ylim(0, n_bins + 0.6)
    cax.set_xticks([])
    cax.set_yticks(range(n_bins + 1))
    cax.set_yticklabels([f"{b:g}" for b in boundaries], fontsize=17)
    cax.yaxis.tick_right()
    cax.set_ylabel(label, fontsize=20)
    cax.yaxis.set_label_position("right")
    for spine in cax.spines.values():
        spine.set_visible(False)


# ---- Figura completa: annuale + trimestri + eventi (3 pannelli) ----
fig = plt.figure(figsize=(18, 24))
gs = fig.add_gridspec(3, 1, right=0.85)
ax0 = fig.add_subplot(gs[0])
ax1 = fig.add_subplot(gs[1])
ax2 = fig.add_subplot(gs[2])

draw_manual_panel(ax0, annual, years, "Mean seiches amplitude per year", "Years")
draw_manual_panel(ax1, data_mean, trimesters, "Mean seiches amplitude per trimester (2020-2023)", "Trimesters")
draw_manual_panel(ax2, events, events_labels, "Seiches amplitude during extreme sea level events", "Events",
                   rotate_x=45, ha="right")
draw_manual_colorbar(fig, [0.87, 0.35, 0.03, 0.35], "Seiches amplitude (cm)")

plt.savefig("heatmap_annual_trimester_events.png", dpi=300, bbox_inches="tight")
print("Salvato: heatmap_annual_trimester_events.png")
plt.close(fig)

# ---- Seconda figura: solo trimestri + eventi ----
fig2 = plt.figure(figsize=(18, 16))
gs2 = fig2.add_gridspec(2, 1, right=0.85)
ax20 = fig2.add_subplot(gs2[0])
ax21 = fig2.add_subplot(gs2[1])

draw_manual_panel(ax20, data_mean, trimesters, "Mean seiches amplitude per trimester (2020-2023)", "Trimesters")
draw_manual_panel(ax21, events, events_labels, "Seiches amplitude during extreme events", "Events",
                   rotate_x=45, ha="right")
draw_manual_colorbar(fig2, [0.87, 0.3, 0.03, 0.4], "Seiches amplitude (cm)")

plt.savefig("heatmap_trimester_events.png", dpi=300, bbox_inches="tight")
print("Salvato: heatmap_trimester_events.png")
plt.close(fig2)

