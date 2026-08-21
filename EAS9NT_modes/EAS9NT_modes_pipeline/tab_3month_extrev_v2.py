"""
Versione di tab_3month_extrev.py che carica i dati da
forced_modes_results.npz (prodotto da merge_forced_modes.py) invece di
averli scritti a mano - stessa identica logica di plotting
dell'originale, cambia solo la sorgente di 'data' ed 'events'.
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap
mpl.use("Agg")

# =========================
# Caricamento risultati calcolati (al posto dei valori hardcoded)
# =========================
RESULTS_FILE = "/work/cmcc/ag15419/basin_modes_new/basin_modes_forced/merged/forced_modes_results.npz"  # DA VERIFICARE: allinea a work_dir

npz = np.load(RESULTS_FILE, allow_pickle=True)

data = npz["data"]      # shape (10, 4, 4) modo x trimestre x anno - NON ancora mediato sugli anni
events = npz["events"]  # shape (10, 7) modo x evento

modes = [f"Mode {i}" for i in range(10)]
trimesters = list(npz["trimesters"])
events_names = list(npz["events_names"])

# QA: segnala se sono rimasti buchi (finestre non trovate/non processate)
if np.isnan(data).any():
    n_missing = np.isnan(data).sum()
    print(f"ATTENZIONE: {n_missing} valori mancanti (NaN) in 'data' (trimestrale) - "
          f"controlla i box/finestre corrispondenti prima di usare la figura.")
if np.isnan(events).any():
    n_missing = np.isnan(events).sum()
    print(f"ATTENZIONE: {n_missing} valori mancanti (NaN) in 'events' - "
          f"controlla i box/finestre corrispondenti prima di usare la figura.")

# Mean over years (axis 2) - IGNORANDO i NaN (invece di lasciare che un
# singolo anno mancante corrompa la media di tutto il trimestre)
data_mean = np.nanmean(data, axis=2)  # shape 10 x 4

# ----------------------
# Plot side-by-side: trimesters + extreme events
# (invariato rispetto all'originale)
# ----------------------
fig, axes = plt.subplots(1, 2, figsize=(16, 6), gridspec_kw={'width_ratios': [1, 0.7]})

cmap_trimesters = plt.get_cmap("cubehelix_r")

vmin = 0
vmax = 8
transition_value = 1.0

n_low = 128
n_high = 128

low_colors = plt.get_cmap("cubehelix_r")(np.linspace(0, 1, n_low))
high_colors = plt.get_cmap("gnuplot")(np.linspace(0, 1, n_high))
colors_combined = np.vstack([low_colors, high_colors])

transition_norm = (transition_value - vmin) / (vmax - vmin)
positions_low = np.linspace(0, transition_norm, n_low)
positions_high = np.linspace(transition_norm, 1.0, n_high)
positions = np.concatenate([positions_low, positions_high])

cmap_events = LinearSegmentedColormap.from_list("cmap_events", list(zip(positions, colors_combined)))
cmap_events = plt.get_cmap("gist_stern_r")
cmap_events = LinearSegmentedColormap.from_list(
    "gist_stern_r_no_black",
    cmap_events(np.linspace(0, 0.95, 256)))

sns.heatmap(data_mean, ax=axes[0], cmap=cmap_trimesters, annot=True, fmt=".1f",
            xticklabels=trimesters, yticklabels=modes, cbar_kws={'label': 'Mean amplitude (cm)'},
            vmin=0, vmax=transition_value)
axes[0].set_title("Mean seiches amplitude per trimester (2020-2023)")
axes[0].set_xlabel("Trimesters")
axes[0].set_ylabel("Modes")
axes[0].tick_params(axis='x', rotation=0)

sns.heatmap(events, ax=axes[1], cmap=cmap_events, annot=True, fmt=".1f",
            xticklabels=events_names, yticklabels=modes, cbar_kws={'label': '99th perc. amplitude (cm)'},
            vmin=vmin, vmax=vmax)
axes[1].set_title("Seiches amplitude (99th percentile) during extreme events")
axes[1].set_xlabel("Events")
axes[1].set_ylabel("Modes")
axes[1].tick_params(axis='x', rotation=45)
plt.setp(axes[1].get_xticklabels(), ha='right')

plt.tight_layout()
plt.savefig("heatmap_trimesters_events.png", dpi=300, bbox_inches='tight')
print("Salvato: heatmap_trimesters_events.png")
