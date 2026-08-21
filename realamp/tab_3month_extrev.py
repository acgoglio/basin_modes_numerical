import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap
mpl.use("Agg")  # For non-interactive backend

# =========================
# Modi e trimestri
# =========================
modes = [f"Mode {i}" for i in range(10)]
trimesters = ["JFM","AMJ","JAS","OND"]
events_names = [
    "Storm Gloria",
    "Medicane Ianos",
    "Medicane Apollo",
    "Storm Blas",
    "Acqua Alta 2022",
    "Cyclone Helios",
    "Storm Daniel"
]

# ----------------------
# Trimestral data: modes x trimesters x years
# ----------------------
data = np.zeros((10,4,4))

# JFM 2020-2023
data[:,0,0] = [0.67,1.31,0.19,2.64,0.47,0.40,0.52,0.00,0.00,0.00]
data[:,0,1] = [0.32,0.63,0.12,1.28,0.06,0.12,0.12,0.32,0.27,0.28]
data[:,0,2] = [0.21,0.65,0.16,1.12,0.10,0.08,0.11,0.48,0.54,0.59]
data[:,0,3] = [0.35,0.83,0.18,1.33,0.24,0.13,0.12,0.12,0.20,0.64]

# AMJ 2020-2023
data[:,1,0] = [0.37,0.28,0.24,1.57,0.37,0.42,0.42,0.00,0.00,0.00]
data[:,1,1] = [0.21,0.24,0.09,0.87,0.16,0.05,0.09,0.33,0.21,0.33]
data[:,1,2] = [0.17,0.44,0.46,0.93,0.11,0.14,0.40,0.43,0.00,0.44]
data[:,1,3] = [0.23,0.12,0.07,0.74,0.18,0.04,0.13,0.13,0.17,0.17]


# JAS 2020-2023
data[:,2,0] = [0.14,0.47,0.22,2.01,0.15,0.23,0.23,0.00,0.00,0.00]
data[:,2,1] = [0.48,0.17,0.55,1.66,0.10,0.14,0.07,0.20,0.00,0.21]
data[:,2,2] = [0.32,0.34,0.07,1.41,0.17,0.22,0.21,0.26,0.00,0.34]
data[:,2,3] = [0.23,0.53,0.00,1.14,0.18,0.07,0.11,0.46,0.31,0.36]

# OND 2020-2023
data[:,3,0] = [0.61,2.35,0.24,0.68,0.62,1.24,2.23,0.00,0.00,0.00]
data[:,3,1] = [0.23,0.33,0.04,0.79,0.12,0.02,0.21,0.23,0.23,0.22]
data[:,3,2] = [0.12,0.42,0.08,0.58,0.08,0.08,0.07,0.30,0.41,0.38]
data[:,3,3] = [0.17,0.64,0.08,0.88,0.21,0.14,0.07,0.14,0.12,0.54]

# Mean over years (axis 2)
data_mean = data.mean(axis=2)  # shape 10 x 4

# ----------------------
# Extreme events matrix: modes x 7 events
# ----------------------
events = np.full((10,7), np.nan)
events[:,0] = [2.40,0.98,1.60,9.72,1.47,0.51,0.89,1.32,1.75,1.75]  # Storm Gloria
events[:,1] = [1.92,1.80,0.55,8.34,1.90,0.71,1.08,1.48,1.51,1.59]  # Medicane Ianos
events[:,2] = [1.15,2.08,0.94,5.41,0.55,0.69,0.69,1.02,2.05,1.52]  # Medicane Apollo
events[:,3] = [1.36,3.28,0.35,4.90,3.94,0.92,1.01,2.57,2.57,2.57]  # Storm Blas
events[:,4] = [11.93,15.91,1.84,8.81,1.74,1.81,1.25,2.64,15.91,10.37]  # Acqua Alta 2022
events[:,5] = [10.70,5.32,1.20,4.72,2.42,1.36,1.32,1.37,1.57,1.94] # Cyclone Helios 
events[:,6] = [7.03,2.44,0.60,7.01,2.83,1.31,3.25,3.25,3.25,3.25] # Storm Daniel


# Round to millimeter
#data_mean = np.round(data_mean, 1)
#events = np.round(events, 1)

# ----------------------
# Plot side-by-side: trimesters + extreme events
# ----------------------
fig, axes = plt.subplots(1,2, figsize=(16,6), gridspec_kw={'width_ratios':[1,0.7]})

# ----------------------
# Heatmap 1: Trimesters 0-2
# ----------------------
cmap_trimesters = plt.get_cmap("cubehelix_r")

# ----------------------
# Heatmap 2: Events 0-10, concatenated cmap
# ----------------------
vmin = 0
vmax = 8
transition_value = 1.0  # valore dove parte la seconda parte della colormap

n_low = 128
n_high = 128

# prima parte: 0 -> transition_value con cubehelix_r
low_colors = plt.get_cmap("cubehelix_r")(np.linspace(0,1,n_low))

# seconda parte: transition_value -> vmax con gnuplot
high_colors = plt.get_cmap("gnuplot")(np.linspace(0,1,n_high))

# concateno le due parti
colors_combined = np.vstack([low_colors, high_colors])

# posizioni normalizzate 0-1
transition_norm = (transition_value - vmin) / (vmax - vmin)
positions_low = np.linspace(0, transition_norm, n_low)
positions_high = np.linspace(transition_norm, 1.0, n_high)
positions = np.concatenate([positions_low, positions_high])

# colormap finale
cmap_events = LinearSegmentedColormap.from_list("cmap_events", list(zip(positions, colors_combined)))
cmap_events =  plt.get_cmap("gist_stern_r")
# Prendo solo il 95% della colormap (elimino la parte finale nera)
cmap_events = LinearSegmentedColormap.from_list(
    "gist_stern_r_no_black",
    cmap_events(np.linspace(0, 0.95, 256)))

# Trimestral heatmap
sns.heatmap(data_mean, ax=axes[0], cmap=cmap_trimesters, annot=True, fmt=".1f",
            xticklabels=trimesters, yticklabels=modes, cbar_kws={'label':'Mean amplitude (cm)'},vmin=0, vmax=transition_value)
axes[0].set_title("Mean seiches amplitude per trimester (2020-2023)")
axes[0].set_xlabel("Trimesters")
axes[0].set_ylabel("Modes")
axes[0].tick_params(axis='x', rotation=0)  # horizontal mode labels

# Events heatmap
sns.heatmap(events, ax=axes[1], cmap=cmap_events, annot=True, fmt=".1f",
            xticklabels=events_names, yticklabels=modes, cbar_kws={'label':'Extreme event amplitude (cm)'},vmin=vmin, vmax=vmax)
axes[1].set_title("Mean seiches amplitude during extreme events")
axes[1].set_xlabel("Events")
axes[1].set_ylabel("Modes")
axes[1].tick_params(axis='x', rotation=45)  # solo rotazione
plt.setp(axes[1].get_xticklabels(), ha='right')  # allinea le etichette a destra

plt.tight_layout()
plt.savefig("heatmap_trimesters_events.png", dpi=300, bbox_inches='tight')
