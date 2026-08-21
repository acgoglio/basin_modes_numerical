import numpy as np
import re
import xarray as xr
import pandas as pd
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap
import netCDF4 as ncdf
from matplotlib.colors import LogNorm
import shutil
from area_ini import *

mpl.use("Agg")  # For non-interactive backend

# --- Paths ---
indir = work_dir
infile = os.path.join(indir, "basin_modes_amp_med.nc")
csvfile = os.path.join(indir, "periods_grouped_amp.csv")
outfile = os.path.join(indir, "mode_groups_amp.nc")
mesh_mask_file = mesh_mask
bathy_file = bathy_meter
output_plot_dir = os.path.join(indir, "mode_plots_amp")

# If the output directory exists, clear its contents
if os.path.exists(output_plot_dir):
    shutil.rmtree(output_plot_dir)
# Recreate the output  directory
os.makedirs(output_plot_dir, exist_ok=True)

# ----------------------------------------------------------------------
# FIX A: l'indice del modo veniva letto con int(var[1]), cioe' UN SOLO
# carattere. Per 'm10_Amp' restituiva '1', quindi m1/m10/m11... collidevano
# sulla stessa chiave del dizionario e vinceva l'ultima letta: le mappe dei
# modi con indice >= 10 erano attribuite al modo sbagliato.
# ----------------------------------------------------------------------
def mode_index(var):
    m = re.match(r'm(\d+)_', var)
    if m is None:
        raise ValueError(f"Unexpected variable name: {var}")
    return int(m.group(1))

# ----------------------------------------------------------------------
# Compatibilita' colormap: mpl.cm.get_cmap e' deprecato da matplotlib 3.7
# ma mpl.colormaps esiste solo da 3.5. Questa funzione usa l'API nuova se
# disponibile e ricade su quella vecchia altrimenti.
# ----------------------------------------------------------------------
def get_cmap(name):
    try:
        return mpl.colormaps[name]          # matplotlib >= 3.5
    except AttributeError:
        return mpl.cm.get_cmap(name)        # matplotlib < 3.5

# Truncate the colormap to exclude the lightest part (e.g. bottom 20%)
def truncate_colormap(cmap, minval=0.2, maxval=1.0, n=256):
    new_cmap = LinearSegmentedColormap.from_list(
        f"trunc({cmap.name},{minval:.2f},{maxval:.2f})",
        cmap(np.linspace(minval, maxval, n))
    )
    return new_cmap

# ----------------------------------------------------------------------
# Soglia di copertura: i gruppi rilevati su meno di MIN_COV punti mare del
# Mediterraneo non vengono plottati. Stesso valore usato in
# mode_period_tab_amp.py per la Fig. 3.5, cosi' le mappe prodotte sono
# esattamente i modi che quella figura tiene.
# Il netCDF mode_groups_amp.nc resta COMPLETO: il filtro agisce solo sui
# plot, quindi la scelta e' reversibile senza rigenerare i campi.
# ----------------------------------------------------------------------
MIN_COV = 7.0            # [%] copertura minima per produrre le mappe

# ----------------------------------------------------------------------
# Precisione dei periodi nei titoli, legata alla tolleranza del gruppo:
#     prec = -floor(log10(tol))
# cioe' tante cifre decimali quante ne servono per la propria tolleranza.
# Stessa logica delle etichette della Fig. 3.5 (mode_period_tab_amp.py),
# cosi' titoli e figura riportano i periodi con la stessa precisione, e
# tutto si adatta da solo se cambia segment_len_days.
# Serviva anche perche' con una cifra fissa i gruppi a 3.58 e 3.63 h
# diventavano entrambi "3.6 h": due mappe diverse con lo stesso titolo.
# ----------------------------------------------------------------------
def _prec_from_tol(tol):
    """Cifre decimali sufficienti a rappresentare la tolleranza."""
    try:
        tol = float(tol)
    except (ValueError, TypeError):
        return 2
    if not np.isfinite(tol) or tol <= 0:
        return 2
    return max(0, -int(np.floor(np.log10(tol))))

def fmt_period(gp, tol):
    return f"{gp:.{_prec_from_tol(tol)}f}"

def fmt_tol(tol):
    """Tolleranza con la sua stessa precisione (evita '0.03000000000001')."""
    return f"{tol:.{_prec_from_tol(tol)}f}"

# --- Load datasets ---
ds = xr.open_dataset(infile)
grouped_df = pd.read_csv(csvfile)
try:
    group_centers = grouped_df["Grouped_Period"].values
    tolerance = grouped_df["Tolerance"].values
except:
    group_centers = grouped_df["Grouped_Period"].values
    tolerance = np.full(len(group_centers), fixed_uncertainty)

# --- Copertura di ciascun gruppo, dal CSV prodotto da mode_period_tab_amp ---
coverage_pct = None
for _col in ("Coverage_pct", "Percentage", "%"):
    if _col in grouped_df.columns:
        coverage_pct = grouped_df[_col].values.astype(float)
        print(f"Copertura letta dalla colonna '{_col}' di {os.path.basename(csvfile)}")
        break

if coverage_pct is None:
    print("ATTENZIONE: nessuna colonna di copertura nel CSV, "
          "nessun gruppo viene escluso.")
    keep_group = np.ones(len(group_centers), dtype=bool)
else:
    keep_group = coverage_pct >= MIN_COV
    _skip = np.where(~keep_group)[0]
    print(f"Gruppi con copertura >= {MIN_COV}%: "
          f"{int(keep_group.sum())} su {len(group_centers)}")
    if len(_skip) > 0:
        print("   esclusi dai plot: " + ", ".join(
            f"{group_centers[i]:.2f} h ({coverage_pct[i]:.1f}%)" for i in _skip))

# --- Load grid and mask ---
nc_nemogrid = ncdf.Dataset(mesh_mask_file)
nav_lat = nc_nemogrid.variables["nav_lat"][:]
nav_lon = nc_nemogrid.variables["nav_lon"][:]
tmask = nc_nemogrid.variables["tmask"][0, 0, :, :]  # shape (y, x)
field_shape = nav_lat.shape

# Load bathymetry
ds_bathy = xr.open_dataset(bathy_file, decode_times=False)
bathy = ds_bathy["Bathymetry"][0, :, :]

# --- Initialize output fields ---
fields = {f"mode_{gp:.2f}h": np.full(field_shape, np.nan) for gp in group_centers}
amp_fields = {f"amp_{gp:.2f}h": np.full(field_shape, np.nan) for gp in group_centers}

# --- Extract variables ---
modes_vars = [var for var in ds.data_vars if var.startswith("m") and "_T" in var]
amp_vars = {mode_index(var): var for var in ds.data_vars if var.startswith("m") and "_Amp" in var}

# --- Assign values ---
for var in modes_vars:
    mode_num = mode_index(var)
    period_data = ds[var].values

    if np.issubdtype(period_data.dtype, np.timedelta64):
        period_data = period_data.astype("timedelta64[s]").astype(float) / 3600
    else:
        period_data = period_data.astype(float)
    period_data = np.round(period_data, 2)

    if mode_num not in amp_vars:
        continue
    amp_data = ds[amp_vars[mode_num]].values

    written_mask = np.zeros(field_shape, dtype=bool)

    for idx_gp,gp in enumerate(group_centers):
        match = np.abs(period_data - gp) <= tolerance[idx_gp]
        match = np.logical_and(match, ~written_mask)

        if np.any(match):
            print(f"Mode {mode_num} matches group {gp:.2f}h in {np.sum(match)} points")

        fields[f"mode_{gp:.2f}h"][match] = mode_num
        amp_fields[f"amp_{gp:.2f}h"][match] = np.fmax(
            amp_fields[f"amp_{gp:.2f}h"][match], amp_data[match]
        )

        written_mask[match] = True

# --- Save all fields to NetCDF ---
combined_fields = {
    **{name: (("y", "x"), data) for name, data in fields.items()},
    **{name: (("y", "x"), data) for name, data in amp_fields.items()}
}
output_ds = xr.Dataset(
    combined_fields,
    coords={"nav_lat": (("y", "x"), nav_lat), "nav_lon": (("y", "x"), nav_lon)}
)
output_ds.to_netcdf(outfile)
print(f"Saved grouped mode and amplitude maps to: {outfile}")

all_amp_vals=[]

# --- Plot each grouped mode field ---
for idx_gp, gp in enumerate(group_centers):
    if not keep_group[idx_gp]:
        continue
    mode_field  = fields[f"mode_{gp:.2f}h"]
    amp_field   = amp_fields[f"amp_{gp:.2f}h"]
    uncertainty = tolerance[idx_gp] 

    # Modes
    #plt.figure(figsize=(10, 4))
    #cmap = get_cmap("Reds_r")
    #cmap.set_bad("white")
    #masked_field = np.ma.masked_invalid(mode_field)
    #plt.contourf(nav_lon, nav_lat, masked_field, levels=np.arange(0, 9), cmap=cmap, extend="neither")
    #plt.colorbar(label="Mode Number")
    #plt.contourf(nav_lon, nav_lat, tmask, levels=[-1000, 0.05], colors="gray")
    #plt.contour(nav_lon, nav_lat, tmask, levels=[0.05], colors="black", linewidths=0.8)
    #plt.title(f"Modes with Period: {fmt_period(gp, uncertainty)} h ± {fmt_tol(uncertainty)} h")
    #plt.xlabel("Longitude")
    #plt.ylabel("Latitude")
    #plt.xlim(-6, 36.3)
    #plt.ylim(30, 46)
    #plt.savefig(os.path.join(output_plot_dir, f"mode_amp_{idx_gp}_{gp:.2f}h.png"), dpi=300, bbox_inches="tight")
    #plt.close()

    # Presence 
    presence_mask = (~np.isnan(mode_field)).astype(int)
    #plt.figure(figsize=(10, 4))
    cmap_presence = mpl.colors.ListedColormap(["white", "tab:blue"])
    bounds = [0, 0.5, 1.5]
    norm = mpl.colors.BoundaryNorm(bounds, cmap_presence.N)
    #plt.contourf(nav_lon, nav_lat, presence_mask, levels=bounds, cmap=cmap_presence, norm=norm)
    #plt.contourf(nav_lon, nav_lat, tmask, levels=[-1000, 0.05], colors="gray")
    #plt.contour(nav_lon, nav_lat, tmask, levels=[0.05], colors="black", linewidths=0.8)
    #plt.xlabel("Longitude")
    #plt.ylabel("Latitude")
    #plt.xlim(-6, 36.3)
    #plt.ylim(30, 46)

    # Grid points count
    # Define spatial boundaries for the Mediterranean domain
    lon_min, lon_max = -6, 36.3
    lat_min, lat_max = 30, 46
    # Select only points within the Mediterranean bounding box
    spatial_mask = (nav_lon >= lon_min) & (nav_lon <= lon_max) & \
                   (nav_lat >= lat_min) & (nav_lat <= lat_max)
    # Exclude the Atlantic Ocean area (west of 0°E and north of 42°N)
    atlantic_exclusion = (nav_lon < 0) & (nav_lat > 42)
    # Select valid sea points in the Mediterranean, excluding the Atlantic
    sea_mask = (tmask == 1) & spatial_mask & (~atlantic_exclusion)
    # Mask where the mode field is valid (not NaN) and corresponds to sea points
    presence_mask_valid = (~np.isnan(mode_field)) & sea_mask
    # Count total sea points (excluding Atlantic)
    total_sea_points = np.sum(sea_mask)
    # Count number of valid (non-NaN) mode points in the sea
    presence_count = np.sum(presence_mask_valid)
    # %
    if total_sea_points > 0:
       presence_percentage = 100 * presence_count / total_sea_points
    else:
       presence_percentage = np.nan
    # Output
    print(f"Presence points: {presence_count}")
    print(f"Total sea points: {total_sea_points}")
    print(f"Presence percentage: {presence_percentage:.1f}%") 
   
    #plt.title(f"Presence Map for Period: {fmt_period(gp, uncertainty)} h ± {fmt_tol(uncertainty)} h ({presence_percentage:.1f}%)")
    #plt.savefig(os.path.join(output_plot_dir, f"mode_flag_amp_{idx_gp}_{gp:.2f}h.png"), dpi=300, bbox_inches="tight") #{idx_gp}_{gp:.2f}h.png"), dpi=300, bbox_inches="tight")
    #plt.close()


    # Rel amplitude Med basin
    plt.figure(figsize=(10, 4))
    all_amp_vals.append(amp_field)
    # Normalize to the max or to the 99th perc (or to 100)
    #amp_percent = amp_field / np.nanmax(amp_field) * 100
    # Normalizzazione sul 99o percentile invece che sul massimo: rende le
    # mappe leggibili anche quando pochi punti hanno ampiezze molto maggiori
    # del resto del bacino. I valori sopra il 99o percentile sono poi
    # riportati a 100 dal clip piu' sotto.
    amp_percent = amp_field / np.nanpercentile(amp_field, 99) * 100
    masked_amp_pct = np.ma.masked_invalid(amp_percent)
    cmap_amp_pct = get_cmap("gist_stern_r")
    cmap_amp_pct = truncate_colormap(cmap_amp_pct, 0.05, 0.95)
    cmap_amp_pct.set_bad("white")

    # Handle values > the 99th perc
    clipped_masked_amp_pct = np.clip(masked_amp_pct, 0, 100)
    # Plot
    im = plt.contourf(nav_lon, nav_lat, clipped_masked_amp_pct, levels=np.linspace(0, 100, 41), cmap=cmap_amp_pct)
    plt.colorbar(im, label="Amplitude (% of 99th percentile)")
    #plt.contour(nav_lon, nav_lat, masked_amp_pct, levels=[0.00], colors="magenta", linewidths=1.0)
    plt.contourf(nav_lon, nav_lat, tmask, levels=[-1000, 0.05], colors="gray")
    plt.contour(nav_lon, nav_lat, tmask, levels=[0.05], colors="black", linewidths=0.8)

    # Add % dashed contours
    contour_levels_perc = np.arange(0, 101, 10)
    valid_levels = [lev for lev in contour_levels_perc
        if np.nanmin(masked_amp_pct) <= lev <= np.nanmax(masked_amp_pct)]
    if valid_levels:
        CS_perc = plt.contour(nav_lon, nav_lat, masked_amp_pct,
            levels=valid_levels, colors='k', linestyles='dashed', linewidths=0.5)
        # FIX D: clabel stava fuori dalla guardia `if valid_levels`, quindi
        # con valid_levels vuoto usava un CS_perc stantio (o dava NameError).
        plt.clabel(CS_perc, inline=True, fontsize=7, fmt='%d%%')

    #contour_levels = [100, 200, 300]
    #plt.contour(nav_lon, nav_lat, bathy, levels=contour_levels, colors="black", linewidths=0.5, linestyles="dashed")
    #contour_levels = [400, 500]
    #plt.contour(nav_lon, nav_lat, bathy, levels=contour_levels, colors="black", linewidths=0.5, linestyles="dotted")

    plt.title(f"Amplitude for Mode with Period: {fmt_period(gp, uncertainty)} h ± {fmt_tol(uncertainty)} h")
    plt.xlabel("Longitude")
    plt.ylabel("Latitude")
    plt.xlim(-6, 36.3)
    plt.ylim(30, 46)
    plt.savefig(os.path.join(output_plot_dir, f"mode_amprelval_{idx_gp}_{gp:.2f}h.png"), dpi=300, bbox_inches="tight") #{idx_gp}_{gp:.2f}h.png"), dpi=300, bbox_inches="tight")
    plt.close()

   # Rel amplitude Adriatic Sea
    plt.figure(figsize=(5, 4))
    im = plt.contourf(nav_lon, nav_lat, clipped_masked_amp_pct, levels=np.linspace(0, 100, 41), cmap=cmap_amp_pct)
    plt.colorbar(im, label="Amplitude (% of 99th percentile)")
    plt.contourf(nav_lon, nav_lat, tmask, levels=[-1000, 0.05], colors="gray")
    plt.contour(nav_lon, nav_lat, tmask, levels=[0.05], colors="black", linewidths=0.8)
    if valid_levels:
        CS_perc = plt.contour(nav_lon, nav_lat, masked_amp_pct,
            levels=valid_levels, colors='k', linestyles='dashed', linewidths=0.5)
        plt.clabel(CS_perc, inline=True, fontsize=7, fmt='%d%%', inline_spacing=0.2, manual=False)
    plt.title(f"Amplitude for Mode with Period: {fmt_period(gp, uncertainty)} h ± {fmt_tol(uncertainty)} h")
    plt.xlabel("Longitude")
    plt.ylabel("Latitude")
    plt.xlim(12.0, 20.0)
    plt.ylim(39, 46)
    plt.savefig(os.path.join(output_plot_dir, f"mode_amprelval_{idx_gp}_{gp:.2f}h_AdrSea.png"), dpi=300, bbox_inches="tight") #{idx_gp}_{gp:.2f}h.png"), dpi=300, bbox_inches="tight")
    plt.close()


# ----------------------------------------------------------------------
# Massimo assoluto usato per normalizzare le mappe "absolute amplitude".
# I valori oltre AMP_OUTLIER_M sono esclusi: senza questo filtro un singolo
# punto anomalo (tipicamente in bacini chiusi o punti isolati) schiaccerebbe
# la scala colori di TUTTE le mappe.
# La soglia era hardcoded a 2: qui e' un parametro, e si stampa quanti punti
# vengono esclusi. Se il numero non e' trascurabile non e' un dettaglio
# grafico ma il sintomo di un problema, e va guardato.
# NB: con la normalizzazione corretta (N_real invece di N_fft) le ampiezze
# crescono di circa il 10%, quindi valori appena sotto la soglia possono
# ora superarla.
# ----------------------------------------------------------------------
AMP_OUTLIER_M = 2.0        # [m] ampiezze oltre questa soglia non fissano la scala

all_amp_vals = np.array(all_amp_vals)
mask_all_amp_vals = all_amp_vals < AMP_OUTLIER_M

n_valid   = int(np.sum(np.isfinite(all_amp_vals)))
n_outlier = int(np.sum(np.isfinite(all_amp_vals) & ~mask_all_amp_vals))

if not np.any(mask_all_amp_vals):
    raise RuntimeError(
        f"Nessuna ampiezza sotto {AMP_OUTLIER_M} m: soglia da rivedere.")

all_amp_vals_max = np.nanmax(all_amp_vals[mask_all_amp_vals])
print(f'Abs Amp Max: {all_amp_vals_max:.4f} m  '
      f'(su {n_valid} valori validi)')
if n_outlier > 0:
    _out = all_amp_vals[np.isfinite(all_amp_vals) & ~mask_all_amp_vals]
    print(f'   esclusi {n_outlier} valori >= {AMP_OUTLIER_M} m '
          f'({100.0*n_outlier/max(n_valid,1):.3f}% del totale), '
          f'massimo escluso {np.max(_out):.2f} m')
else:
    print(f'   nessun valore >= {AMP_OUTLIER_M} m')

# --- Plot each grouped mode field for abs plots---
for idx_gp, gp in enumerate(group_centers):
    if not keep_group[idx_gp]:
        continue
    mode_field = fields[f"mode_{gp:.2f}h"]
    amp_field = amp_fields[f"amp_{gp:.2f}h"]
    # FIX B: `uncertainty` non veniva ridefinita in questo loop, quindi tutti
    # i titoli riportavano la tolleranza dell'ULTIMO gruppo del loop
    # precedente invece della propria.
    uncertainty = tolerance[idx_gp]

    # Abs amplitude
    plt.figure(figsize=(10, 4))
    amp_percent = amp_field / all_amp_vals_max * 100
    masked_amp_pct = np.ma.masked_invalid(amp_percent)
    cmap_amp_pct = get_cmap("gist_stern_r")
    cmap_amp_pct = truncate_colormap(cmap_amp_pct, 0.05, 0.95)
    cmap_amp_pct.set_bad("white")

    im = plt.contourf(nav_lon, nav_lat, masked_amp_pct, levels=np.linspace(0, 100, 41), cmap=cmap_amp_pct)
    plt.colorbar(im, label="Amplitude (% of abs max)")
    plt.contour(nav_lon, nav_lat, masked_amp_pct, levels=[0.00], colors="magenta", linewidths=1.0)
    plt.contourf(nav_lon, nav_lat, tmask, levels=[-1000, 0.05], colors="gray")
    plt.contour(nav_lon, nav_lat, tmask, levels=[0.05], colors="black", linewidths=0.8)

    contour_levels = [100, 200, 300]
    #plt.contour(nav_lon, nav_lat, bathy, levels=contour_levels, colors="black", linewidths=0.5, linestyles="dashed")
    contour_levels = [400, 500]
    #plt.contour(nav_lon, nav_lat, bathy, levels=contour_levels, colors="black", linewidths=0.5, linestyles="dotted")

    plt.title(f"Absolute Amplitude for Mode with Period: {fmt_period(gp, uncertainty)} h ± {fmt_tol(uncertainty)} h")
    plt.xlabel("Longitude")
    plt.ylabel("Latitude")
    plt.xlim(-6, 36.3)
    plt.ylim(30, 46)
    plt.savefig(os.path.join(output_plot_dir, f"mode_absampval_{idx_gp}_{gp:.2f}h.png"), dpi=300, bbox_inches="tight") #{idx_gp}_{gp:.2f}h.png"), dpi=300, bbox_inches="tight")
    plt.close()

# --- Plot each grouped mode field for plots of real amplitude---
for idx_gp, gp in enumerate(group_centers):
    if not keep_group[idx_gp]:
        continue
    mode_field = fields[f"mode_{gp:.2f}h"]
    amp_field = amp_fields[f"amp_{gp:.2f}h"]
    # FIX B: `uncertainty` non veniva ridefinita in questo loop, quindi tutti
    # i titoli riportavano la tolleranza dell'ULTIMO gruppo del loop
    # precedente invece della propria.
    uncertainty = tolerance[idx_gp]

    # Abs amplitude
    plt.figure(figsize=(10, 4))
    amp_percent = amp_field * 100
    masked_amp_pct = np.ma.masked_invalid(amp_percent)
    cmap_amp_pct = get_cmap("cubehelix_r")
    cmap_amp_pct = truncate_colormap(cmap_amp_pct, 0.05, 0.95)
    cmap_amp_pct.set_bad("white")

    im = plt.contourf(nav_lon, nav_lat, masked_amp_pct, levels=np.linspace(0, np.nanmax(amp_field*100), 41), cmap=cmap_amp_pct)
    plt.colorbar(im, label="Amplitude (cm)")
    plt.contour(nav_lon, nav_lat, masked_amp_pct, levels=[0.00], colors="magenta", linewidths=1.0)
    plt.contourf(nav_lon, nav_lat, tmask, levels=[-1000, 0.05], colors="gray")
    plt.contour(nav_lon, nav_lat, tmask, levels=[0.05], colors="black", linewidths=0.8)

    contour_levels = [100, 200, 300]
    #plt.contour(nav_lon, nav_lat, bathy, levels=contour_levels, colors="black", linewidths=0.5, linestyles="dashed")
    contour_levels = [400, 500]
    #plt.contour(nav_lon, nav_lat, bathy, levels=contour_levels, colors="black", linewidths=0.5, linestyles="dotted")

    plt.title(f"Amplitude for Mode with Period: {fmt_period(gp, uncertainty)} h ± {fmt_tol(uncertainty)} h")
    plt.xlabel("Longitude")
    plt.ylabel("Latitude")
    plt.xlim(-6, 36.3)
    plt.ylim(30, 46)
    plt.savefig(os.path.join(output_plot_dir, f"mode_realampval_{idx_gp}_{gp:.2f}h.png"), dpi=300, bbox_inches="tight") #{idx_gp}_{gp:.2f}h.png"), dpi=300, bbox_inches="tight")
    plt.close()
