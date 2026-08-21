#!/usr/bin/env python
# ---------------------------------------------------------------------------
# Sea level perturbation (tilting) usata come condizione iniziale per la
# procedura numerica: mappa di sossheig, un pannello per ogni istante
# disponibile nel/nei file di input.
#
# Uso:
#   python plot_sl_tilting.py [file1.nc file2.nc ...]
# senza argomenti usa la lista INFILES definita qui sotto.
# ---------------------------------------------------------------------------

import os
import sys
import glob

import numpy as np
import netCDF4 as nc
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm, BoundaryNorm, LinearSegmentedColormap

# ============================ INPUTS =======================================

# File di input (accetta wildcard). Sovrascrivibili da riga di comando.
INFILES = [
    "/work/cmcc/ag15419/exp/EAS9BT_med-modes_atmp_BF3/EXP00/20150104/"
    "model/medfs-eas9_1h_20150104_2D_grid_T.nc",
]

#INFILES = [
#     "/work/cmcc/ag15419/exp/EAS9BT_med-modes_atmp/EXP00/20150103/"
#     "model/medfs-eas9_1h_20150103_2D_grid_T.nc",
#]

#INFILES = ["/work/cmcc/ag15419/exp/EAS9BT_med-modes_atmp/EXP00/20150101/"
#       "model/medfs-eas9_1h_20150101_2D_grid_T.nc",
#]

#INFILES = [
#      "/data/cmcc/ag15419/exp_juno/fix_mfseas9_longrun_hmslp_2NT_AB/EXP00/20150104/"
#      "model/medfs-eas9_1h_20150104_2D_grid_T.nc",
#]

# Mesh mask NEMO (per lat/lon e maschera terra/mare)
MESH_MASK = "/work/cmcc/ag15419/VAA_paper/DATA0/mesh_mask.nc"

# Directory di output
OUTDIR = "/work/cmcc/ag15419/basin_modes_ssh_tilting/"

# Nome della variabile di livello del mare
SSH_VAR = "sossheig"

# Limiti della mappa (lon_min, lon_max, lat_min, lat_max).
# Metti None per usare l'estensione completa della griglia.
MAP_EXTENT = (-6.0, 36.5, 30.0, 46.5)

# Scala colori:
#   'global'  -> stessi limiti per tutti gli istanti (pannelli confrontabili)
#   'per_time'-> limiti ricalcolati a ogni istante
COLOR_SCALE = 'global'

# Limite massimo imposto a mano [m]. None -> calcolato dai dati.
VMAX_FIXED = None

# Percentile usato per definire il limite (evita che un singolo punto
# costiero anomalo schiacci tutta la palette). 100 = massimo assoluto.
VMAX_PERCENTILE = 99.5

# Numero di livelli discreti della palette (dispari -> livello centrato su 0)
N_LEVELS = 21

# Disegnare le isolinee
DRAW_CONTOURS = True
CONTOUR_STEP_RATIO = 0.2      # passo isolinee come frazione di vmax

# Mascherare la box atlantica (come nella procedura semi-analitica)
MASK_ATLANTIC = True

# Estetica
FIGSIZE = (11, 5.0)
DPI = 200
FS_TITLE = 15
FS_LABEL = 13
FS_TICK = 11
CMAP = 'RdBu_r'               # divergente, leggibile anche in daltonismo

# ===========================================================================


def mask_atlantic_box(mask):
    """Stessa maschera usata nella procedura semi-analitica."""
    ny, nx = mask.shape
    J, I = np.indices((ny, nx))
    atlantic = (I < 300) | ((I < 415) & (J > 250))
    out = mask.copy()
    out[atlantic] = False
    return out


def read_grid(mesh_mask_path):
    with nc.Dataset(mesh_mask_path) as ds:
        tmask = ds.variables['tmask'][0, 0, :, :].astype(bool)
        lat = np.array(ds.variables['nav_lat'][:])
        lon = np.array(ds.variables['nav_lon'][:])
    if MASK_ATLANTIC:
        tmask = mask_atlantic_box(tmask)
    return tmask, lat, lon


def read_ssh(files, ssh_var):
    """Legge sossheig da uno o piu' file, restituendo (ssh, labels)."""
    fields = []
    labels = []

    for f in files:
        with nc.Dataset(f) as ds:
            if ssh_var not in ds.variables:
                raise KeyError(f"Variabile '{ssh_var}' assente in {f}. "
                               f"Disponibili: {list(ds.variables)}")
            data = np.array(ds.variables[ssh_var][:])
            if data.ndim == 2:            # singolo istante senza asse tempo
                data = data[np.newaxis, ...]

            # Etichette temporali
            tlab = None
            for tname in ('time_counter', 'time'):
                if tname in ds.variables:
                    tvar = ds.variables[tname]
                    try:
                        dates = nc.num2date(tvar[:], tvar.units,
                                            only_use_cftime_datetimes=False)
                        tlab = [d.strftime('%Y-%m-%d %H:%M') for d in
                                np.atleast_1d(dates)]
                    except Exception:
                        tlab = None
                    break
            if tlab is None:
                base = os.path.basename(f)
                tlab = [f"{base} step {i}" for i in range(data.shape[0])]

        fields.append(data)
        labels.extend(tlab)

    ssh = np.concatenate(fields, axis=0)
    return ssh, labels


def compute_vmax(ssh, tmask):
    """Limite simmetrico della palette, robusto agli outlier."""
    vals = np.abs(ssh[:, tmask])
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return 1.0
    if VMAX_PERCENTILE >= 100:
        v = np.max(vals)
    else:
        v = np.percentile(vals, VMAX_PERCENTILE)
    # arrotonda a 2 cifre significative, verso l'alto
    if v <= 0:
        return 1.0
    expo = np.floor(np.log10(v))
    return float(np.ceil(v / 10 ** (expo - 1)) * 10 ** (expo - 1))


def make_cmap(name, bad_color='0.75'):
    """Copia di una colormap con colore per i valori mascherati.

    Costruita ricampionando la palette invece di usare cmap.copy(), che
    esiste solo da Matplotlib 3.4 in poi.
    """
    base = plt.get_cmap(name)
    cmap = LinearSegmentedColormap.from_list(
        base.name + '_bad', base(np.linspace(0, 1, 256)))
    cmap.set_bad(bad_color)
    return cmap


def plot_one(field, tmask, lat, lon, vmax, title, outpath):

    masked = np.ma.masked_where(~tmask, field)

    levels = np.linspace(-vmax, vmax, N_LEVELS)
    norm = BoundaryNorm(levels, ncolors=256, clip=True)

    cmap = make_cmap(CMAP, '0.75')          # terra in grigio

    fig, ax = plt.subplots(figsize=FIGSIZE)

    try:
        im = ax.pcolormesh(lon, lat, masked, cmap=cmap, norm=norm,
                           shading='auto')
    except (TypeError, ValueError):
        # shading='auto' richiede Matplotlib >= 3.3
        im = ax.pcolormesh(lon, lat, masked, cmap=cmap, norm=norm)

    if DRAW_CONTOURS:
        step = CONTOUR_STEP_RATIO * vmax
        clev = np.arange(-vmax, vmax + step / 2, step)
        clev = clev[np.abs(clev) > 1e-12]     # salta lo zero
        cs = ax.contour(lon, lat, masked, levels=clev, colors='k',
                        linewidths=0.4, linestyles='solid')
        ax.clabel(cs, inline=True, fontsize=7, fmt='%.02f')

    # Linea di costa: contorno della maschera
    ax.contour(lon, lat, tmask.astype(float), levels=[0.5],
               colors='black', linewidths=0.7)

    cbar = fig.colorbar(im, ax=ax, pad=0.02, extend='both')
    cbar.set_label('SSH (m)', fontsize=FS_LABEL)
    cbar.ax.tick_params(labelsize=FS_TICK)

    if MAP_EXTENT is not None:
        ax.set_xlim(MAP_EXTENT[0], MAP_EXTENT[1])
        ax.set_ylim(MAP_EXTENT[2], MAP_EXTENT[3])

    ax.set_xlabel('Longitude', fontsize=FS_LABEL)
    ax.set_ylabel('Latitude', fontsize=FS_LABEL)
    ax.tick_params(labelsize=FS_TICK)
    ax.set_title(title, fontsize=FS_TITLE)
    ax.set_facecolor('0.75')

    fig.tight_layout()
    fig.savefig(outpath, dpi=DPI, bbox_inches='tight')
    plt.close(fig)


def main():
    files = sys.argv[1:] if len(sys.argv) > 1 else INFILES

    expanded = []
    for f in files:
        hits = sorted(glob.glob(f))
        if not hits:
            print(f"ATTENZIONE: nessun file per {f}")
        expanded.extend(hits)

    if not expanded:
        raise SystemExit("Nessun file di input trovato.")

    print(f"File letti: {len(expanded)}")
    for f in expanded:
        print("  ", f)

    os.makedirs(OUTDIR, exist_ok=True)

    tmask, lat, lon = read_grid(MESH_MASK)
    ssh, labels = read_ssh(expanded, SSH_VAR)
    print(f"Istanti disponibili: {ssh.shape[0]}  (griglia {ssh.shape[1:]})")

    if ssh.shape[1:] != tmask.shape:
        raise ValueError(f"Shape incoerenti: ssh {ssh.shape[1:]} "
                         f"vs mesh_mask {tmask.shape}")

    if VMAX_FIXED is not None:
        vmax_global = VMAX_FIXED
    else:
        vmax_global = compute_vmax(ssh, tmask)
    print(f"Limite palette (global): +/- {vmax_global:g} m")

    for k in range(ssh.shape[0]):
        if COLOR_SCALE == 'global' or VMAX_FIXED is not None:
            vmax = vmax_global
        else:
            vmax = compute_vmax(ssh[k:k + 1], tmask)

        title = f"Sea Level Perturbation" # - {labels[k]}"
        # Etichetta temporale ripulita per il nome del file
        stamp = (labels[k].replace('-', '').replace(':', '')
                          .replace(' ', '_').replace('.', ''))
        stamp = "".join(c for c in stamp if c.isalnum() or c == '_')
        outpath = os.path.join(OUTDIR, f"sl_tilting_{k:03d}_{stamp}.png")

        plot_one(ssh[k], tmask, lat, lon, vmax, title, outpath)
        print(f"  salvato {outpath}   (max |SSH| = "
              f"{np.nanmax(np.abs(ssh[k][tmask])):.4f} m)")

    print("Fatto.")


if __name__ == '__main__':
    main()
