#!/usr/bin/env python3
# ======================================================================
# Verifica di sea_gp_num (valore hardcoded in mode_period_tab_*.py) e
# ripartizione per regione dei punti mare nel dominio dei box.
#
# Perche' serve:
#  - i box di run_area_BM.sh coprono lon idx 300..1305, cioe' da ~5.6W:
#    escludono quasi tutta la box atlantica MA includono la porzione
#    sud-orientale del Golfo di Biscaglia, che sta a nord di 43N e a ovest
#    di ~1W ed e' dentro il dominio (il modello si ferma a 45.98N);
#  - se sea_gp_num conta il solo Mediterraneo mentre i box cercano i modi
#    anche in Biscaglia, numeratore e denominatore stanno su domini
#    diversi e una copertura puo' superare il 100%.
#
# Uso:  python check_sea_gp_num.py
# ======================================================================

import numpy as np
import xarray as xr
from area_ini import mesh_mask

SEA_GP_NUM_USED = 145023        # valore attualmente in mode_period_tab_*.py

# indici dei box, come in run_area_BM.sh
LON_MIN, LON_MAX = 300, 1306
LAT_MIN, LAT_MAX = 0, 379

# confini indicativi per la classificazione geografica
LAT_BISCAY  = 43.0     # a nord di questa e a ovest di LON_BISCAY -> Biscaglia
LON_BISCAY  = -1.0
LON_GIBR    = -5.35    # longitudine dello Stretto di Gibilterra
LAT_GIBR    = 40.0     # a sud di questa e a ovest di LON_GIBR -> Atlantico SW

ds = xr.open_dataset(mesh_mask)

name = next((v for v in ("tmask", "tmaskutil", "top_level") if v in ds), None)
if name is None:
    raise KeyError(f"Nessuna maschera trovata. Variabili: {list(ds.data_vars)}")
print(f"Maschera: {name}   dims={ds[name].dims}   shape={ds[name].shape}")

m = ds[name].values
while m.ndim > 2:               # scarta time/depth, tiene il livello di superficie
    m = m[0]
m = (m > 0)
ny, nx = m.shape
print(f"Dominio del modello: {nx} x {ny} punti\n")

# coordinate
lon = ds["nav_lon"].values if "nav_lon" in ds else None
lat = ds["nav_lat"].values if "nav_lat" in ds else None
if lon is None or lat is None:
    print("ATTENZIONE: nav_lon/nav_lat assenti, ripartizione geografica saltata.\n")
else:
    while lon.ndim > 2: lon = lon[0]
    while lat.ndim > 2: lat = lat[0]

# ---------------------------------------------------------------- totali
tot_dom   = int(m.sum())
sub       = m[LAT_MIN:LAT_MAX, LON_MIN:LON_MAX]
tot_boxes = int(sub.sum())
tot_atl   = int(m[LAT_MIN:LAT_MAX, 0:LON_MIN].sum())

print(f"  punti mare, dominio intero              : {tot_dom:>8,}")
print(f"  punti mare, box atlantica (lon idx<{LON_MIN}) : {tot_atl:>8,}")
print(f"  punti mare, dominio dei box             : {tot_boxes:>8,}")
print(f"  valore usato in mode_period_tab_*.py    : {SEA_GP_NUM_USED:>8,}")

# ------------------------------------------------- ripartizione regionale
if lon is not None:
    slon = lon[LAT_MIN:LAT_MAX, LON_MIN:LON_MAX]
    slat = lat[LAT_MIN:LAT_MAX, LON_MIN:LON_MAX]

    biscay = sub & (slat > LAT_BISCAY) & (slon < LON_BISCAY)
    atl_sw = sub & (slat < LAT_GIBR)   & (slon < LON_GIBR)
    med    = sub & ~biscay & ~atl_sw

    n_bis, n_asw, n_med = int(biscay.sum()), int(atl_sw.sum()), int(med.sum())
    print("\n  Ripartizione del dominio dei box:")
    print(f"    Mediterraneo                          : {n_med:>8,}  ({100*n_med/tot_boxes:5.2f}%)")
    print(f"    Golfo di Biscaglia (lat>{LAT_BISCAY:.0f}N, lon<{LON_BISCAY:.0f}) : {n_bis:>8,}  ({100*n_bis/tot_boxes:5.2f}%)")
    print(f"    Atlantico a ovest di Gibilterra       : {n_asw:>8,}  ({100*n_asw/tot_boxes:5.2f}%)")

    if n_bis > 0:
        print("\n  --> I box includono punti del Golfo di Biscaglia.")
        print("      L'etichetta \"% of Med. sea grid points\" non e' esatta:")
        print("      o si usa il totale dei box come denominatore, oppure si")
        print("      escludono quei punti anche dal numeratore.")

# ------------------------------------------------------------ conclusione
print()
if tot_boxes == SEA_GP_NUM_USED:
    print("  OK: sea_gp_num coincide con il totale del dominio dei box.")
elif lon is not None and n_med == SEA_GP_NUM_USED:
    print("  sea_gp_num corrisponde al SOLO Mediterraneo, ma i box coprono")
    print(f"  {tot_boxes:,} punti: le coperture possono superare il 100%.")
else:
    print(f"  Differenza rispetto al dominio dei box: {tot_boxes-SEA_GP_NUM_USED:+,}")
    if abs(tot_dom - SEA_GP_NUM_USED) < abs(tot_boxes - SEA_GP_NUM_USED):
        print("  Il valore usato somiglia al conteggio sul DOMINIO INTERO:")
        print("  includerebbe la box atlantica -> percentuali sottostimate.")
    print(f"\n  Valore coerente col numeratore: sea_gp_num = {tot_boxes}")
