"""
Driver principale per l'analisi dei modi sulla run FORZATA (§3.3.2).

Uso (stesso pattern CLI di run_basin_modes_amp_idx.py, per compatibilita'
con l'infrastruttura di lancio a box sul cluster):

    python run_forced_modes.py <min_lon> <max_lon> <min_lat> <max_lat> <box_idx>

Cicla INTERNAMENTE su tutte le 27 finestre (consolidamento deciso in
chat: da 27 lanci separati a pochi lanci, ciascuno che itera sulle
finestre) - per il box di griglia assegnato a questo job.

Per ciascuna finestra produce:
  - una mappa 2D (lat,lon) di ampiezza assoluta per ciascuno dei 10
    modi di riferimento -> netCDF work_dir/forced_modes_<key>_<box_idx>.nc
  - un istogramma diagnostico di copertura dei picchi trovati (validazione
    visiva, stile Fig. 3.5) -> PNG (fatto SOLO se run senza splitting a
    box, o accumulando dati su tutti i box - vedi nota sotto)

ATTENZIONE - punto NON deciso, da verificare con te:
Se il job gira per singolo box (come l'infrastruttura esistente), la
mappa 2D per il 99th percentile serve calcolata sull'INTERO dominio Med,
non sul singolo box. Questo driver quindi scrive le mappe per box (come
la pipeline esistente, che poi le fondeva con merge_amp_idx.py), ma il
calcolo del 99th percentile finale e dell'istogramma diagnostico vanno
fatti DOPO il merge di tutti i box, non qui dentro. Ho lasciato un
placeholder per il merge (vedi merge_forced_modes.py, da scrivere) -
questo script si occupa solo della fase per-box.
"""

import os
import sys
import numpy as np
import xarray as xr
import netCDF4 as nc

import forced_windows_ini as ini
from f_forced_ampspt import (
    compute_point_spectrum,
    find_all_peaks,
    extract_reference_mode_amplitudes,
)
from build_file_list import build_file_list
from rt_stats_tools import get_med_mask

########################################################################
# CLI: box di griglia assegnato a questo job
########################################################################

min_lon = int(sys.argv[1])
max_lon = int(sys.argv[2])
min_lat = int(sys.argv[3])
max_lat = int(sys.argv[4])
box_idx = str(sys.argv[5])

os.makedirs(ini.work_dir, exist_ok=True)

########################################################################
# Griglia e maschera (una volta sola, riusata per tutte le finestre)
########################################################################

mesh_nemo = nc.Dataset(ini.mesh_mask, "r")
tmask = mesh_nemo.variables["tmask"][0, 0, :, :]
nav_lat_full = mesh_nemo.variables["nav_lat"][:]
nav_lon_full = mesh_nemo.variables["nav_lon"][:]
mesh_nemo.close()

# Maschera Mediterraneo (esclude Baia di Biscaglia / Mar Nero), come nel
# resto della pipeline - usata per delimitare i punti su cui e' valido
# calcolare/accumulare il 99th percentile e l'istogramma diagnostico.
med_mask_full = get_med_mask(nav_lon_full, nav_lat_full, tmask.astype(bool))

########################################################################
# Loop sulle 27 finestre
########################################################################

for kind, key, start_date, end_date in ini.WINDOWS:

    print(f"\n=== Finestra {kind}/{key}: {start_date}-{end_date} ===")

    spec_settings = ini.SPECTRUM_SETTINGS[kind]

    # ---- Lista esplicita dei file (niente glob) ----
    try:
        files, missing = build_file_list(
            start_date, end_date, ini.forced_run_daily_template
        )
    except FileNotFoundError as e:
        print(f"  SALTO finestra {key}: {e}")
        continue

    # ---- Apertura SSH solo per il box assegnato ----
    ds = xr.open_mfdataset(files, combine="by_coords", parallel=True)
    ssh_box = ds[ini.ssh_varname][:, min_lat:max_lat, min_lon:max_lon]

    ny_box = max_lat - min_lat
    nx_box = max_lon - min_lon

    # Output: una mappa per modo, per questo box
    amp_maps = np.full((ini.N_MODES, ny_box, nx_box), np.nan)
    fallback_maps = np.zeros((ini.N_MODES, ny_box, nx_box), dtype=bool)

    # Accumulo picchi per l'istogramma diagnostico (solo per questo box;
    # il merge finale dovra' concatenare tra i box - vedi nota in testa)
    all_peak_periods_box = []

    for j in range(ny_box):
        for i in range(nx_box):
            lat_idx = min_lat + j
            lon_idx = min_lon + i

            if tmask[lat_idx, lon_idx] != 1:
                continue  # punto di terra

            ssh_ts = ssh_box[:, j, i].values * 100.0  # m -> cm (sossheig e' in metri)

            freq_pos, periods, amplitudes = compute_point_spectrum(
                ssh_ts, ini.dt,
                flag_hanning=ini.flag_hanning,
                flag_nfft=ini.flag_nfft,
                N_fft=ini.N_fft,
                flag_segmented_spectrum=spec_settings["flag_segmented_spectrum"],
                segment_len_days=spec_settings["segment_len_days"],
                flag_filter=ini.flag_filter,
                th_filter=ini.th_filter,
            )

            if periods is None:
                # serie troppo corta/tutta NaN per questo punto/finestra
                continue

            peak_periods, peak_amps = find_all_peaks(periods, amplitudes)

            if med_mask_full[lat_idx, lon_idx]:
                all_peak_periods_box.append(peak_periods)

            mode_amps, used_fallback = extract_reference_mode_amplitudes(
                periods, amplitudes, peak_periods, peak_amps,
                ini.REFERENCE_MODES,
                flag_use_fallback=ini.flag_use_fallback,
            )

            amp_maps[:, j, i] = mode_amps
            fallback_maps[:, j, i] = used_fallback

    ds.close()

    # ---- Scrittura output per-box ----
    outfile = os.path.join(ini.work_dir, f"forced_modes_{kind}_{key}_{box_idx}.nc")
    with nc.Dataset(outfile, "w") as out:
        out.createDimension("mode", ini.N_MODES)
        out.createDimension("y", ny_box)
        out.createDimension("x", nx_box)

        v_amp = out.createVariable("amplitude", "f4", ("mode", "y", "x"), fill_value=np.nan)
        v_amp[:] = amp_maps
        v_amp.units = "cm"  # sossheig e' in metri, convertito a cm in lettura

        v_fb = out.createVariable("used_fallback", "i1", ("mode", "y", "x"))
        v_fb[:] = fallback_maps.astype(np.int8)

        out.setncattr("kind", kind)
        out.setncattr("key", key)
        out.setncattr("start_date", start_date)
        out.setncattr("end_date", end_date)
        out.setncattr("min_lon", min_lon)
        out.setncattr("min_lat", min_lat)

        for m, (label, T_ref, tol) in enumerate(ini.REFERENCE_MODES):
            out.setncattr(f"mode{m}_label", label)
            out.setncattr(f"mode{m}_T_ref_h", T_ref)
            out.setncattr(f"mode{m}_tolerance_h", tol)

    # Picchi grezzi del box, per il merge dell'istogramma diagnostico
    peaks_outfile = os.path.join(ini.work_dir, f"peaks_{kind}_{key}_{box_idx}.npz")
    if len(all_peak_periods_box) > 0:
        np.savez(peaks_outfile,
                 peak_periods=np.concatenate(all_peak_periods_box))
    else:
        np.savez(peaks_outfile, peak_periods=np.array([]))

    print(f"  Scritto: {outfile}")
    print(f"  Scritto: {peaks_outfile}")

print("\nFatto.")
