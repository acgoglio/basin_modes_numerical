"""
Driver principale per l'analisi dei modi sulla run FORZATA (§3.3.2).

Uso (stesso pattern CLI di run_basin_modes_amp_idx.py, per compatibilita'
con l'infrastruttura di lancio a box sul cluster):

    python run_forced_modes.py <min_lon> <max_lon> <min_lat> <max_lat> <box_idx>

Cicla INTERNAMENTE su tutte le 27 finestre per il box di griglia
assegnato a questo job.

Scrittura atomica: ogni output viene scritto prima su un file
temporaneo, poi rinominato al nome definitivo (rename atomico su
Linux) - cosi' merge_forced_modes.py, anche se legge nello stesso
istante, non trova mai un file a meta'.
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

# Categoria/e da processare (deciso in chat, dopo i kill per
# wall-time): 6' argomento opzionale, lista separata da virgola tra
# 'event','trimester','annual'. Default: tutte (comportamento
# originale, un solo job per tutte le 27 finestre).
if len(sys.argv) > 6:
    categories = sys.argv[6].split(",")
else:
    categories = ini.ALL_CATEGORIES
windows_to_run = ini.get_windows_for_categories(categories)
print(f"Categorie richieste: {categories} -> {len(windows_to_run)} finestre")

os.makedirs(ini.work_dir, exist_ok=True)

########################################################################
# Griglia e maschera (una volta sola, riusata per tutte le finestre)
########################################################################

mesh_nemo = nc.Dataset(ini.mesh_mask, "r")
tmask = mesh_nemo.variables["tmask"][0, 0, :, :]
nav_lat_full = mesh_nemo.variables["nav_lat"][:]
nav_lon_full = mesh_nemo.variables["nav_lon"][:]
mesh_nemo.close()

med_mask_full = get_med_mask(nav_lon_full, nav_lat_full, tmask.astype(bool))

########################################################################
# Loop sulle 27 finestre
########################################################################

for kind, key, start_date, end_date in windows_to_run:

    print(f"\n=== Finestra {kind}/{key}: {start_date}-{end_date} ===")

    # Resume: se il file di output per QUESTO box e QUESTA finestra
    # esiste gia' (scritto in modo atomico, quindi garantito completo -
    # vedi nota in testa al file), salta il ricalcolo. Utile per
    # riprendere un box ucciso per wall-time senza rifare da capo le
    # finestre gia' completate prima del kill.
    outfile_check = os.path.join(ini.work_dir, f"forced_modes_{kind}_{key}_{box_idx}.nc")
    if os.path.exists(outfile_check):
        print(f"  Gia' presente ({outfile_check}), salto (resume).")
        continue

    spec_settings = ini.SPECTRUM_SETTINGS[kind]

    try:
        files, missing = build_file_list(
            start_date, end_date, ini.forced_run_daily_template
        )
    except FileNotFoundError as e:
        print(f"  SALTO finestra {key}: {e}")
        continue

    ds = xr.open_mfdataset(files, combine="by_coords", parallel=True)
    ssh_box_lazy = ds[ini.ssh_varname][:, min_lat:max_lat, min_lon:max_lon]

    # FIX PERFORMANCE: carica l'intero box in RAM UNA SOLA VOLTA, invece
    # di forzare una lettura separata per ciascun punto griglia dentro il
    # loop (che con xr.open_mfdataset su centinaia di file giornalieri
    # significa centinaia di letture ripetute per ogni singolo punto -
    # causa identificata dei kill per wall-time su annual/trimester).
    # Il box e' piccolo anche per la finestra piu' pesante (annuale,
    # 8760 timestep): dell'ordine di poche centinaia di MB, ben dentro
    # i 100G gia' richiesti per job - vedi verifica in chat.
    ssh_box = ssh_box_lazy.values * 100.0  # m -> cm (sossheig e' in metri), letto una volta
    ds.close()  # i dati sono gia' in RAM (ssh_box), il file puo' essere chiuso subito

    ny_box = max_lat - min_lat
    nx_box = max_lon - min_lon

    amp_maps = np.full((ini.N_MODES, ny_box, nx_box), np.nan)
    fallback_maps = np.zeros((ini.N_MODES, ny_box, nx_box), dtype=bool)

    all_peak_periods_box = []

    for j in range(ny_box):
        for i in range(nx_box):
            lat_idx = min_lat + j
            lon_idx = min_lon + i

            if tmask[lat_idx, lon_idx] != 1:
                continue  # punto di terra

            ssh_ts = ssh_box[:, j, i]  # gia' in RAM, gia' in cm (nessuna nuova lettura da disco)

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

    # ---- Scrittura output per-box (atomica: temp + rename) ----
    outfile = os.path.join(ini.work_dir, f"forced_modes_{kind}_{key}_{box_idx}.nc")
    tmp_outfile = outfile + f".tmp{os.getpid()}"
    with nc.Dataset(tmp_outfile, "w") as out:
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
    os.replace(tmp_outfile, outfile)  # rename atomico

    # Picchi grezzi del box, per il merge dell'istogramma diagnostico
    peaks_outfile = os.path.join(ini.work_dir, f"peaks_{kind}_{key}_{box_idx}.npz")
    tmp_peaks_outfile = peaks_outfile[:-4] + f".tmp{os.getpid()}.npz"  # mantiene .npz finale
    if len(all_peak_periods_box) > 0:
        np.savez(tmp_peaks_outfile,
                 peak_periods=np.concatenate(all_peak_periods_box))
    else:
        np.savez(tmp_peaks_outfile, peak_periods=np.array([]))
    os.replace(tmp_peaks_outfile, peaks_outfile)

    print(f"  Scritto: {outfile}")
    print(f"  Scritto: {peaks_outfile}")

print("\nFatto.")
