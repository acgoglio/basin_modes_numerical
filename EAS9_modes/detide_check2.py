"""
Verifica del detiding sui punti di idx_pt.coo (righe non commentate),
da lanciare PRIMA del run completo con flag_detide=True.

Per ogni punto e per ogni finestra in CHECK_WINDOWS:
  - figura con (a) serie temporale dei primi TS_DAYS giorni: sossheig
    originale, marea ricostruita, segnale detidato; (b) spettro di
    ampiezza prima/dopo, calcolato con LA STESSA funzione e gli STESSI
    parametri della pipeline per quella categoria di finestra, con le
    bande dei modi (T_ref +/- tol) e i periodi dei costituenti;
  - riga di riepilogo (CSV + stampa): varianza nelle bande diurna e
    semidiurna prima/dopo e ampiezze dei 10 modi prima/dopo, estratte
    come nella pipeline.

Come riferimento, in entrambi i pannelli e nel riepilogo c'e' anche la
run senza maree (EAS9_minr_nt) negli stessi punti e finestre: se il
detiding funziona, il segnale detidato dovrebbe somigliarle (a meno
dell'interazione marea-surge). Se la run NT non e' disponibile per una
finestra, la curva viene omessa con un messaggio.

Le finestre di default sono una DENTRO il periodo del fit armonico
(OND_2022, fit su lug-dic 2022) e una FUORI (JFM_2020), per verificare
che le costanti del 2022 funzionino anche sugli altri anni.

Il detiding qui e' sempre applicato (e' il confronto con/senza),
indipendentemente da flag_detide; richiede EXPERIMENT = "EAS9_simu".
Gli output vanno in <work_dir>/check_detide/ (non tocca altro).

Uso:  python check_detide.py
"""

import os
import csv
import numpy as np
import pandas as pd
import xarray as xr
import netCDF4 as nc
import matplotlib
matplotlib.use("Agg")
import matplotlib.ticker
import matplotlib.pyplot as plt

# Font piu' grandi in tutte le figure
plt.rcParams.update({
    "font.size": 14,
    "axes.titlesize": 15,
    "axes.labelsize": 14,
    "xtick.labelsize": 13,
    "ytick.labelsize": 13,
    "legend.fontsize": 12,
})

import forced_windows_ini as ini
from build_file_list import build_file_list
from f_forced_ampspt import (
    compute_point_spectrum,
    find_all_peaks,
    extract_reference_mode_amplitudes,
)
from f_detide import (
    get_tidal_frequencies,
    load_amppha_box,
    times_to_ttide,
    synth_tide_point,
)

########################################################################
# Parametri della verifica
########################################################################

COO_FILE = "idx_pt.coo"                   # formato: "i j nome", i = x, j = y; '#' = commento
CHECK_WINDOWS = ["OND_2022", "JFM_2020"]  # chiavi di ini.WINDOWS
TS_DAYS = 7                               # giorni di serie temporale nel pannello (a)
PERIOD_RANGE_PLOT = (4.0, 40.0)           # [h] asse x dello spettro

# Bande [h] per la varianza prima/dopo. Nota: la banda diurna contiene
# anche il modo 0 (27 +/- 2 h), quindi dopo il detiding la varianza non
# va a zero: ci resta l'energia della sessa.
VAR_BANDS = {"diurna": (22.0, 28.0), "semidiurna": (11.0, 13.0)}

OUT_DIR = os.path.join(ini.work_dir, "check_detide")

# Run senza maree usata come riferimento
NT_TEMPLATE = ini._EXPERIMENTS["EAS9_minr_nt"]["template"]

assert ini.EXPERIMENT == "EAS9_simu", \
    f"check_detide ha senso solo sulla run con maree (EXPERIMENT={ini.EXPERIMENT})"


########################################################################
# Funzioni di supporto
########################################################################

def read_coo(path):
    """Punti non commentati di idx_pt.coo -> lista di (i, j, nome)."""
    pts = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            pts.append((int(parts[0]), int(parts[1]), " ".join(parts[2:])))
    return pts


def band_variance(x, dt_h, T_min, T_max):
    """
    Varianza del segnale nella banda di periodi [T_min, T_max] h, da
    periodogramma (Parseval) della serie senza NaN e senza trend lineare.
    """
    x = np.asarray(x, dtype=np.float64)
    x = x[np.isfinite(x)]
    n = len(x)
    t = np.arange(n)
    x = x - np.polyval(np.polyfit(t, x, 1), t)
    X = np.fft.rfft(x)
    f = np.fft.rfftfreq(n, d=dt_h)  # [cph]
    with np.errstate(divide="ignore"):
        T = 1.0 / f
    sel = (T >= T_min) & (T <= T_max)
    # rfft: i bin interni contano due volte (frequenze negative)
    w = np.where((np.arange(len(f)) == 0) | ((n % 2 == 0) & (np.arange(len(f)) == len(f) - 1)), 1.0, 2.0)
    return float(np.sum(w[sel] * np.abs(X[sel]) ** 2) / n ** 2)


def read_points(template, d0, d1, pts):
    """sossheig [cm] (time, pt) e tempi (datetime) nei punti, per la finestra."""
    files, missing = build_file_list(d0, d1, template)
    ds = xr.open_mfdataset(files, combine="by_coords", parallel=True)
    i_idx = xr.DataArray([p[0] for p in pts], dims="pt")
    j_idx = xr.DataArray([p[1] for p in pts], dims="pt")
    vals = ds[ini.ssh_varname].isel(y=j_idx, x=i_idx).values * 100.0
    times = pd.to_datetime(ds["time_counter"].values).to_pydatetime()
    ds.close()
    return vals, times


def rms_diff(a, b):
    """rms della differenza tra due serie, tolte le rispettive medie (cm)."""
    ok = np.isfinite(a) & np.isfinite(b)
    if not ok.any():
        return np.nan
    return float(np.sqrt(np.mean(((a[ok] - a[ok].mean()) - (b[ok] - b[ok].mean())) ** 2)))


def spectrum_and_modes(ts, spec_settings):
    """Spettro e ampiezze dei modi con le stesse funzioni della pipeline."""
    freq_pos, periods, amplitudes = compute_point_spectrum(
        ts, ini.dt,
        flag_hanning=ini.flag_hanning,
        flag_nfft=ini.flag_nfft,
        N_fft=ini.N_fft,
        flag_segmented_spectrum=spec_settings["flag_segmented_spectrum"],
        segment_len_days=spec_settings["segment_len_days"],
        flag_filter=ini.flag_filter,
        th_filter=ini.th_filter,
    )
    if periods is None:
        return None, None, np.full(ini.N_MODES, np.nan)
    peak_periods, peak_amps = find_all_peaks(periods, amplitudes)
    mode_amps, _ = extract_reference_mode_amplitudes(
        periods, amplitudes, peak_periods, peak_amps,
        ini.REFERENCE_MODES, flag_use_fallback=ini.flag_use_fallback,
    )
    return periods, amplitudes, np.asarray(mode_amps, dtype=np.float64)


########################################################################
# Preparazione
########################################################################

os.makedirs(OUT_DIR, exist_ok=True)

points = read_coo(COO_FILE)
print(f"Punti da {COO_FILE}: {len(points)}")

mesh = nc.Dataset(ini.mesh_mask)
tmask = mesh.variables["tmask"][0, 0, :, :]
nav_lat = np.array(mesh.variables["nav_lat"][:])
nav_lon = np.array(mesh.variables["nav_lon"][:])
mesh.close()
ny, nx = tmask.shape

tidal_names_raw, tidal_freq, tidal_names = get_tidal_frequencies(ini.TIDE_CONSTITUENTS)
tidal_periods = {n: 1.0 / f for n, f in zip(tidal_names, tidal_freq)}

# Punti validi: dentro la griglia, mare, ampiezze/fasi finite
valid_pts = []
for (i, j, name) in points:
    if not (0 <= i < nx and 0 <= j < ny):
        print(f"  SALTO {name} ({i},{j}): fuori griglia ({nx} x {ny})")
        continue
    if tmask[j, i] != 1:
        print(f"  SALTO {name} ({i},{j}): punto di terra")
        continue
    amp_pt, pha_pt = load_amppha_box(ini.AMPPHA_FILE, ini.TIDE_CONSTITUENTS, j, j + 1, i, i + 1)
    if not all(np.isfinite(amp_pt[c][0, 0]) and np.isfinite(pha_pt[c][0, 0]) for c in tidal_names):
        print(f"  SALTO {name} ({i},{j}): ampiezze/fasi non valide")
        continue
    print(f"  {name:18s} i={i:4d} j={j:3d}  lon={nav_lon[j, i]:7.3f}  lat={nav_lat[j, i]:6.3f}  "
          f"M2={100 * amp_pt['M2'][0, 0]:.1f} cm  K1={100 * amp_pt['K1'][0, 0]:.1f} cm")
    valid_pts.append((i, j, name, amp_pt, pha_pt))

if not valid_pts:
    raise SystemExit("Nessun punto valido.")

windows = {key: (kind, d0, d1) for (kind, key, d0, d1) in ini.WINDOWS}
mode_labels = [f"{T_ref}h" for (_, T_ref, _) in ini.REFERENCE_MODES]

rows = []

########################################################################
# Loop sulle finestre
########################################################################

for wkey in CHECK_WINDOWS:
    if wkey not in windows:
        print(f"\nSALTO finestra {wkey}: non presente in ini.WINDOWS "
              f"(TEST_SINGLE_WINDOW_KEY = {getattr(ini, 'TEST_SINGLE_WINDOW_KEY', None)!r}: "
              f"se e' attivo, ini.WINDOWS contiene solo quella finestra)")
        continue
    kind, d0, d1 = windows[wkey]
    spec_settings = ini.SPECTRUM_SETTINGS[kind]
    print(f"\n=== Finestra {kind}/{wkey}: {d0}-{d1} ===")

    ssh_pts, times_dt = read_points(ini.forced_run_daily_template, d0, d1, valid_pts)

    # Run senza maree, stessi punti e finestra (solo se stessi tempi)
    ssh_nt = None
    try:
        ssh_nt, times_nt = read_points(NT_TEMPLATE, d0, d1, valid_pts)
        if len(times_nt) != len(times_dt) or any(a != b for a, b in zip(times_nt, times_dt)):
            print(f"  ATTENZIONE: tempi della run NT diversi da quelli della run con maree "
                  f"(NT: {len(times_nt)} passi, da {times_nt[0]}; maree: {len(times_dt)} passi, "
                  f"da {times_dt[0]}): curva NT omessa.")
            ssh_nt = None
    except FileNotFoundError as e:
        print(f"  Run NT non disponibile per {wkey}: {e} -> curva NT omessa.")

    t_num = times_to_ttide(times_dt)  # si ferma se non tutti a HH:30
    if len(t_num) != ssh_pts.shape[0]:
        raise ValueError(f"time_counter ({len(t_num)}) e sossheig ({ssh_pts.shape[0]}) diversi")

    n_ts = min(len(times_dt), TS_DAYS * 24)

    for p, (i, j, name, amp_pt, pha_pt) in enumerate(valid_pts):
        orig = ssh_pts[:, p]
        tide = 100.0 * synth_tide_point(t_num, tidal_names_raw, tidal_freq, tidal_names,
                                        amp_pt, pha_pt, 0, 0, float(nav_lat[j, i]))
        det = orig - tide

        per_o, amp_o, modes_o = spectrum_and_modes(orig, spec_settings)
        per_d, amp_d, modes_d = spectrum_and_modes(det, spec_settings)
        if ssh_nt is not None:
            nt = ssh_nt[:, p]
            per_n, amp_n, modes_n = spectrum_and_modes(nt, spec_settings)
        else:
            nt, per_n, amp_n, modes_n = None, None, None, np.full(ini.N_MODES, np.nan)

        row = {"window": wkey, "point": name, "i": i, "j": j}
        for bname, (tmin, tmax) in VAR_BANDS.items():
            vo = band_variance(orig, ini.dt / 3600.0, tmin, tmax)
            vd = band_variance(det, ini.dt / 3600.0, tmin, tmax)
            row[f"var_{bname}_prima_cm2"] = round(vo, 3)
            row[f"var_{bname}_dopo_cm2"] = round(vd, 3)
            row[f"rid_{bname}_%"] = round(100.0 * (1.0 - vd / vo), 1) if vo > 0 else np.nan
            row[f"var_{bname}_NT_cm2"] = (round(band_variance(nt, ini.dt / 3600.0, tmin, tmax), 3)
                                          if nt is not None else np.nan)
        # scarto rms (medie tolte) rispetto alla run NT, prima e dopo il detiding
        row["rms_orig_meno_NT_cm"] = round(rms_diff(orig, nt), 2) if nt is not None else np.nan
        row["rms_detid_meno_NT_cm"] = round(rms_diff(det, nt), 2) if nt is not None else np.nan
        for lab, ao, ad, an in zip(mode_labels, modes_o, modes_d, modes_n):
            row[f"A_{lab}_prima"] = round(float(ao), 2) if np.isfinite(ao) else np.nan
            row[f"A_{lab}_dopo"] = round(float(ad), 2) if np.isfinite(ad) else np.nan
            row[f"A_{lab}_NT"] = round(float(an), 2) if np.isfinite(an) else np.nan
        rows.append(row)

        # ---- Figura ----
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 10))

        tt = times_dt[:n_ts]
        ax1.plot(tt, orig[:n_ts], color="0.5", lw=1.0, label="Original sossheig")
        ax1.plot(tt, tide[:n_ts], color="tab:blue", lw=0.9, label="Reconstructed tide")
        ax1.plot(tt, det[:n_ts], color="tab:red", lw=1.0, label="Detided")
        if nt is not None:
            ax1.plot(tt, nt[:n_ts], color="tab:green", lw=1.0, ls="--", label="No-tide run (NT)")
        ax1.set_ylabel("Sea level [cm]")
        ax1.set_title(f"{name} (i={i}, j={j}; {nav_lon[j, i]:.2f}E, {nav_lat[j, i]:.2f}N) "
                      f"- {wkey}: first {TS_DAYS} days")
        ax1.legend(loc="upper right")
        ax1.grid(ls="--", lw=0.4, alpha=0.5)

        for m, (_, T_ref, tol) in enumerate(ini.REFERENCE_MODES):
            ax2.axvspan(T_ref - tol, T_ref + tol, color="gold", alpha=0.25, lw=0)
            # etichette alternate su due altezze per non sovrapporre bande vicine (7.2/6.9)
            ax2.text(T_ref, 1.01 + 0.05 * (m % 2), f"{T_ref:g}", transform=ax2.get_xaxis_transform(),
                     ha="center", va="bottom", fontsize=11)
        # costituenti vicini (P1/K1, M2/S2/K2): etichette alternate su due altezze
        for k, (cname, Tc) in enumerate(sorted(tidal_periods.items(), key=lambda kv: kv[1])):
            ax2.axvline(Tc, color="tab:blue", ls=":", lw=0.8)
            ax2.text(Tc, 0.02 + 0.13 * (k % 2), cname, transform=ax2.get_xaxis_transform(),
                     rotation=90, ha="right", va="bottom", fontsize=11, color="tab:blue")
        if per_o is not None:
            ax2.plot(per_o, amp_o, color="0.5", lw=1.0, label="Original")
        if per_d is not None:
            ax2.plot(per_d, amp_d, color="tab:red", lw=1.0, label="Detided")
        if per_n is not None:
            ax2.plot(per_n, amp_n, color="tab:green", lw=1.0, ls="--", label="No-tide run (NT)")
        ax2.set_xscale("log")
        ax2.set_yscale("log")
        ax2.set_xlim(*PERIOD_RANGE_PLOT)
        ax2.invert_xaxis()
        ticks = [t for t in (4, 5, 6, 8, 10, 12, 15, 20, 25, 30, 40)
                 if PERIOD_RANGE_PLOT[0] <= t <= PERIOD_RANGE_PLOT[1]]
        ax2.set_xticks(ticks)
        ax2.set_xticklabels([str(t) for t in ticks])
        ax2.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax2.set_xlabel("Period [h]")
        ax2.set_ylabel("Amplitude [cm]")
        ax2.set_title(f"Amplitude spectrum ({kind} settings, as in the pipeline): "
                      f"mode bands in yellow, tidal constituents in blue", pad=38)
        ax2.legend(loc="upper right")
        ax2.grid(ls="--", lw=0.4, alpha=0.5, which="both")

        plt.tight_layout()
        fname = os.path.join(OUT_DIR, f"check_detide_{wkey}_{name}_{i}_{j}.png")
        plt.savefig(fname, dpi=130)
        plt.close(fig)

        print(f"  {name:18s} riduzione varianza diurna {row['rid_diurna_%']:6.1f}%  "
              f"semidiurna {row['rid_semidiurna_%']:6.1f}%  |  "
              f"A27h {row.get('A_27.0h_prima')} -> {row.get('A_27.0h_dopo')}  "
              f"A11.9h {row.get('A_11.9h_prima')} -> {row.get('A_11.9h_dopo')}  |  "
              f"rms vs NT: {row['rms_orig_meno_NT_cm']} -> {row['rms_detid_meno_NT_cm']} cm")

########################################################################
# Riepilogo
########################################################################

if rows:
    csv_path = os.path.join(OUT_DIR, "check_detide_summary.csv")
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"\nRiepilogo: {csv_path}")
    print(f"Figure:    {OUT_DIR}")
