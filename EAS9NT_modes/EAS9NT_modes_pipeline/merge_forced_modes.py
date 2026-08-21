"""
Merge dei box per ciascuna delle 27 finestre. Da lanciare DOPO che
run_forced_modes.py e' stato eseguito su tutti i 108 box.

Per ciascuna finestra:
  1. Fonde le mappe (mode,y,x) dei 108 box in un'unica mappa full-domain
  2. Esclude i 4 box della Baia di Biscaglia (73,74,91,92), stesso
     schema di merge_amp_idx.py
  3. Calcola il 99th percentile per modo sui punti mare rimanenti
  4. Concatena i picchi grezzi di tutti i box -> istogramma diagnostico
     (PNG, stile Fig. 3.5) per validazione visiva
  5. Calcola % di used_fallback per modo (QA)

Alla fine assembla:
  - data[10 modi, 4 trimestri, 4 anni]  (NON ancora mediato sugli anni -
    la media la fa tab_3month_extrev.py, come richiesto)
  - events[10 modi, 7 eventi]
  - annual[10 modi, 4 anni]  (extra, non usato da tab_3month_extrev.py
    ma utile per QA/confronto - i periodi di riferimento erano gia'
    presi da qui in origine)

e li scrive in forced_modes_results.npz, pronto per essere caricato da
una versione modificata di tab_3month_extrev.py.
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import forced_windows_ini as ini
import netCDF4 as nc

# Schema a box: unica fonte di verita' in forced_windows_ini.py
x_edges = ini.BOX_X_EDGES
y_edges = ini.BOX_Y_EDGES
N_ROWS, N_COLS = ini.BOX_N_ROWS, ini.BOX_N_COLS
EXCLUDE_BOXES = ini.BOX_EXCLUDE

NX_FULL = x_edges[-1] - x_edges[0]
NY_FULL = y_edges[-1] - y_edges[0]


def box_bounds(idx):
    idx0 = idx - 1
    r = idx0 // N_COLS
    c = idx0 % N_COLS
    return x_edges[c], x_edges[c + 1], y_edges[r], y_edges[r + 1]


def merge_window(kind, key):
    """Fonde i 108 box per una finestra. Ritorna (amp_full, fallback_full)
    con shape (N_MODES, NY_FULL, NX_FULL), o (None, None) se nessun box
    trovato."""
    amp_full = np.full((ini.N_MODES, NY_FULL, NX_FULL), np.nan)
    fb_full = np.zeros((ini.N_MODES, NY_FULL, NX_FULL), dtype=bool)
    n_found = 0

    for box_idx in range(1, N_ROWS * N_COLS + 1):
        fpath = os.path.join(ini.work_dir, f"forced_modes_{kind}_{key}_{box_idx}.nc")
        if not os.path.exists(fpath):
            continue
        n_found += 1
        x0, x1, y0, y1 = box_bounds(box_idx)
        with nc.Dataset(fpath, "r") as ds:
            amp_full[:, y0:y1, x0:x1] = ds.variables["amplitude"][:]
            fb_full[:, y0:y1, x0:x1] = ds.variables["used_fallback"][:].astype(bool)

    if n_found == 0:
        print(f"  ATTENZIONE: nessun box trovato per {kind}/{key}")
        return None, None

    if n_found < N_ROWS * N_COLS:
        print(f"  ATTENZIONE: solo {n_found}/{N_ROWS*N_COLS} box trovati per {kind}/{key}")

    # Escludi Baia di Biscaglia
    for b in EXCLUDE_BOXES:
        x0, x1, y0, y1 = box_bounds(b)
        amp_full[:, y0:y1, x0:x1] = np.nan
        fb_full[:, y0:y1, x0:x1] = False

    return amp_full, fb_full


def merge_peaks(kind, key):
    """Concatena i picchi grezzi di tutti i box per una finestra."""
    all_peaks = []
    for box_idx in range(1, N_ROWS * N_COLS + 1):
        fpath = os.path.join(ini.work_dir, f"peaks_{kind}_{key}_{box_idx}.npz")
        if not os.path.exists(fpath):
            continue
        d = np.load(fpath)
        if d["peak_periods"].size > 0:
            all_peaks.append(d["peak_periods"])
    if len(all_peaks) == 0:
        return np.array([])
    return np.concatenate(all_peaks)


def diagnostic_histogram(peak_periods, kind, key, outdir):
    """Istogramma di copertura dei picchi, con le bande di riferimento
    T_i+/-tol_i sovrapposte - permette il controllo visivo che i cluster
    di picchi coincidano con i modi attesi (stile Fig. 3.5)."""
    fig, ax = plt.subplots(figsize=(10, 4))
    if len(peak_periods) > 0:
        ax.hist(peak_periods, bins=200, range=(3, 30), color="0.6")
    for label, T_ref, tol in ini.REFERENCE_MODES:
        ax.axvspan(T_ref - tol, T_ref + tol, alpha=0.25, color="tab:orange")
        ax.axvline(T_ref, color="tab:red", lw=0.8)
    ax.set_xlabel("Periodo (h)")
    ax.set_ylabel("Numero di picchi (tutti i box)")
    ax.set_title(f"Istogramma diagnostico - {kind}/{key}")
    fig.tight_layout()
    outfile = os.path.join(outdir, f"hist_diag_{kind}_{key}.png")
    fig.savefig(outfile, dpi=150)
    plt.close(fig)
    return outfile


########################################################################

def main():
    outdir = os.path.join(ini.work_dir, "merged")
    os.makedirs(outdir, exist_ok=True)

    # Contenitori finali
    data = np.full((ini.N_MODES, 4, 4), np.nan)     # modo, trimestre, anno
    events = np.full((ini.N_MODES, 7), np.nan)       # modo, evento
    annual = np.full((ini.N_MODES, 4), np.nan)       # modo, anno (extra/QA)
    fallback_pct = {}  # (kind,key) -> array (N_MODES,) % fallback

    year_index = {y: i for i, y in enumerate(ini.YEARS)}
    trimester_index = {t: i for i, t in enumerate(ini.TRIMESTER_ORDER)}
    event_index = {name: i for i, name in enumerate(ini.EVENT_ORDER)}

    for kind, key, start_date, end_date in ini.WINDOWS:
        print(f"Merging {kind}/{key} ...")
        amp_full, fb_full = merge_window(kind, key)
        if amp_full is None:
            continue

        # 99th percentile per modo sui punti mare validi (non-NaN)
        p99 = np.full(ini.N_MODES, np.nan)
        fb_pct = np.full(ini.N_MODES, np.nan)
        for m in range(ini.N_MODES):
            vals = amp_full[m][~np.isnan(amp_full[m])]
            if vals.size > 0:
                p99[m] = np.percentile(vals, 99)
            fb_vals = fb_full[m][~np.isnan(amp_full[m])]
            if fb_vals.size > 0:
                fb_pct[m] = 100.0 * np.mean(fb_vals)
        fallback_pct[(kind, key)] = fb_pct

        # Istogramma diagnostico
        peaks = merge_peaks(kind, key)
        diagnostic_histogram(peaks, kind, key, outdir)

        # Smista nel contenitore finale giusto
        if kind == "annual":
            y = int(key)
            annual[:, year_index[y]] = p99
        elif kind == "trimester":
            tname, y = key.split("_")
            y = int(y)
            data[:, trimester_index[tname], year_index[y]] = p99
        elif kind == "event":
            events[:, event_index[key]] = p99

    # Salva risultati
    outfile = os.path.join(outdir, "forced_modes_results.npz")
    np.savez(
        outfile,
        data=data,
        events=events,
        annual=annual,
        modes=np.array([label for label, _, _ in ini.REFERENCE_MODES]),
        trimesters=np.array(ini.TRIMESTER_ORDER),
        years=np.array(ini.YEARS),
        events_names=np.array(ini.EVENT_ORDER),
    )
    print(f"\nScritto: {outfile}")

    # Riepilogo % fallback (QA) - stampato, non salvato in automatico:
    # decidi tu come/se usarlo (soglia di esclusione, nota a piè di
    # tabella, o nient'altro)
    print("\n=== Riepilogo %% fallback per finestra/modo (QA) ===")
    for (kind, key), fb_pct in fallback_pct.items():
        high = [(ini.REFERENCE_MODES[m][0], fb_pct[m])
                for m in range(ini.N_MODES) if fb_pct[m] > 50]
        if high:
            print(f"  {kind}/{key}: fallback >50% per: "
                  + ", ".join(f"{lbl} ({v:.0f}%)" for lbl, v in high))


if __name__ == "__main__":
    main()
