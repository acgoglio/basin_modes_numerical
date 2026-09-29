"""
Merge dei box per ciascuna delle 27 finestre. Da lanciare DOPO che
run_forced_modes.py e' stato eseguito su tutti i 108 box (o anche a
run parziale, per un controllo di avanzamento).

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
    ma utile per QA/confronto)

e li scrive in forced_modes_results.npz.
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import forced_windows_ini as ini
import netCDF4 as nc

# Schema a box: unica fonte di verita' in forced_windows_ini.py.
# Il dominio fisico totale (NX_FULL/NY_FULL) e' identico in tutti gli
# schemi (coarse/fine coprono la stessa area, solo suddivisa in modo
# diverso) - lo prendo dallo schema coarse come riferimento.
_ref_scheme = ini.BOX_SCHEMES["coarse"]
NX_FULL = _ref_scheme["x_edges"][-1] - _ref_scheme["x_edges"][0]
NY_FULL = _ref_scheme["y_edges"][-1] - _ref_scheme["y_edges"][0]

# Offset del dominio: box_bounds ritorna coordinate ASSOLUTE della
# griglia (es. 300-1306 in x), ma amp_full/fb_full sono array
# indicizzati da 0 (dimensione NX_FULL/NY_FULL) - va sottratto
# l'offset prima di usarle come indici, altrimenti si scrive fuori
# posizione (o fuori dai limiti dell'array, per le colonne piu' a
# destra) invece che nella cella corretta.
DOMAIN_X0 = _ref_scheme["x_edges"][0]
DOMAIN_Y0 = _ref_scheme["y_edges"][0]


def box_bounds(scheme, idx):
    idx0 = idx - 1
    r = idx0 // scheme["n_cols"]
    c = idx0 % scheme["n_cols"]
    return (scheme["x_edges"][c], scheme["x_edges"][c + 1],
            scheme["y_edges"][r], scheme["y_edges"][r + 1])


def merge_window(kind, key):
    """Fonde i box (108 o 432, a seconda della categoria 'kind') per
    una finestra. Ritorna (amp_full, fallback_full) con shape
    (N_MODES, NY_FULL, NX_FULL), o (None, None) se nessun box trovato."""
    scheme = ini.get_box_scheme_for_category(kind)
    n_boxes = scheme["n_rows"] * scheme["n_cols"]

    amp_full = np.full((ini.N_MODES, NY_FULL, NX_FULL), np.nan)
    fb_full = np.zeros((ini.N_MODES, NY_FULL, NX_FULL), dtype=bool)
    n_found = 0

    for box_idx in range(1, n_boxes + 1):
        fpath = os.path.join(ini.work_dir, f"forced_modes_{kind}_{key}_{box_idx}.nc")
        if not os.path.exists(fpath):
            continue
        x0, x1, y0, y1 = box_bounds(scheme, box_idx)
        # Indici LOCALI (relativi all'origine del dominio), non assoluti -
        # vedi nota sopra DOMAIN_X0/DOMAIN_Y0.
        lx0, lx1 = x0 - DOMAIN_X0, x1 - DOMAIN_X0
        ly0, ly1 = y0 - DOMAIN_Y0, y1 - DOMAIN_Y0
        try:
            with nc.Dataset(fpath, "r") as ds:
                amp_var = ds.variables["amplitude"][:]
                fb_var = ds.variables["used_fallback"][:]
                expected_shape = (ini.N_MODES, ly1 - ly0, lx1 - lx0)
                if amp_var.shape != expected_shape:
                    print(f"  ATTENZIONE: box {box_idx} ({kind}/{key}) ha shape "
                          f"{amp_var.shape} ma lo schema box si aspetta {expected_shape} "
                          f"(file: {fpath}) - box SALTATO, non incluso nel merge.")
                    continue
                amp_full[:, ly0:ly1, lx0:lx1] = amp_var
                fb_full[:, ly0:ly1, lx0:lx1] = fb_var.astype(bool)
        except (OSError, KeyError, RuntimeError, ValueError) as e:
            print(f"  Box {box_idx} ({kind}/{key}) non leggibile o inconsistente, salto: {e}")
            continue
        n_found += 1

    if n_found == 0:
        print(f"  ATTENZIONE: nessun box trovato per {kind}/{key}")
        return None, None

    if n_found < n_boxes:
        print(f"  ATTENZIONE: solo {n_found}/{n_boxes} box trovati per {kind}/{key}")

    # Escludi Baia di Biscaglia (rettangolo fisico identico in
    # qualunque schema, solo gli indici di box cambiano) - stesso
    # offset locale applicato sopra.
    for b in scheme["exclude"]:
        x0, x1, y0, y1 = box_bounds(scheme, b)
        lx0, lx1 = x0 - DOMAIN_X0, x1 - DOMAIN_X0
        ly0, ly1 = y0 - DOMAIN_Y0, y1 - DOMAIN_Y0
        amp_full[:, ly0:ly1, lx0:lx1] = np.nan
        fb_full[:, ly0:ly1, lx0:lx1] = False

    return amp_full, fb_full


def merge_peaks(kind, key):
    """Concatena i picchi grezzi di tutti i box per una finestra."""
    scheme = ini.get_box_scheme_for_category(kind)
    n_boxes = scheme["n_rows"] * scheme["n_cols"]
    all_peaks = []
    for box_idx in range(1, n_boxes + 1):
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
    T_i+/-tol_i sovrapposte."""
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

    data = np.full((ini.N_MODES, 4, 4), np.nan)     # modo, trimestre, anno
    events = np.full((ini.N_MODES, 7), np.nan)       # modo, evento
    annual = np.full((ini.N_MODES, 4), np.nan)       # modo, anno (extra/QA)
    fallback_pct = {}

    year_index = {y: i for i, y in enumerate(ini.YEARS)}
    trimester_index = {t: i for i, t in enumerate(ini.TRIMESTER_ORDER)}
    event_index = {name: i for i, name in enumerate(ini.EVENT_ORDER)}

    for kind, key, start_date, end_date in ini.WINDOWS:
        print(f"Merging {kind}/{key} ...")
        amp_full, fb_full = merge_window(kind, key)
        if amp_full is None:
            continue

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

        peaks = merge_peaks(kind, key)
        diagnostic_histogram(peaks, kind, key, outdir)

        if kind == "annual":
            y = int(key)
            annual[:, year_index[y]] = p99
        elif kind == "trimester":
            tname, y = key.split("_")
            y = int(y)
            data[:, trimester_index[tname], year_index[y]] = p99
        elif kind == "event":
            events[:, event_index[key]] = p99

        print(f"  99th percentile per modo ({key}):")
        for m, (label, T_ref, tol) in enumerate(ini.REFERENCE_MODES):
            val = "n/d" if np.isnan(p99[m]) else f"{p99[m]:.2f} cm"
            print(f"    {label:30s} (T={T_ref:>5.2f}h): {val}")

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

    print("\n=== Riepilogo %% fallback per finestra/modo (QA) ===")
    for (kind, key), fb_pct in fallback_pct.items():
        high = [(ini.REFERENCE_MODES[m][0], fb_pct[m])
                for m in range(ini.N_MODES) if fb_pct[m] > 50]
        if high:
            print(f"  {kind}/{key}: fallback >50% per: "
                  + ", ".join(f"{lbl} ({v:.0f}%)" for lbl, v in high))


if __name__ == "__main__":
    main()
