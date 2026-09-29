"""
Lettura rapida di forced_modes_results.npz per controllare, a run in
corso o finito, come stanno venendo i risultati.

Uso:
    python peek_results.py
"""

import numpy as np
import forced_windows_ini as ini
import os

RESULTS_FILE = os.path.join(ini.work_dir, "merged", "forced_modes_results.npz")

if not os.path.exists(RESULTS_FILE):
    print(f"Nessun risultato ancora disponibile: {RESULTS_FILE} non esiste.")
    print("Lancia prima almeno una volta: ./run_forced_BM.sh merge")
    raise SystemExit(1)

npz = np.load(RESULTS_FILE, allow_pickle=True)
data = npz["data"]        # (10 modi, 4 trimestri, 4 anni)
events = npz["events"]    # (10 modi, 7 eventi)
annual = npz["annual"]    # (10 modi, 4 anni)
modes = list(npz["modes"])
trimesters = list(npz["trimesters"])
years = list(npz["years"])
events_names = list(npz["events_names"])

n_modes = len(modes)


def fmt(v):
    return "   .  " if np.isnan(v) else f"{v:6.2f}"


print("=" * 70)
print("EVENTI (99th percentile, cm) - una colonna per evento")
print("=" * 70)
header = f"{'Modo':30s} " + " ".join(f"{n[:10]:>10s}" for n in events_names)
print(header)
for m in range(n_modes):
    row = f"{modes[m]:30s} " + " ".join(f"{fmt(events[m, e]):>10s}" for e in range(len(events_names)))
    print(row)
n_done_events = np.sum(~np.isnan(events).all(axis=0))
print(f"\n-> {n_done_events}/{len(events_names)} eventi con almeno un valore calcolato")

print()
print("=" * 70)
print("TRIMESTRI (99th percentile, cm) - media sugli anni gia' disponibili")
print("=" * 70)
with np.errstate(invalid="ignore"):
    data_mean_partial = np.nanmean(data, axis=2)  # ignora i NaN (anni non ancora fatti)
header = f"{'Modo':30s} " + " ".join(f"{t:>8s}" for t in trimesters)
print(header)
for m in range(n_modes):
    row = f"{modes[m]:30s} " + " ".join(f"{fmt(data_mean_partial[m, t]):>8s}" for t in range(4))
    print(row)

n_windows_total = 4 * 4 + 4
n_windows_done = int(np.sum(~np.isnan(data).all(axis=0)) + np.sum(~np.isnan(annual).all(axis=0)))
print(f"\n-> {n_windows_done}/{n_windows_total} finestre trimestre/anno con almeno un valore calcolato")

print()
print("=" * 70)
print("RUN ANNUALI (99th percentile, cm) - riferimento/QA")
print("=" * 70)
header = f"{'Modo':30s} " + " ".join(f"{y:>8d}" for y in years)
print(header)
for m in range(n_modes):
    row = f"{modes[m]:30s} " + " ".join(f"{fmt(annual[m, y]):>8s}" for y in range(4))
    print(row)

print()
print("NB: '.' = ancora nessun dato per quella cella (finestra/box non completati).")
print("Rilancia './run_forced_BM.sh merge' per aggiornare questo snapshot con l'avanzamento piu' recente.")
