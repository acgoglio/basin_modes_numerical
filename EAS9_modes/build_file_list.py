"""
Costruisce la lista ESPLICITA dei file netCDF giornalieri necessari per
una finestra [start_date, end_date], secondo il tuo consiglio di non
aprire piu' file di quelli strettamente necessari (niente glob su '20*').

Assunzione DA VERIFICARE: un file per ogni giorno del calendario tra
start_date e end_date inclusi (nessun giorno mancante nell'archivio).
Se ci sono buchi noti nella run forzata, va gestito qui (es. con un
controllo os.path.exists + warning invece di errore).
"""

import os
from datetime import datetime, timedelta


def daterange(start_date, end_date):
    """start_date, end_date: stringhe 'YYYYMMDD' (inclusive)"""
    d0 = datetime.strptime(start_date, "%Y%m%d")
    d1 = datetime.strptime(end_date, "%Y%m%d")
    n_days = (d1 - d0).days
    if n_days < 0:
        raise ValueError(f"end_date {end_date} precede start_date {start_date}")
    for i in range(n_days + 1):
        yield (d0 + timedelta(days=i)).strftime("%Y%m%d")


def build_file_list(start_date, end_date, daily_template, verbose=True):
    """
    Ritorna la lista dei path esistenti tra start_date e end_date.

    daily_template: stringa con placeholder {yyyymmdd}, es.
        "/work/.../{yyyymmdd}/model/medfs-eas9_1h_{yyyymmdd}_2D_grid_T.nc"
    """
    files = []
    missing = []
    for yyyymmdd in daterange(start_date, end_date):
        f = daily_template.format(yyyymmdd=yyyymmdd)
        if os.path.exists(f):
            files.append(f)
        else:
            missing.append(f)

    if verbose:
        print(f"[build_file_list] {start_date}-{end_date}: "
              f"{len(files)} file trovati, {len(missing)} mancanti")
        if missing:
            print("  Primi file mancanti:")
            for m in missing[:5]:
                print("   ", m)

    if len(files) == 0:
        raise FileNotFoundError(
            f"Nessun file trovato per la finestra {start_date}-{end_date} "
            f"con template {daily_template}"
        )

    return files, missing
