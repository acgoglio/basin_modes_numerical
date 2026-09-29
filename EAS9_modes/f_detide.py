"""
Detiding armonico per la pipeline dei modi sulla run con maree
(EAS9-simu/EXP00). Attivato da forced_windows_ini.flag_detide.

NON ricalcola ampiezze e fasi: le legge dal file di analisi armonica
gia' esistente (prodotto con fit_marea.py) e ricostruisce la marea con
ttide.t_predic, come in 4-add-rm_tides.py (OMI).

Convenzione dei tempi:
  - fit_marea.py ha passato a t_tide stime = datetime(..., 0, 30, 0),
    cioe' il centro della prima media oraria (HH:30);
  - time_counter di EAS9-simu e' gia' a HH:30 (verificato con ncdump);
  => la sintesi si fa direttamente sui tempi di time_counter, senza lo
     spostamento di +30 min che serviva nello script OMI (dove l'asse
     era etichettato HH:00).
  I tempi vengono convertiti con ttide.time.date2num (convenzione
  ordinale di ttide), NON con matplotlib.dates.date2num, che da
  matplotlib 3.3 usa l'epoca 1970 ed e' incompatibile con ttide.
"""

import warnings
import numpy as np
import netCDF4 as nc
import ttide
from ttide.t_predic import t_predic
from ttide import time as ttide_time

warnings.filterwarnings("ignore", category=RuntimeWarning, module="ttide")


def get_tidal_frequencies(constituents):
    """
    Nomi e frequenze [cph] nel formato e nell'ORDINE attesi da
    t_predic (t_tide riordina per frequenza crescente: l'ordine di
    tidal_names NON e' quello di 'constituents'). Stesso trucco della
    serie fittizia usato in 4-add-rm_tides.py. Da chiamare una volta.
    """
    dummy = np.zeros(24 * 365)
    out = ttide.t_tide(dummy, dt=1, constitnames=constituents,
                       out_style=None, outfile=None)
    names = [n.decode().strip() if isinstance(n, (bytes, bytearray)) else n.strip()
             for n in out["nameu"]]
    return out["nameu"], out["fu"], names


def load_amppha_box(amppha_file, constituents, min_lat, max_lat, min_lon, max_lon):
    """
    Legge le sottomappe (y, x) di ampiezza [m] e fase [deg] per il box,
    con la stessa indicizzazione di sossheig (il file ha la griglia
    completa 380 x 1307, identica a quella della run).
    """
    amp, pha = {}, {}
    with nc.Dataset(amppha_file, "r") as ds:
        for c in constituents:
            amp[c] = np.array(ds.variables[f"{c}_Amp"][min_lat:max_lat, min_lon:max_lon], dtype=np.float64)
            pha[c] = np.array(ds.variables[f"{c}_Pha"][min_lat:max_lat, min_lon:max_lon], dtype=np.float64)
    return amp, pha


def times_to_ttide(times_dt, check_half_hour=True):
    """
    Converte una sequenza di datetime nella convenzione ordinale di
    ttide. Se check_half_hour, si ferma se i tempi non sono tutti a
    HH:30 (convenzione con cui sono state stimate le fasi).
    """
    times_dt = list(times_dt)
    if check_half_hour:
        bad = [t for t in times_dt if not (t.minute == 30 and t.second == 0)]
        if bad:
            raise ValueError(
                f"{len(bad)} tempi non a HH:30 (primo: {bad[0]}). Le fasi in "
                f"AMPPHA_FILE sono riferite a HH:30: controllare time_counter."
            )
    return np.array([ttide_time.date2num(t) for t in times_dt], dtype=np.float64)


def synth_tide_point(t_num, tidal_names_raw, tidal_freq, tidal_names,
                     amp_box, pha_box, j, i, lat):
    """
    Marea ricostruita [m] nel punto (j, i) del box, ai tempi t_num
    (gia' convertiti con times_to_ttide). tidecon costruito
    nell'ordine di tidal_names (ordine di t_tide), non di CONSTITUENTS.
    Ritorna None se ampiezze/fasi non sono valide nel punto.
    """
    tidecon = np.zeros((len(tidal_names), 4), dtype=np.float64)
    for k, name in enumerate(tidal_names):
        tidecon[k, 0] = amp_box[name][j, i]  # ampiezza [m]
        tidecon[k, 1] = 1.0                  # errore fittizio (>0; non usato con synth=0)
        tidecon[k, 2] = pha_box[name][j, i]  # fase [deg]
        tidecon[k, 3] = 0.0                  # errore fittizio
    if not np.all(np.isfinite(tidecon[:, [0, 2]])):
        return None
    tide = t_predic(t_num, tidal_names_raw, tidal_freq, tidecon,
                    lat=lat, ltype="nodal", synth=0)
    return np.asarray(tide, dtype=np.float64).ravel()
