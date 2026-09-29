# =====================================================================
# Configurazione per la pipeline di analisi dei modi sulla run FORZATA
# (MedFS + ECMWF, 2020-2023) - §3.3.2 della tesi.
#
# Sostituisce, per questo caso d'uso, area_ini.py: qui i modi NON
# vengono ri-scoperti alla cieca (niente peak-finding+grouping libero),
# ma estratti in banda attorno ai 10 periodi di riferimento gia' noti
# dall'esperimento free-oscillation (Tab. 3.1, colonna Numerico).
#
# ATTENZIONE (da confermare): i valori/percorsi segnati con "DA
# VERIFICARE" sono la mia migliore ricostruzione da quanto discusso in
# chat, ma non li ho mai visti scritti da te in un file - controllali
# prima del primo lancio vero.
# =====================================================================

import numpy as np

# ---------------------------------------------------------------
# Schema a box (griglia NEMO). Due schemi disponibili:
#   'coarse' - 108 box (18x6), usato per event/trimester
#   'fine'   - 432 box (36x12), usato per annual (le finestre piu'
#              pesanti: 365gg x 18 segmenti Welch per punto griglia -
#              box piu' piccoli = meno lavoro per job, per stare
#              dentro il wall-time di 24h dopo i kill TERM_RUNLIMIT)
#
# Il rettangolo della Baia di Biscaglia da escludere e' definito UNA
# VOLTA in coordinate griglia fisiche (verificato dai 4 box esclusi
# originali: x in [300,412], y in [252,379]) - da qui si ricalcola
# automaticamente quali box escludere in QUALUNQUE schema, invece di
# avere liste di indici scritte a mano che rischiano di disallinearsi
# quando cambia la suddivisione.
# ---------------------------------------------------------------

BISCAY_X_RANGE = (300, 412)
BISCAY_Y_RANGE = (252, 379)


def _subdivide_edges(edges, factor):
    """Inserisce (factor-1) punti equispaziati tra ogni coppia di edge
    consecutivi, per ottenere una griglia piu' fine mantenendo gli
    stessi confini fisici totali."""
    new_edges = [edges[0]]
    for a, b in zip(edges[:-1], edges[1:]):
        step = (b - a) / factor
        for k in range(1, factor + 1):
            new_edges.append(round(a + k * step))
    return new_edges


def _compute_exclude_boxes(x_edges, y_edges, n_cols, x_range, y_range):
    """Trova gli indici (1-based, row-major) dei box interamente
    contenuti nel rettangolo [x_range] x [y_range]."""
    excluded = []
    for r in range(len(y_edges) - 1):
        for c in range(len(x_edges) - 1):
            x0, x1 = x_edges[c], x_edges[c + 1]
            y0, y1 = y_edges[r], y_edges[r + 1]
            if x0 >= x_range[0] and x1 <= x_range[1] and y0 >= y_range[0] and y1 <= y_range[1]:
                excluded.append(r * n_cols + c + 1)
    return excluded


# Schema 'coarse' (108 box) - invariato rispetto a prima
_COARSE_X_EDGES = [300, 356, 412, 468, 524, 580, 636, 692, 748, 804,
                    860, 916, 972, 1028, 1084, 1140, 1196, 1252, 1306]
_COARSE_Y_EDGES = [0, 63, 126, 189, 252, 315, 379]
_COARSE_N_ROWS = len(_COARSE_Y_EDGES) - 1  # 6
_COARSE_N_COLS = len(_COARSE_X_EDGES) - 1  # 18
_COARSE_EXCLUDE = _compute_exclude_boxes(_COARSE_X_EDGES, _COARSE_Y_EDGES,
                                          _COARSE_N_COLS, BISCAY_X_RANGE, BISCAY_Y_RANGE)
assert _COARSE_EXCLUDE == [73, 74, 91, 92], (
    f"Lo schema coarse ricalcolato non riproduce i 4 box originali: {_COARSE_EXCLUDE}"
)

# Schema 'fine' (432 box = 4x108, ogni box coarse diviso in 2x2)
_FINE_X_EDGES = _subdivide_edges(_COARSE_X_EDGES, 2)
_FINE_Y_EDGES = _subdivide_edges(_COARSE_Y_EDGES, 2)
_FINE_N_ROWS = len(_FINE_Y_EDGES) - 1  # 12
_FINE_N_COLS = len(_FINE_X_EDGES) - 1  # 36
_FINE_EXCLUDE = _compute_exclude_boxes(_FINE_X_EDGES, _FINE_Y_EDGES,
                                        _FINE_N_COLS, BISCAY_X_RANGE, BISCAY_Y_RANGE)
assert _FINE_N_ROWS * _FINE_N_COLS == 432, f"Attesi 432 box fine, trovati {_FINE_N_ROWS*_FINE_N_COLS}"

BOX_SCHEMES = {
    "coarse": {"x_edges": _COARSE_X_EDGES, "y_edges": _COARSE_Y_EDGES,
               "n_rows": _COARSE_N_ROWS, "n_cols": _COARSE_N_COLS, "exclude": _COARSE_EXCLUDE},
    "fine":   {"x_edges": _FINE_X_EDGES, "y_edges": _FINE_Y_EDGES,
               "n_rows": _FINE_N_ROWS, "n_cols": _FINE_N_COLS, "exclude": _FINE_EXCLUDE},
}

# Quale schema usare per ciascuna categoria di finestra (deciso in chat)
CATEGORY_BOX_SCHEME = {"event": "coarse", "trimester": "fine", "annual": "fine"}


def get_box_scheme_for_category(category):
    scheme_name = CATEGORY_BOX_SCHEME[category]
    return BOX_SCHEMES[scheme_name]


# Retro-compatibilita' (codice esistente che si aspetta lo schema
# 'coarse' come default, es. per event/trimester)
BOX_X_EDGES = _COARSE_X_EDGES
BOX_Y_EDGES = _COARSE_Y_EDGES
BOX_N_ROWS = _COARSE_N_ROWS
BOX_N_COLS = _COARSE_N_COLS
BOX_EXCLUDE = _COARSE_EXCLUDE

# ---------------------------------------------------------------
# Parametri di sottomissione job (bsub / LSF)
# ---------------------------------------------------------------

QUEUE = "s_long"
QUEUE_SHORT = "s_long"
QUEUE_MEDIUM = "s_medium"
QMEM = "400G"
QPRJ = "0723"

# ---------------------------------------------------------------
# Percorsi
# ---------------------------------------------------------------

# Run forzata: organizzata per cartella giornaliera YYYYMMDD/model/...
# NB: qui e' un TEMPLATE per costruire il path di un singolo giorno,
# non un pattern da glob - vedi build_file_list.py, che genera la
# lista esplicita di file invece di aprire tutto con '20*'.
# ---------------------------------------------------------------
# Esperimento da analizzare (template e work_dir insieme, per non
# scrivere gli output della run con maree nella cartella della NT)
# ---------------------------------------------------------------
EXPERIMENT = "EAS9_simu"   # "EAS9_minr_nt" (senza maree) | "EAS9_simu" (con maree)

_EXPERIMENTS = {
    "EAS9_minr_nt": {
        "template": "/work/cmcc/med-dev/exp/EAS9_minr_nt/{yyyymmdd}/model/"
                    "medfs-eas9_1h_{yyyymmdd}_2D_grid_T.nc",
        "work_dir": "/work/cmcc/ag15419/basin_modes_new/basin_modes_EAS9NT/",
    },
    "EAS9_simu": {
        "template": "/work/cmcc/med-dev/exp/EAS9-simu/EXP00/{yyyymmdd}/model/"
                    "medfs-eas9_1h_{yyyymmdd}_2D_grid_T.nc",
        "work_dir": "/work/cmcc/ag15419/basin_modes_new/basin_modes_EAS9/",
    },
}
forced_run_daily_template = _EXPERIMENTS[EXPERIMENT]["template"]
work_dir = _EXPERIMENTS[EXPERIMENT]["work_dir"]

# ---------------------------------------------------------------
# Detiding (solo per la run con maree)
# ---------------------------------------------------------------
flag_detide = True
AMPPHA_FILE = ("/data/cmcc/ag15419/harmonic_analysis/quid_EAS9/area_EAS9_inserr_2022/"
               "amppha2D_0_sossheig_20220701_20221231_mod_medfs-eas9.nc")
TIDE_CONSTITUENTS = ["M2", "S2", "K1", "O1", "N2", "P1", "Q1", "K2"]

assert not (flag_detide and EXPERIMENT == "EAS9_minr_nt"), \
    "flag_detide=True su una run senza maree: si sottrarrebbe una marea che non c'e'."


# SSH variable name - DA VERIFICARE: assumo sia lo stesso 'sossheig'
# usato ovunque nella pipeline free-oscillation/BF-3. Se la run
# forzata usa un nome diverso va corretto qui.
ssh_varname = "sossheig"

# SSH time-series frequency [s]
dt = 3600

# Mesh mask / bathimetria - stessi file usati in area_ini.py
mesh_mask = "/work/cmcc/ag15419/VAA_paper/DATA0/mesh_mask.nc"
bathy_meter = "/work/cmcc/ag15419/VAA_paper/DATA0/bathy_meter.nc"

# Directory di output
work_dir = "/work/cmcc/ag15419/basin_modes_new/basin_modes_EAS9/" 

# ---------------------------------------------------------------
# Parametri spettro (stessi flag di area_ini.py, per compatibilita'
# con f_point_ampspt.py / f_point_powspt.py corretti)
# ---------------------------------------------------------------

flag_nfft = 1
N_fft = 512

flag_hanning = 1

flag_filter = "true"
th_filter = 40

# Impostazioni spettro DIVERSE per tipo di finestra (deciso in chat):
# - annual/trimester: Welch segmentato, 20gg/segmento (come BF-3)
# - event: nessuna segmentazione, un solo FFT su tutta la finestra
#   (massimizza la risoluzione spettrale per finestre corte, invece
#   di suddividerle ulteriormente in sotto-segmenti)
SPECTRUM_SETTINGS = {
    "annual":    {"flag_segmented_spectrum": True,  "segment_len_days": 20},
    "trimester": {"flag_segmented_spectrum": True,  "segment_len_days": 20},
    "event":     {"flag_segmented_spectrum": False, "segment_len_days": None},
}

# Qui NON servono piu' i filtri per soglia di ampiezza/energia sui
# picchi (quelli servivano alla discovery cieca): il matching e'
# in banda fissa sui 10 modi di riferimento. Li lascio a 0 per
# chiarezza - il modulo core non li usa piu'.
amplitude_threshold_ratio = 0.0
energy_threshold_ratio = 0.0

# ---------------------------------------------------------------
# Modi di riferimento: i primi 10 modi Mediterranei, periodo e
# tolleranza dalla Tab. 3.1 (colonna Numerico), tolleranza ricalcolata
# con extra_unc=0.10 (confermato in chat) sulla stessa formula di
# mode_period_tab_amp.py: tol = round_1sigfig(T^2 * delta_f * (1+extra_unc))
# con delta_f = 1/(20*24) cph (risoluzione BF-3, 20 giorni).
#
# NB: questa tolleranza NON e' ricalcolata sulla risoluzione spettrale
# reale di ciascuna finestra (10gg / trimestre / anno) - e' rimasto un
# punto aperto in chat (l'aumento a extra_unc=10% la allarga di poco).
# La tengo fissa per ora come deciso, ma il flag qui sotto permette di
# passare in futuro a una tolleranza per-finestra senza toccare il
# resto del codice.
# ---------------------------------------------------------------

REFERENCE_MODES = [
    # (label, T_ref [h], tolerance [h])
    ("Mediterranean Basin",        27.00, 2.00),
    ("Adriatic Sea",                21.00, 1.00),
    ("Aegean Sea/Gulf of Gabes",    13.80, 0.40),
    ("Gulf of Gabes/Aegean Sea",    11.90, 0.30),
    ("Adriatic Sea/Gulf of Gabes",  10.90, 0.30),
    ("Adriatic Sea/Gulf of Gabes",   9.50, 0.20),
    ("Tyrrhenian Sea/Alboran Sea",   8.50, 0.20),
    ("Adriatic Sea/Gulf of Gabes",   7.20, 0.10),
    ("Gulf of Gabes/Adriatic Sea",   6.90, 0.10),
    ("Alboran Sea",                  5.75, 0.08),
]
N_MODES = len(REFERENCE_MODES)  # 10

# Flag: se True, usa una tolleranza per-finestra ricalcolata sulla
# risoluzione spettrale reale della finestra invece della tolleranza
# fissa sopra. NON ANCORA IMPLEMENTATO - lasciato qui come promemoria
# esplicito del punto aperto, non silenziosamente ignorato.
flag_window_specific_tolerance = False

# Flag: se nessun picco cade in banda T_i+/-tol_i per un punto/modo:
#   True  -> fallback sul valore grezzo dello spettro al bin piu'
#            vicino a T_i (nessun NaN, ma valore meno affidabile)
#   False -> NaN (il punto non contribuisce a quel modo/finestra)
# Deciso in chat: si parte con False (NaN); 'used_fallback' resta
# comunque calcolato per QA, per capire quanti/dove sono i buchi.
flag_use_fallback = False

# ---------------------------------------------------------------
# Finestre da processare (27 totali)
# Ogni finestra: (kind, key, start_date 'YYYYMMDD', end_date 'YYYYMMDD')
# kind in {'annual','trimester','event'}
# key: usato per indicizzare l'output (es. per year+trimester, o year
#      per annuale, o nome evento)
# ---------------------------------------------------------------

YEARS = [2020, 2021, 2022, 2023]

TRIMESTERS = {
    "JFM": ("0101", "0331"),
    "AMJ": ("0401", "0630"),
    "JAS": ("0701", "0930"),
    "OND": ("1001", "1231"),
}

WINDOWS = []

# Ordine deciso in chat: dal piu' leggero (evento, 10-17gg, un solo
# FFT) al piu' pesante (annuale, 365gg, 18 segmenti Welch PER PUNTO
# GRIGLIA) - cosi', se un box viene ucciso per wall-time (RUNLIMIT),
# perde le finestre pesanti ma ha gia' completato quelle leggere,
# invece del contrario.

# 7 eventi (date confermate in chat)
EVENTS = [
    ("Storm Gloria",     "20200117", "20200126"),
    ("Medicane Ianos",   "20200915", "20200924"),
    ("Medicane Apollo",  "20211023", "20211102"),
    ("Storm Blas",       "20211102", "20211118"),  # 17 gg, non 10 - confermato cosi'
    ("Acqua Alta 2022",  "20221117", "20221127"),
    ("Cyclone Helios",   "20230206", "20230215"),
    ("Storm Daniel",     "20230904", "20230914"),  # corretto: Settembre, non Novembre
]
for name, d0, d1 in EVENTS:
    WINDOWS.append(("event", name, d0, d1))

# 16 trimestri (calendario) x anno
for y in YEARS:
    for tname, (m0d0, m1d1) in TRIMESTERS.items():
        WINDOWS.append(("trimester", f"{tname}_{y}", f"{y}{m0d0}", f"{y}{m1d1}"))

# 4 run annuali (le piu' pesanti, per ultime)
for y in YEARS:
    WINDOWS.append(("annual", str(y), f"{y}0101", f"{y}1231"))

assert len(WINDOWS) == 27, f"Attese 27 finestre, trovate {len(WINDOWS)}"

# ---------------------------------------------------------------
# Disaccoppiamento per categoria (deciso in chat, dopo i kill per
# wall-time sugli annuali): permette di lanciare separatamente
# 'event', 'trimester', 'annual' - cosi' un limite di 24h sugli
# annuali non blocca piu' anche eventi/trimestri nello stesso job.
#
# Uso: run_forced_modes.py accetta un 6' argomento opzionale, lista
# di categorie separate da virgola (es. "event" o "event,trimester").
# Default: tutte e tre insieme (comportamento originale).
# ---------------------------------------------------------------

ALL_CATEGORIES = ["event", "trimester", "annual"]


def get_windows_for_categories(categories):
    """categories: lista di stringhe tra ALL_CATEGORIES, o None/vuota
    per tutte e tre."""
    cats = categories if categories else ALL_CATEGORIES
    invalid = set(cats) - set(ALL_CATEGORIES)
    if invalid:
        raise ValueError(f"Categorie non valide: {invalid}. Valide: {ALL_CATEGORIES}")
    return [w for w in WINDOWS if w[0] in cats]

# ---------------------------------------------------------------
# Test rapido su UNA sola finestra (per validare I/O/formule prima
# del lancio completo sui 108 box x 27 finestre).
#
# Per attivarlo: metti qui la 'key' della finestra che vuoi testare
# (es. "Storm Gloria" per un evento, "JFM_2020" per un trimestre,
# "2020" per un anno), poi lancia run_forced_modes.py A MANO su un
# box piccolo (NON tramite run_forced_BM.sh, che lancia tutti i 108
# box). Esempio:
#   python run_forced_modes.py 300 310 0 10 999
#
# Lascia TEST_SINGLE_WINDOW_KEY = None per il run completo (27
# finestre, lanciato con run_forced_BM.sh run).
# ---------------------------------------------------------------
TEST_SINGLE_WINDOW_KEY = None  # es: "Storm Gloria"

if TEST_SINGLE_WINDOW_KEY is not None:
    _filtered = [w for w in WINDOWS if w[1] == TEST_SINGLE_WINDOW_KEY]
    if len(_filtered) != 1:
        raise ValueError(
            f"TEST_SINGLE_WINDOW_KEY='{TEST_SINGLE_WINDOW_KEY}' non trovata. "
            f"Chiavi valide: {[w[1] for w in WINDOWS]}"
        )
    WINDOWS = _filtered
    print(f"[TEST MODE] Solo la finestra '{TEST_SINGLE_WINDOW_KEY}' verra' processata.")

# Ordine fisso di modi/trimestri/eventi per l'output finale, IDENTICO
# a quello di tab_3month_extrev.py (cosi' l'array numpy prodotto qui
# si allinea direttamente, senza bisogno di rimappare indici).
TRIMESTER_ORDER = ["JFM", "AMJ", "JAS", "OND"]
EVENT_ORDER = [e[0] for e in EVENTS]
