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
# Schema a box (griglia NEMO), stesso schema di run_area_BM.sh /
# merge_amp_idx.py: 108 box (18 colonne x 6 righe), 4 esclusi
# (Baia di Biscaglia). Definiti QUI UNA SOLA VOLTA e importati sia da
# submit_forced_BM.py (lancio bsub) sia da merge_forced_modes.py
# (fusione box) - elimina la duplicazione tra file diversi.
# ---------------------------------------------------------------

BOX_X_EDGES = [300, 356, 412, 468, 524, 580, 636, 692, 748, 804,
               860, 916, 972, 1028, 1084, 1140, 1196, 1252, 1306]
BOX_Y_EDGES = [0, 63, 126, 189, 252, 315, 379]
BOX_N_ROWS = len(BOX_Y_EDGES) - 1  # 6
BOX_N_COLS = len(BOX_X_EDGES) - 1  # 18
BOX_EXCLUDE = [73, 74, 91, 92]     # Baia di Biscaglia

# ---------------------------------------------------------------
# Parametri di sottomissione job (bsub / LSF)
# ---------------------------------------------------------------

QUEUE = "s_long"
QUEUE_SHORT = "s_long"
QMEM = "400G"
QPRJ = "0723"

# ---------------------------------------------------------------
# Percorsi
# ---------------------------------------------------------------

# Run forzata: organizzata per cartella giornaliera YYYYMMDD/model/...
# NB: qui e' un TEMPLATE per costruire il path di un singolo giorno,
# non un pattern da glob - vedi build_file_list.py, che genera la
# lista esplicita di file invece di aprire tutto con '20*'.
forced_run_daily_template = (
    "/work/cmcc/med-dev/exp/EAS9_minr_nt/{yyyymmdd}/model/"
    "medfs-eas9_1h_{yyyymmdd}_2D_grid_T.nc"
)

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
work_dir = "/work/cmcc/ag15419/basin_modes_new/basin_modes_EAS9NT/" 

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

# 4 run annuali
for y in YEARS:
    WINDOWS.append(("annual", str(y), f"{y}0101", f"{y}1231"))

# 16 trimestri (calendario) x anno
for y in YEARS:
    for tname, (m0d0, m1d1) in TRIMESTERS.items():
        WINDOWS.append(("trimester", f"{tname}_{y}", f"{y}{m0d0}", f"{y}{m1d1}"))

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

assert len(WINDOWS) == 27, f"Attese 27 finestre, trovate {len(WINDOWS)}"

# Ordine fisso di modi/trimestri/eventi per l'output finale, IDENTICO
# a quello di tab_3month_extrev.py (cosi' l'array numpy prodotto qui
# si allinea direttamente, senza bisogno di rimappare indici).
TRIMESTER_ORDER = ["JFM", "AMJ", "JAS", "OND"]
EVENT_ORDER = [e[0] for e in EVENTS]
