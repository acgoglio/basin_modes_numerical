# Pipeline modi - run forzata (§3.3.2)

Analisi dei 10 modi di riferimento sulla run MedFS forzata (EAS9_minr_nt,
2020-2023): 4 run annuali, 16 trimestri, 7 eventi estremi.

## File

| File | Ruolo |
|---|---|
| `forced_windows_ini.py` | Configurazione (path, box, modi di riferimento, 27 finestre) |
| `build_file_list.py` | Lista esplicita dei file giornalieri per finestra |
| `f_forced_ampspt.py` | Spettro (bug fix filtro) + matching modi in banda |
| `run_forced_modes.py` | Calcolo per un box di griglia (CLI) |
| `merge_forced_modes.py` | Fusione dei 108 box + 99° percentile + istogrammi |
| `run_forced_BM.sh` | Lancio bsub (`run`/`status`/`merge`) |
| `peek_results.py` | Lettura rapida dei risultati (anche parziali) |
| `tab_3month_extrev_v2.py` | Plot finale (heatmap trimestri + eventi) |
| `rt_stats_tools.py` | Libreria utilizzata da `get_med_mask` |

## Setup (una volta sola)

Apri `forced_windows_ini.py` e verifica/imposta:
- `work_dir` - cartella di output sul cluster
- `ssh_varname` - nome variabile SSH nei netCDF (default `sossheig`)
- `forced_run_daily_template` - path dei file giornalieri della run forzata
- `mesh_mask`, `bathy_meter` - path griglia NEMO

## Lancio

Le 27 finestre sono divise in 3 categorie, lanciabili separatamente
(deciso dopo dei kill per wall-time sugli annuali): `event` (7,
leggeri), `trimester` (16), `annual` (4, i piu' pesanti - usano una
griglia piu' fine, 432 box invece di 108, per stare dentro le 24h).

```bash
./run_forced_BM.sh run event        # sottomette 108 job bsub, solo i 7 eventi
./run_forced_BM.sh status event      # controlla l'avanzamento (solo eventi)

./run_forced_BM.sh run trimester     # 108 job, solo i 16 trimestri
./run_forced_BM.sh status trimester

./run_forced_BM.sh run annual        # 432 job (griglia fine), solo i 4 annuali
./run_forced_BM.sh status annual

./run_forced_BM.sh merge             # fonde TUTTO quello che trova (tutte le categorie insieme)
python tab_3month_extrev_v2.py       # genera heatmap_trimesters_events.png
```

Puoi lanciare le tre categorie in un ordine qualsiasi, anche in
parallelo. `run_forced_BM.sh run` (senza categoria) lancia tutte e 27
le finestre in un solo job per box, come nella versione originale -
usalo solo se il wall-time non e' un problema per il tuo caso.

`status` e `merge` si possono rilanciare in qualunque momento, anche a
run in corso: leggono solo i file gia' scritti (scrittura atomica lato
`run_forced_modes.py`, quindi mai un file letto a meta').

## Controllare i risultati parziali durante il run

```bash
./run_forced_BM.sh merge     # snapshot con quello che c'e' finora
python peek_results.py       # tabella leggibile modi x trimestri/eventi/anni
```

Le celle non ancora calcolate mostrano `.` invece di un numero.

## Test su una sola finestra (prima del lancio completo)

In `forced_windows_ini.py`:
```python
TEST_SINGLE_WINDOW_KEY = "Storm Gloria"   # o altro evento/trimestre/anno
```
poi, a mano, su un box piccolo (non tramite `run_forced_BM.sh`):
```bash
python run_forced_modes.py 300 310 0 10 999
```
Ricorda di rimettere `TEST_SINGLE_WINDOW_KEY = None` prima del run vero.

## Se un job sembra bloccato

```bash
bjobs -a                     # RUN/PEND = ancora attivo, EXIT/DONE = terminato
grep -l "Error\|Traceback" ${WORK_DIR}/err_forced_*   # cerca errori reali
ls -la --time-style=full-iso ${WORK_DIR}/err_forced_<N>   # data dell'errore: vecchio o nuovo?
```
