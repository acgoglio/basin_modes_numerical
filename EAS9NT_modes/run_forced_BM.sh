#!/bin/bash
#
# Lancio della pipeline sulla run forzata (§3.3.2).
# Stessa sintassi bsub di run_area_BM.sh. I valori di work_dir/queue/box
# sono letti da forced_windows_ini.py (unica fonte di verita', nessuna
# duplicazione tra bash e python).
#
# Uso:
#   ./run_forced_BM.sh run [categoria]       # sottomette 108 job bsub
#   ./run_forced_BM.sh status [categoria]    # controlla l'avanzamento
#   ./run_forced_BM.sh merge                 # fonde i box, calcola 99th perc + istogrammi
#
# categoria: "event" | "trimester" | "annual" | "all" (default: "all")
#
# Disaccoppiamento deciso in chat dopo i kill per wall-time (24h) sulle
# finestre annuali: puoi lanciare separatamente
#   ./run_forced_BM.sh run event
#   ./run_forced_BM.sh run trimester
#   ./run_forced_BM.sh run annual
# cosi' un timeout sugli annuali (i piu' pesanti) non si porta via
# anche eventi/trimestri gia' completati nello stesso job.
#
set -e

STEP="${1:?Uso: $0 [run|status|merge] [event|trimester|annual|all]}"
CATEGORY="${2:-all}"

if [[ "${CATEGORY}" == "all" ]]; then
   CATEGORY_ARG=""       # run_forced_modes.py: nessun 6' argomento -> tutte le categorie
   CATEGORY_TAG="all"    # per i nomi dei file di log
else
   CATEGORY_ARG="${CATEGORY}"
   CATEGORY_TAG="${CATEGORY}"
fi

# Valori letti da forced_windows_ini.py (unica fonte di verita')
WORK_DIR=$(python3 -c "import forced_windows_ini as ini; print(ini.work_dir)")
QUEUE=$(python3 -c "import forced_windows_ini as ini; print(ini.QUEUE)")
QUEUE_M=$(python3 -c "import forced_windows_ini as ini; print(ini.QUEUE_MEDIUM)")
QUEUE_S=$(python3 -c "import forced_windows_ini as ini; print(ini.QUEUE_SHORT)")
QMEM=$(python3 -c "import forced_windows_ini as ini; print(ini.QMEM)")
QPRJ=$(python3 -c "import forced_windows_ini as ini; print(ini.QPRJ)")

mkdir -p "${WORK_DIR}"
echo "Work dir: ${WORK_DIR}"
echo "Categoria: ${CATEGORY_TAG}"

if [[ "${STEP}" == "run" ]]; then
   echo "Lancio i box per la categoria '${CATEGORY_TAG}'.."

   read -ra min_lon_list <<< "$(python3 -c "
import forced_windows_ini as ini
s = ini.get_box_scheme_for_category('${CATEGORY}') if '${CATEGORY}' != 'all' else ini.BOX_SCHEMES['coarse']
vals=[]
for r in range(s['n_rows']):
    for c in range(s['n_cols']):
        vals.append(s['x_edges'][c])
print(' '.join(str(v) for v in vals))
")"
   read -ra max_lon_list <<< "$(python3 -c "
import forced_windows_ini as ini
s = ini.get_box_scheme_for_category('${CATEGORY}') if '${CATEGORY}' != 'all' else ini.BOX_SCHEMES['coarse']
vals=[]
for r in range(s['n_rows']):
    for c in range(s['n_cols']):
        vals.append(s['x_edges'][c+1])
print(' '.join(str(v) for v in vals))
")"
   read -ra min_lat_list <<< "$(python3 -c "
import forced_windows_ini as ini
s = ini.get_box_scheme_for_category('${CATEGORY}') if '${CATEGORY}' != 'all' else ini.BOX_SCHEMES['coarse']
vals=[]
for r in range(s['n_rows']):
    for c in range(s['n_cols']):
        vals.append(s['y_edges'][r])
print(' '.join(str(v) for v in vals))
")"
   read -ra max_lat_list <<< "$(python3 -c "
import forced_windows_ini as ini
s = ini.get_box_scheme_for_category('${CATEGORY}') if '${CATEGORY}' != 'all' else ini.BOX_SCHEMES['coarse']
vals=[]
for r in range(s['n_rows']):
    for c in range(s['n_cols']):
        vals.append(s['y_edges'][r+1])
print(' '.join(str(v) for v in vals))
")"

   n_boxes=${#min_lon_list[@]}
   echo "Schema box per '${CATEGORY_TAG}': ${n_boxes} box"

   for i in "${!min_lon_list[@]}"; do
       min_lon=${min_lon_list[$i]}
       max_lon=${max_lon_list[$i]}
       min_lat=${min_lat_list[$i]}
       max_lat=${max_lat_list[$i]}
       box_idx=$((i + 1))

       bsub -n 1 -q ${QUEUE} -P ${QPRJ} -M ${QMEM} \
            -o ${WORK_DIR}/out_forced_${CATEGORY_TAG}_${box_idx} -e ${WORK_DIR}/err_forced_${CATEGORY_TAG}_${box_idx} \
            python run_forced_modes.py $min_lon $max_lon $min_lat $max_lat $box_idx ${CATEGORY_ARG}
   done

   echo "Sottomessi ${n_boxes} job (categoria: ${CATEGORY_TAG}). Controlla con: $0 status ${CATEGORY_TAG}"

elif [[ "${STEP}" == "status" ]]; then
   CATEGORY_TAG="${CATEGORY_TAG}" python3 - <<'PYEOF'
import os
import forced_windows_ini as ini

category = os.environ.get("CATEGORY_TAG", "all")
categories = None if category == "all" else [category]
windows = ini.get_windows_for_categories(categories)

# Il numero di box attesi dipende dalla categoria (108 per
# event/trimester, 432 per annual) - per 'all' non ha senso un unico
# conteggio univoco (categorie miste), quindi in quel caso mostriamo
# il conteggio per-finestra usando lo schema specifico di ciascuna.
def n_boxes_for(kind):
    s = ini.get_box_scheme_for_category(kind)
    return s["n_rows"] * s["n_cols"]

print(f"Categoria: {category} -> {len(windows)} finestre\n")

missing_by_window = {}
for kind, key, _, _ in windows:
    n_boxes = n_boxes_for(kind)
    n_present = sum(
        os.path.exists(os.path.join(ini.work_dir, f"forced_modes_{kind}_{key}_{b}.nc"))
        for b in range(1, n_boxes + 1)
    )
    if n_present < n_boxes:
        missing_by_window[(kind, key)] = (n_present, n_boxes)

if not missing_by_window:
    print("Tutti i box completati per tutte le finestre di questa categoria. Puoi lanciare il merge.")
else:
    print("Finestre non ancora complete (box presenti / attesi):")
    for (kind, key), (n_present, n_boxes) in missing_by_window.items():
        print(f"  {kind}/{key}: {n_present}/{n_boxes}")
PYEOF

elif [[ "${STEP}" == "merge" ]]; then
   echo "Lancio il merge (fai questo quando le categorie che ti interessano sono complete - controlla con '$0 status <categoria>').."
   echo "NB: il merge processa SEMPRE tutte le 27 finestre (quelle non ancora pronte vengono saltate con un avviso, non servono lanci separati)."
   echo "prova bsub -n 1 -q ${QUEUE_M} -P ${QPRJ} -M ${QMEM} "
   bsub -n 1 -q ${QUEUE_M} -P ${QPRJ} -M ${QMEM} \
        -o ${WORK_DIR}/out_forced_merge -e ${WORK_DIR}/err_forced_merge \
        python merge_forced_modes.py
   echo "Merge sottomesso."

else
   echo "Uso: $0 [run|status|merge] [event|trimester|annual|all]"
   exit 1
fi
