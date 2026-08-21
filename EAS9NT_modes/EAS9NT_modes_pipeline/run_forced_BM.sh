#!/bin/bash
#
# Lancio della pipeline sulla run forzata (§3.3.2).
# Stessa sintassi bsub di run_area_BM.sh. I valori di work_dir/queue/box
# sono letti da forced_windows_ini.py (unica fonte di verita', nessuna
# duplicazione tra bash e python).
#
# Uso:
#   ./run_forced_BM.sh run       # sottomette i 108 box (ciascuno cicla
#                                 # internamente sulle 27 finestre)
#   ./run_forced_BM.sh status    # controlla quanti box/finestre sono completi
#   ./run_forced_BM.sh merge     # fonde i box, calcola 99th perc + istogrammi
#
set -e

STEP="${1:?Uso: $0 [run|status|merge]}"

# Valori letti da forced_windows_ini.py (unica fonte di verita')
WORK_DIR=$(python3 -c "import forced_windows_ini as ini; print(ini.work_dir)")
QUEUE=$(python3 -c "import forced_windows_ini as ini; print(ini.QUEUE)")
QUEUE_S=$(python3 -c "import forced_windows_ini as ini; print(ini.QUEUE_SHORT)")
QMEM=$(python3 -c "import forced_windows_ini as ini; print(ini.QMEM)")
QPRJ=$(python3 -c "import forced_windows_ini as ini; print(ini.QPRJ)")

mkdir -p "${WORK_DIR}"
echo "Work dir: ${WORK_DIR}"

if [[ "${STEP}" == "run" ]]; then
   echo "Lancio 108 box, ciascuno con le 27 finestre interne.."

   # Box edges letti da forced_windows_ini.py, non riscritti qui
   read -ra min_lon_list <<< "$(python3 -c "
import forced_windows_ini as ini
vals=[]
for r in range(ini.BOX_N_ROWS):
    for c in range(ini.BOX_N_COLS):
        vals.append(ini.BOX_X_EDGES[c])
print(' '.join(str(v) for v in vals))
")"
   read -ra max_lon_list <<< "$(python3 -c "
import forced_windows_ini as ini
vals=[]
for r in range(ini.BOX_N_ROWS):
    for c in range(ini.BOX_N_COLS):
        vals.append(ini.BOX_X_EDGES[c+1])
print(' '.join(str(v) for v in vals))
")"
   read -ra min_lat_list <<< "$(python3 -c "
import forced_windows_ini as ini
vals=[]
for r in range(ini.BOX_N_ROWS):
    for c in range(ini.BOX_N_COLS):
        vals.append(ini.BOX_Y_EDGES[r])
print(' '.join(str(v) for v in vals))
")"
   read -ra max_lat_list <<< "$(python3 -c "
import forced_windows_ini as ini
vals=[]
for r in range(ini.BOX_N_ROWS):
    for c in range(ini.BOX_N_COLS):
        vals.append(ini.BOX_Y_EDGES[r+1])
print(' '.join(str(v) for v in vals))
")"

   for i in "${!min_lon_list[@]}"; do
       min_lon=${min_lon_list[$i]}
       max_lon=${max_lon_list[$i]}
       min_lat=${min_lat_list[$i]}
       max_lat=${max_lat_list[$i]}
       box_idx=$((i + 1))

       bsub -n 1 -q ${QUEUE} -P ${QPRJ} -M ${QMEM} \
            -o ${WORK_DIR}/out_forced_${box_idx} -e ${WORK_DIR}/err_forced_${box_idx} \
            python run_forced_modes.py $min_lon $max_lon $min_lat $max_lat $box_idx
   done

   echo "Sottomessi 108 job. Controlla con: $0 status"

elif [[ "${STEP}" == "status" ]]; then
   # Riusa la logica di conteggio gia' scritta in submit_forced_BM.py,
   # invece di riscriverla in bash
   python3 submit_forced_BM.py status

elif [[ "${STEP}" == "merge" ]]; then
   echo "Lancio il merge (fai questo SOLO dopo che tutti i 108 box sono completati - controlla con '$0 status').."
   bsub -n 1 -q ${QUEUE_S} -P ${QPRJ} -M ${QMEM} \
        -o ${WORK_DIR}/out_forced_merge -e ${WORK_DIR}/err_forced_merge \
        python merge_forced_modes.py
   echo "Merge sottomesso."

else
   echo "Uso: $0 [run|status|merge]"
   exit 1
fi
