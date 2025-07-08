#!/bin/bash
#
# by AC Goglio (CMCC)
# annachiara.goglio@cmcc.it
#
# Written: 04/07/2025
#
#set -u
set -e
#set -x
########################
# Read ini file
source $(pwd)/area_ini.py

echo "Work dir: ${work_dir}"
echo "Input files: ${file_template}"

# Submission parameters
QUEUE="s_long"
QUEUE_S="s_short"
QMEM="100G"
QPRJ="0723"

# 1) Define and Analyze 27 sub-domains
if [[ ${flag_compute_modes} == 1 ]]; then
   echo "I am going to compute the modes on 27 sub-domains.."   

   # Defines the sub-domains
   
   # +--+--+--+--+--+--+--+--+--+
   # |19|20|21|22|23|24|25|26|27|
   # +--+--+--+--+--+--+--+--+--+
   # |10|11|12|13|14|15|16|17|18|
   # +--+--+--+--+--+--+--+--+--+
   # |1 |2 |3 |4 |5 |6 |7 |8 |9 |
   # +--+--+--+--+--+--+--+--+--+
   
   min_lon_list=(300 421 541 661 781 901 1021 1141 1261 300 421 541 661 781 901 1021 1141 1261 300 421 541 661 781 901 1021 1141 1261)
   max_lon_list=(421 541 661 781 901 1021 1141 1261 1306 421 541 661 781 901 1021 1141 1261 1306 421 541 661 781 901 1021 1141 1261 1306)
   min_lat_list=(0   0   0   0   0   0    0    0    0    128 128 128 128 128 128 128 128 128 255 255 255 255 255 255 255 255 255)
   max_lat_list=(128 128 128 128 128 128 128 128 128 255 255 255 255 255 255 255 255 255 379 379 379 379 379 379 379 379 379)
   
   # Loop on the subdomains
   for i in "${!min_lon_list[@]}"; do
       min_lon=${min_lon_list[$i]}
       max_lon=${max_lon_list[$i]}
       min_lat=${min_lat_list[$i]}
       max_lat=${max_lat_list[$i]}
       box_idx=$((i + 1))
   
       echo "Running box $box_idx: lon $min_lon-$max_lon, lat $min_lat-$max_lat"
   
       bsub -n 1 -q ${QUEUE} -P ${QPRJ} -M ${QMEM} -o ${work_dir}/out_ampAspt_${box_idx} -e ${work_dir}/err_ampAspt_${box_idx} python run_basin_modes_amp_idx.py $min_lon $max_lon $min_lat $max_lat $box_idx
       bsub -n 1 -q ${QUEUE} -P ${QPRJ} -M ${QMEM} -o ${work_dir}/out_powAspt_${box_idx} -e ${work_dir}/err_powAspt_${box_idx} python run_basin_modes_pow_idx.py $min_lon $max_lon $min_lat $max_lat $box_idx
   done
fi

# 2) Merge the subdomains
if [[ ${flag_merge_modes} == 1 ]]; then
   echo "I am going to merge all the subdomains.."
   bsub -n 1 -q ${QUEUE_S} -P ${QPRJ} -M ${QMEM} -o ${work_dir}/out_mergeamp -e ${work_dir}/err_mergeamp python merge_amp_idx.py
   bsub -n 1 -q ${QUEUE_S} -P ${QPRJ} -M ${QMEM} -o ${work_dir}/out_mergepow -e ${work_dir}/err_mergepow python merge_pow_idx.py
fi

# 3) Run the diagnostic
if [[ ${flag_modes_analysis} == 1 ]]; then
   echo "I am going to analyze the modes.."
   bsub -n 1 -q ${QUEUE_S} -P ${QPRJ} -M ${QMEM} -o ${work_dir}/out_tabamp -e ${work_dir}/err_tabamp python mode_period_tab_amp.py
   bsub -n 1 -q ${QUEUE_S} -P ${QPRJ} -M ${QMEM} -o ${work_dir}/out_tabpow -e ${work_dir}/err_tabpow python mode_period_tab_pow.py
fi 
