#!/bin/bash
#
# by AC Goglio (CMCC)
# annachiara.goglio@cmcc.it
#
# Written: 4 Jul 2025
# Last modified: 18 Aug 2025 
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

# 0) Cp ini file in workdir
cp -vf $(pwd)/area_ini.py ${work_dir}/

# 1) Define and Analyze 27 sub-domains
if [[ ${flag_compute_modes} == 1 ]]; then
   echo "I am going to compute the modes on 108 sub-domains.."   

   # Defines the sub-domains
   
   #+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----+
   #|091 |092 |093 |094 |095 |096 |097 |098 |099 |100 |101 |102 |103 |104 |105 |106 |107 |108 |
   #+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----+
   #|073 |074 |075 |076 |077 |078 |079 |080 |081 |082 |083 |084 |085 |086 |087 |088 |089 |090 |
   #+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----+
   #|055 |056 |057 |058 |059 |060 |061 |062 |063 |064 |065 |066 |067 |068 |069 |070 |071 |072 |
   #+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----+
   #|037 |038 |039 |040 |041 |042 |043 |044 |045 |046 |047 |048 |049 |050 |051 |052 |053 |054 |
   #+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----+
   #|019 |020 |021 |022 |023 |024 |025 |026 |027 |028 |029 |030 |031 |032 |033 |034 |035 |036 |
   #+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----+
   #|001 |002 |003 |004 |005 |006 |007 |008 |009 |010 |011 |012 |013 |014 |015 |016 |017 |018 |
   #+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----+
   
   min_lon_list=(300 356 412 468 524 580 636 692 748 804 860 916 972 1028 1084 1140 1196 1252
                 300 356 412 468 524 580 636 692 748 804 860 916 972 1028 1084 1140 1196 1252
                 300 356 412 468 524 580 636 692 748 804 860 916 972 1028 1084 1140 1196 1252
                 300 356 412 468 524 580 636 692 748 804 860 916 972 1028 1084 1140 1196 1252
                 300 356 412 468 524 580 636 692 748 804 860 916 972 1028 1084 1140 1196 1252
                 300 356 412 468 524 580 636 692 748 804 860 916 972 1028 1084 1140 1196 1252)
   
   max_lon_list=(356 412 468 524 580 636 692 748 804 860 916 972 1028 1084 1140 1196 1252 1306
                 356 412 468 524 580 636 692 748 804 860 916 972 1028 1084 1140 1196 1252 1306
                 356 412 468 524 580 636 692 748 804 860 916 972 1028 1084 1140 1196 1252 1306
                 356 412 468 524 580 636 692 748 804 860 916 972 1028 1084 1140 1196 1252 1306
                 356 412 468 524 580 636 692 748 804 860 916 972 1028 1084 1140 1196 1252 1306
                 356 412 468 524 580 636 692 748 804 860 916 972 1028 1084 1140 1196 1252 1306)
   
   min_lat_list=(0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
                 63 63 63 63 63 63 63 63 63 63 63 63 63 63 63 63 63 63
                 126 126 126 126 126 126 126 126 126 126 126 126 126 126 126 126 126 126
                 189 189 189 189 189 189 189 189 189 189 189 189 189 189 189 189 189 189
                 252 252 252 252 252 252 252 252 252 252 252 252 252 252 252 252 252 252
                 315 315 315 315 315 315 315 315 315 315 315 315 315 315 315 315 315 315)
   
   max_lat_list=(63 63 63 63 63 63 63 63 63 63 63 63 63 63 63 63 63 63
                 126 126 126 126 126 126 126 126 126 126 126 126 126 126 126 126 126 126
                 189 189 189 189 189 189 189 189 189 189 189 189 189 189 189 189 189 189
                 252 252 252 252 252 252 252 252 252 252 252 252 252 252 252 252 252 252
                 315 315 315 315 315 315 315 315 315 315 315 315 315 315 315 315 315 315
                 379 379 379 379 379 379 379 379 379 379 379 379 379 379 379 379 379 379)

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
   bsub -n 1 -q ${QUEUE} -P ${QPRJ} -M ${QMEM} -o ${work_dir}/out_mergeamp -e ${work_dir}/err_mergeamp python merge_amp_idx.py
   bsub -n 1 -q ${QUEUE} -P ${QPRJ} -M ${QMEM} -o ${work_dir}/out_mergepow -e ${work_dir}/err_mergepow python merge_pow_idx.py
fi

# 3) Compute the modes
if [[ ${flag_modes_analysis} == 1 ]]; then
   echo "I am going to extract the modes.."
   bsub -n 1 -q ${QUEUE_S} -P ${QPRJ} -M ${QMEM} -o ${work_dir}/out_tabamp -e ${work_dir}/err_tabamp python mode_period_tab_amp.py
   bsub -n 1 -q ${QUEUE_S} -P ${QPRJ} -M ${QMEM} -o ${work_dir}/out_tabpow -e ${work_dir}/err_tabpow python mode_period_tab_pow.py
fi 

# 4) Plot the modes amplitudes
if [[ ${flag_modes_plot} == 1 ]]; then
   echo "I am going to plot the modes.."
bsub -n 1 -q ${QUEUE_S} -P ${QPRJ} -M ${QMEM} -o ${work_dir}/out_plotamp -e ${work_dir}/err_plotamp python mode_period_plot_ampval.py
bsub -n 1 -q ${QUEUE_S} -P ${QPRJ} -M ${QMEM} -o ${work_dir}/out_plotpow -e ${work_dir}/err_plotpow python mode_period_plot_powval.py

fi
