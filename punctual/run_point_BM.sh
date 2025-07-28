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
# Read ini file and copy it in the workdir
source $(pwd)/point_ini.py
echo "Ini file:$(pwd)/point_ini.py"
cp $(pwd)/point_ini.py ${work_dir}/

echo "Work dir: ${work_dir}"
echo "List of points to be analyzed in file: ${coo_file}"

# Submission parameters
QUEUE="s_medium"
QMEM="40G"
QPRJ="0723"

# Loop on points 
while read LINE; do 
 if [[ ${LINE:0:1} != '#' ]]; then 
    echo "Working on point: ${LINE}" 
    # Energy analysis
    bsub -n 1 -q ${QUEUE} -P ${QPRJ} -M ${QMEM} -o ${work_dir}/out_powspt -e ${work_dir}/err_powspt python point_powspt_diag.py $LINE 
    # Amplitude analysis
    bsub -n 1 -q ${QUEUE} -P ${QPRJ} -M ${QMEM} -o ${work_dir}/out_ampspt -e ${work_dir}/err_ampspt python point_ampspt_diag.py $LINE 
 fi 
done < ${coo_file}

# Plot the map with the selected points:
bsub -n 1 -q ${QUEUE} -P ${QPRJ} -M ${QMEM} -o ${work_dir}/out_pplot -e ${work_dir}/err_pplot python bathy_point_map.py
