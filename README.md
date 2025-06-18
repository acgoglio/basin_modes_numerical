# Prapare infiles
cdo selvar,MSL /data/inputs/METOCEAN/historical/model/atmos/ECMWF/IFS_0125/analysis/6h/netcdf/2014/01/20140101-ECMWF---AM0125-MEDATL-b20140102_an-fv10.00.nc /data/cmcc/ag15419/ECMWF_Med_modes/2014/12/20140101-ECMWF---AM0125-MEDATL-b20140102_an-fv10.00.nc_msl

cdo setrtoc,0.9e+05,2e+05,1.01e+05 /data/cmcc/ag15419/ECMWF_Med_modes/2014/12/20140101-ECMWF---AM0125-MEDATL-b20140102_an-fv10.00.nc_msl /data/cmcc/ag15419/ECMWF_Med_modes/2014/12/20140101-ECMWF---AM0125-MEDATL-b20140102_an-fv10.00.nc_cost

cp /data/inputs/METOCEAN/historical/model/atmos/ECMWF/IFS_0125/analysis/6h/netcdf/2014/01/20140101-ECMWF---AM0125-MEDATL-b20140102_an-fv10.00.nc /data/cmcc/ag15419/ECMWF_Med_modes/2014/12/

cdo replace /data/cmcc/ag15419/ECMWF_Med_modes/2014/12/20140101-ECMWF---AM0125-MEDATL-b20140102_an-fv10.00.nc /data/cmcc/ag15419/ECMWF_Med_modes/2014/12/20140101-ECMWF---AM0125-MEDATL-b20140102_an-fv10.00.nc_cost /data/cmcc/ag15419/ECMWF_Med_modes/2014/12/20140101-ECMWF---AM0125-MEDATL-b20140102_an-fv10.00.nc_ok

mv -f /data/cmcc/ag15419/ECMWF_Med_modes/2014/12/20140101-ECMWF---AM0125-MEDATL-b20140102_an-fv10.00.nc_ok /data/cmcc/ag15419/ECMWF_Med_modes/2014/12/20140101-ECMWF---AM0125-MEDATL-b20140102_an-fv10.00.nc 


# Compute Med Basin mode
#
# To run over an area:
# set and run run_basin_modes_cp.sh (it calls point_spt.py, which must be set too before running)
#
# If you just want to analyze a point with plots etc set and run point_spt_diag.py
# point_spt_diag2.py is the same of point_spt.py but includes diagnsotics
#
# NEW version:
f_point_area.py
run_basin_modes.py 

# SPython point_spt_diag_MSLP.py 200 80 prova
#
# LAST puntual VERSION
# while read LINE; do if [[ ${LINE:0:1} != '#' ]]; then echo $LINE ; MPython point_powspt_diag.py $LINE ; MPython point_ampspt_diag.py $LINE ; fi ; done < idx_10pt.coo

# LAST VERSION for the area:
## To run on the whole domain
# 1) Work on single box
LPython run_basin_modes_amp_idx.py 300  421  0 128 1
LPython run_basin_modes_amp_idx.py 421  541  0 128 2
LPython run_basin_modes_amp_idx.py 541  661  0 128 3
LPython run_basin_modes_amp_idx.py 661  781  0 128 4
LPython run_basin_modes_amp_idx.py 781  901  0 128 5
LPython run_basin_modes_amp_idx.py 901  1021 0 128 6
LPython run_basin_modes_amp_idx.py 1021 1141 0 128 7
LPython run_basin_modes_amp_idx.py 1141 1261 0 128 8
LPython run_basin_modes_amp_idx.py 1261 1306 0 128 9
LPython run_basin_modes_amp_idx.py 300  421  128 255 10
LPython run_basin_modes_amp_idx.py 421  541  128 255 11
LPython run_basin_modes_amp_idx.py 541  661  128 255 12
LPython run_basin_modes_amp_idx.py 661  781  128 255 13
LPython run_basin_modes_amp_idx.py 781  901  128 255 14
LPython run_basin_modes_amp_idx.py 901  1021 128 255 15
LPython run_basin_modes_amp_idx.py 1021 1141 128 255 16
LPython run_basin_modes_amp_idx.py 1141 1261 128 255 17
LPython run_basin_modes_amp_idx.py 1261 1306 128 255 18
LPython run_basin_modes_amp_idx.py 300  421  255 379 19
LPython run_basin_modes_amp_idx.py 421  541  255 379 20
LPython run_basin_modes_amp_idx.py 541  661  255 379 21
LPython run_basin_modes_amp_idx.py 661  781  255 379 22
LPython run_basin_modes_amp_idx.py 781  901  255 379 23
LPython run_basin_modes_amp_idx.py 901  1021 255 379 24
LPython run_basin_modes_amp_idx.py 1021 1141 255 379 25
LPython run_basin_modes_amp_idx.py 1141 1261 255 379 26
LPython run_basin_modes_amp_idx.py 1261 1306 255 379 27

LPython run_basin_modes_pow_idx.py 300  421  0 128 1
LPython run_basin_modes_pow_idx.py 421  541  0 128 2
LPython run_basin_modes_pow_idx.py 541  661  0 128 3
LPython run_basin_modes_pow_idx.py 661  781  0 128 4
LPython run_basin_modes_pow_idx.py 781  901  0 128 5
LPython run_basin_modes_pow_idx.py 901  1021 0 128 6
LPython run_basin_modes_pow_idx.py 1021 1141 0 128 7
LPython run_basin_modes_pow_idx.py 1141 1261 0 128 8
LPython run_basin_modes_pow_idx.py 1261 1306 0 128 9
LPython run_basin_modes_pow_idx.py 300  421  128 255 10
LPython run_basin_modes_pow_idx.py 421  541  128 255 11
LPython run_basin_modes_pow_idx.py 541  661  128 255 12
LPython run_basin_modes_pow_idx.py 661  781  128 255 13
LPython run_basin_modes_pow_idx.py 781  901  128 255 14
LPython run_basin_modes_pow_idx.py 901  1021 128 255 15
LPython run_basin_modes_pow_idx.py 1021 1141 128 255 16
LPython run_basin_modes_pow_idx.py 1141 1261 128 255 17
LPython run_basin_modes_pow_idx.py 1261 1306 128 255 18
LPython run_basin_modes_pow_idx.py 300  421  255 379 19
LPython run_basin_modes_pow_idx.py 421  541  255 379 20
LPython run_basin_modes_pow_idx.py 541  661  255 379 21
LPython run_basin_modes_pow_idx.py 661  781  255 379 22
LPython run_basin_modes_pow_idx.py 781  901  255 379 23
LPython run_basin_modes_pow_idx.py 901  1021 255 379 24
LPython run_basin_modes_pow_idx.py 1021 1141 255 379 25
LPython run_basin_modes_pow_idx.py 1141 1261 255 379 26
LPython run_basin_modes_pow_idx.py 1261 1306 255 379 27

# 2) Merge the boxes
SPython merge_amp_idx.py
SPython merge_pow_idx.py

# 3) Run the diagnostic
SPython mode_period_tab_amp.py      # Needed to extract the table of main modes to be analyzed
SPython mode_period_plot_ampval.py  # Amplitude values are analyzed 
SPython mode_period_plot_amp.py     # Map of the modes area 

SPython mode_period_tab_pow.py
SPython mode_period_plot_powval.py
SPython mode_period_plot_pow.p


# NBF case 
LPython run_basin_modes_amp_idx.py 781  901  0 64 5a
LPython run_basin_modes_amp_idx.py 781  901  64 128 5b
LPython run_basin_modes_amp_idx.py 901  1021 0 64 6a
LPython run_basin_modes_amp_idx.py 901  1021 64 128 6b
LPython run_basin_modes_amp_idx.py 1021 1141 0 64 7a
LPython run_basin_modes_amp_idx.py 1021 1141 64 128 7b
LPython run_basin_modes_amp_idx.py 1141 1261 0 64 8a
LPython run_basin_modes_amp_idx.py 1141 1261 64 128 8b
LPython run_basin_modes_amp_idx.py 421  541  128 191 11a
LPython run_basin_modes_amp_idx.py 421  541  191 255 11b
LPython run_basin_modes_amp_idx.py 541  661  128 191 12a
LPython run_basin_modes_amp_idx.py 541  661  191 255 12b
LPython run_basin_modes_amp_idx.py 661  781  128 191 13a
LPython run_basin_modes_amp_idx.py 661  781  191 255 13b
LPython run_basin_modes_amp_idx.py 781  901  128 191 14a
LPython run_basin_modes_amp_idx.py 781  901  191 255 14b
LPython run_basin_modes_amp_idx.py 901  1021 128 191 15a
LPython run_basin_modes_amp_idx.py 901  1021 191 255 15b
LPython run_basin_modes_amp_idx.py 1021 1141 128 191 16a
LPython run_basin_modes_amp_idx.py 1021 1141 191 255 16b

LPython run_basin_modes_pow_idx.py 781  901  0 64 5a
LPython run_basin_modes_pow_idx.py 781  901  64 128 5b
LPython run_basin_modes_pow_idx.py 901  1021 0 64 6a
LPython run_basin_modes_pow_idx.py 901  1021 64 128 6b
LPython run_basin_modes_pow_idx.py 1021 1141 0 64 7a
LPython run_basin_modes_pow_idx.py 1021 1141 64 128 7b
LPython run_basin_modes_pow_idx.py 1141 1261 0 64 8a
LPython run_basin_modes_pow_idx.py 1141 1261 64 128 8b
LPython run_basin_modes_pow_idx.py 421  541  128 191 11a
LPython run_basin_modes_pow_idx.py 421  541  191 255 11b
LPython run_basin_modes_pow_idx.py 541  661  128 191 12a
LPython run_basin_modes_pow_idx.py 541  661  191 255 12b
LPython run_basin_modes_pow_idx.py 661  781  128 191 13a
LPython run_basin_modes_pow_idx.py 661  781  191 255 13b
LPython run_basin_modes_pow_idx.py 781  901  128 191 14a
LPython run_basin_modes_pow_idx.py 781  901  191 255 14b
LPython run_basin_modes_pow_idx.py 901  1021 128 191 15a
LPython run_basin_modes_pow_idx.py 901  1021 191 255 15b
LPython run_basin_modes_pow_idx.py 1021 1141 128 191 16a
LPython run_basin_modes_pow_idx.py 1021 1141 191 255 16b



