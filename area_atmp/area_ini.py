# Initialization file for the med modes areal procedure

# Work directory
#work_dir='/work/cmcc/ag15419/basin_modes_num_wind_z/area/'
#work_dir='/work/cmcc/ag15419/basin_modes_num_wind_m/area/'
work_dir='/work/cmcc/ag15419/basin_modes_num_atmp_long_lw10/area/'

# Flags ( 0->NO; 1->YES):
flag_compute_modes=0
flag_merge_modes=0
flag_modes_analysis=0
flag_modes_plot=1

# Period to be analyzed (date format:YYYYMMDD)
start_date='20150105'
end_date='20150203'

# Input path/name template (use * character for the dates in the string)
# Relaxation after perturbation by Zonal wind 30 m/s
#file_template='/work/cmcc/ag15419/exp/fix_mfseas9_longrun_wind/EXP00/20*/model/medfs-eas9_1h_20*_2D_grid_T.nc'
# Relaxation after perturbation by Meridional wind 40 m/s
#file_template='/work/cmcc/ag15419/exp/fix_mfseas9_longrun_wind_m/EXP00_rel40m/20*/model/medfs-eas9_1h_20*_2D_grid_T.nc'
# Relaxation after perturbation by atmospheric pressure (Pref + 100 hPa)
file_template='/work/cmcc/ag15419/exp/EAS9BT_med-modes_atmp/EXP00/20*/model/medfs-eas9_1h_20*_2D_grid_T.nc'

# SSH time-serie frequency in seconds
dt=3600

# Mesh mask
mesh_mask='/work/cmcc/ag15419/VAA_paper/DATA0/mesh_mask.nc'
# Bathymetry
bathy_meter='/work/cmcc/ag15419/VAA_paper/DATA0/bathy_meter.nc'

######################################################

# Flag and threshold [h] for filtering the spectrum the threshold is also used as plot minimum
flag_filter='true'
th_filter=40

# Filter out modes with low amplitude, e.g. 0.10 means that all the modes with amplitude<10% of the total amplitude are rm (to avoid the filtering set amplitude_threshold_ratio = 0)
amplitude_threshold_ratio=0.0001
# Filter out modes with low energy, e.g. 0.10 means that all the modes with energy<10% of the total energy in the analized point are rm (to avoid the filtering set energy_threshold_ratio = 0)
energy_threshold_ratio=0.0001

# Number of modes to analyze (n_modes = 'auto' if you want to analyze all the modes)
n_modes='auto'

# Flag: use segmented (averaged) spectrum or full time series
flag_segmented_spectrum=True  # False to disable and use full spectrum
segment_len_days=10  # length of each segment (in days) if segmented spectrum is used

# To order by period instead of by amplitude set flag_T_order = 1
flag_T_order=1

# Template for output file
infile_amppha='/data/cmcc/ag15419/basin_modes/basin_modes_ini.nc'

# Uncertainty of the modes periods (1 = variable uncertainty, 0 = fixed uncertainty)
flag_var_unc=0         
# In case of fixed uncertainty = 0 fix the value
fixed_uncertainty=0.0001
# In case of variable uncertainty add a % to the uncertainty due to the spectral resolution (e.g. 0.1 means (1+0.1)*theoretical_uncertainty )
extra_unc=0.10
