# Initialization file for the med modes punctual procedure

# Work directory
work_dir="/work/cmcc/ag15419/basin_modes_num_atmp_long_lw10/point/"

# Period to be analyzed (date format:YYYYMMDD)
start_date="20150104"
end_date="20150203"

# Input path/name template (use * character for the dates in the string)
# Relaxation after perturbation by Zonal wind 30 m/s
#file_template="/work/cmcc/ag15419/exp/EAS9BT_med-modes_wind_z/EXP00/20*/model/medfs-eas9_1h_20*_2D_grid_T.nc"
# Relaxation after perturbation by Meridional wind 40 m/s
#file_template="/work/cmcc/ag15419/exp/EAS9BT_med-modes_wind_m/EXP00_rel40m/20*/model/medfs-eas9_1h_20*_2D_grid_T.nc"
# Relaxation after perturbation by atmospheric pressure (Pref + 100 hPa)
file_template="/work/cmcc/ag15419/exp/EAS9BT_med-modes_atmp/EXP00/20*/model/medfs-eas9_1h_20*_2D_grid_T.nc"

# Exp tag
tag="atmp_"

# SSH time-serie frequency in seconds (e.g. for hourly ts 3600)
dt=3600 

# Bathymetry
bathy_meter="/work/cmcc/ag15419/VAA_paper/DATA0/bathy_meter.nc"

# Mesh_mask
mesh_mask="/work/cmcc/ag15419/VAA_paper/DATA0/mesh_mask.nc"

#######################################

# Number of modes to analyze (n_modes = 'auto' if you want to analyze all the modes)
n_modes='auto'

# Flag and threshold [h] for filtering the spectrum. The threshold is also used as plot minimum
flag_filter='true'
th_filter=40

# Flag for Gaussian smoothing of the spectrum: true, false or plot (to use the original spt but add the plot of the smoothed spt)
flag_smooth='false'
sigma=15

# To order by period instead of by amplitude set flag_T_order = 1
flag_T_order=1

# Filter out modes with low amplitude, e.g. 0.10 means that all the modes with amplitude<10% of the total amplitude are rm (to avoid the filtering set amplitude_threshold_ratio = 0)
amplitude_threshold_ratio=0.002
# Filter out modes with low energy, e.g. 0.10 means that all the modes with energy<10% of the total energy are rm (to avoid the filtering set energy_threshold_ratio = 0)
energy_threshold_ratio=0.002

# Flag: use segmented (averaged) spectrum or full time series
flag_segmented_spectrum=True  # False to disable and use full spectrum
segment_len_days=10  # length of each segment (in days) if segmented spectrum is used

# Coordinate file (list of points to be analyzed)
coo_file='idx_pt.coo'
