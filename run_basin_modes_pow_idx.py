import glob
import sys
import xarray as xr
import netCDF4 as nc
import numpy as np
from scipy.signal import periodogram
import heapq
from functools import partial
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from scipy.signal import detrend
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks
import shutil
import f_point_powspt
mpl.use('Agg')

########

# Lon/lat box indexes
min_lon = int(sys.argv[1])
max_lon = int(sys.argv[2])
min_lat = int(sys.argv[3])
max_lat = int(sys.argv[4])

box_idx = str(sys.argv[5])

# Workdir, otufile template and otfile name
work_dir='/work/cmcc/ag15419/basin_modes_20not/'
infile_amppha='/work/cmcc/ag15419/basin_modes/basin_modes_ini.nc'
outfile=work_dir+'basin_modes_pow_'+box_idx+'.nc'

# Infiles
start_date = "20150121"
end_date = "20150221"
all_files = sorted(glob.glob("/work/cmcc/ag15419/exp/fix_mfseas9_longrun_hmslp_2NT_AB_2/EXP00/20*/model/medfs-eas9_1h_20*_2D_grid_T.nc"))
mesh_mask = "/work/cmcc/ag15419/VAA_paper/DATA0/mesh_mask.nc"

# Model time step in seconds
dt = 3600

###################
# Build the outfile:
shutil.copy(infile_amppha,outfile)

# Read infiles
# Select the period
infile = []
for f in all_files:
    parts = f.split("/")
    file_date = parts[7] # 7 6  
    if start_date <= file_date <= end_date: 
            infile.append(f)

# Read lat and lon in the first file
nav_lat = None
nav_lon = None
first=0
for nc2open in infile:
  if first==0:
    model = nc.Dataset(nc2open, 'r')
    if nav_lat is None:  
        nav_lat = model.variables['nav_lat'][:]
        nav_lon = model.variables['nav_lon'][:]
        first=1
    model.close()

# Read ssh xarray
ds = xr.open_mfdataset(infile, combine='by_coords', parallel=True)
ssh_ts_all = ds['sossheig']

# Read land/sea mask 
mesh_nemo = nc.Dataset(mesh_mask, 'r')
mesh = mesh_nemo.variables['tmask'][0,0,:,:]

########################
# Compute and write values in the netCDF file
modes_outfile = nc.Dataset(outfile, 'a')

# Call the function for each point in the Med
for lon_idx in range (min_lon,max_lon): #(300,len(nav_lon)):
    for lat_idx in range (min_lat,max_lat): # (0,len(nav_lat)):

        # If is sea-point:
        if mesh[lat_idx, lon_idx]==1:

           # Extract the time-series and Call the function
           ssh_ts_point = ssh_ts_all[:, lat_idx, lon_idx].values
           pow_peak_periods_main, pow_peak_amplitudes_main = f_point_powspt.pow_main_modes(lat_idx, lon_idx, ssh_ts_point, dt)

           for i in range(8):
               try:
                  modes_outfile.variables[f'm{i}_Amp'][lat_idx, lon_idx] = pow_peak_amplitudes_main[i]
                  modes_outfile.variables[f'm{i}_T'][lat_idx, lon_idx] = pow_peak_periods_main[i]
               except:
                  modes_outfile.variables[f'm{i}_Amp'][lat_idx, lon_idx] = np.nan
                  modes_outfile.variables[f'm{i}_T'][lat_idx, lon_idx] = np.nan

        # If land point
        else:
           for i in range(8):
                  modes_outfile.variables[f'm{i}_Amp'][lat_idx, lon_idx] = np.nan
                  modes_outfile.variables[f'm{i}_T'][lat_idx, lon_idx] = np.nan
modes_outfile.close()
