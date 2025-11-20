import glob
import sys
import re
import os
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
import f_point_ampspt
from area_ini import *
mpl.use('Agg')

########

# Lon/lat box indexes
min_lon = int(sys.argv[1])
max_lon = int(sys.argv[2])
min_lat = int(sys.argv[3])
max_lat = int(sys.argv[4])

box_idx = str(sys.argv[5])

###################
# Input files
all_files = sorted(glob.glob(file_template))

# Outfiles
outfile=work_dir+'basin_modes_amp_'+box_idx+'.nc'

# Build the outfile:
shutil.copy(infile_amppha,outfile)

# Read infiles
# Select the period
# Cylc archive structure:
#infile = []
#for f in all_files:
#    print ('file',f)
#    parts = f.split("/")
#    file_date = parts[7] # 7 6
#    if start_date <= file_date <= end_date:
#            infile.append(f)
# General archive structure:
print("Num of in files:", len(all_files))
if len(all_files) == 0:
    print("WARNING: template file not found ", file_template)
date_re = re.compile(r'(\d{8})')
infile = []
no_date_files = []
for f in all_files:
    basename = os.path.basename(f)
    m = date_re.search(basename)
    if m:
        file_date = m.group(1)
        if start_date <= file_date <= end_date:
            infile.append(f)
    else:
        no_date_files.append(f)
# Debug
print("Selected files:", len(infile))
if len(no_date_files) > 0:
    print("This files do not have a date in the name..:")
    for ff in no_date_files[:5]:
        print("  ", ff)
if len(infile) == 0:
    print("WARNING: NO files range", start_date, "-", end_date)

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

# Fields in the ini file
initial_modes = 8

#Call the function for each point in the Med
for lon_idx in range(min_lon, max_lon):
    for lat_idx in range(min_lat, max_lat):
        # If is sea-point:
        if mesh[lat_idx, lon_idx] == 1:
            # Extract the time-series and Call the function
            ssh_ts_point = ssh_ts_all[:, lat_idx, lon_idx].values
            amp_peak_periods_main, amp_peak_amplitudes_main = f_point_ampspt.amp_main_modes(
                lat_idx, lon_idx, ssh_ts_point, dt
            )

            nmodes_here = len(amp_peak_periods_main)

            for i in range(nmodes_here):
                # If num field > 8 built it
                amp_varname = f'm{i}_Amp'
                per_varname = f'm{i}_T'

                if amp_varname not in modes_outfile.variables:
                    modes_outfile.createVariable(amp_varname, 'f4', ('y', 'x'), fill_value=np.nan)
                if per_varname not in modes_outfile.variables:
                    modes_outfile.createVariable(per_varname, 'f4', ('y', 'x'), fill_value=np.nan)

                # Write the value
                modes_outfile.variables[amp_varname][lat_idx, lon_idx] = amp_peak_amplitudes_main[i]
                modes_outfile.variables[per_varname][lat_idx, lon_idx] = amp_peak_periods_main[i]

            # Fill extra fields with nan
            for i in range(nmodes_here, initial_modes):
                amp_varname = f'm{i}_Amp'
                per_varname = f'm{i}_T'

                if amp_varname in modes_outfile.variables:
                    modes_outfile.variables[amp_varname][lat_idx, lon_idx] = np.nan
                if per_varname in modes_outfile.variables:
                    modes_outfile.variables[per_varname][lat_idx, lon_idx] = np.nan

        # If land point
        else:
           for i in range(initial_modes):
               amp_varname = f'm{i}_Amp'
               per_varname = f'm{i}_T'
               if amp_varname in modes_outfile.variables:
                   modes_outfile.variables[amp_varname][lat_idx, lon_idx] = np.nan
               if per_varname in modes_outfile.variables:
                   modes_outfile.variables[per_varname][lat_idx, lon_idx] = np.nan

modes_outfile.close()
