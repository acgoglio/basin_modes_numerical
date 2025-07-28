from netCDF4 import Dataset
import sys
import netCDF4 as ncdf
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt
from copy import copy, deepcopy
import glob
import dateutil.parser
import rt_stats_tools
import rt_plotbox
import warnings
import xarray as xr
import netCDF4 as nc
from scipy.signal import periodogram
import heapq
from functools import partial
from matplotlib.colors import ListedColormap

warnings.filterwarnings("ignore")
mpl.use('Agg')

############################
exp_name      = "Med Modes excited by meridional wind forcing"
NEMO_GRID     = "/work/cmcc/ag15419/VAA_paper/DATA0/mesh_mask.nc"
var           = "sossheig"
infiles       = sorted(glob.glob("/work/cmcc/ag15419/exp/EAS9BT_med-modes_wind_z/EXP00/20*/model/medfs-eas9_1h_20*_2D_grid_T.nc"))
outfile_temp  = "/work/cmcc/ag15419/basin_modes_num_wind_m/ssh/ssh_"
minV          = "-0.5"
maxV          = "0.5"
thV           = "-1000"
cmap          = "seismic"
inidate       = "20150105"
enddate       = "20150107"
var_name      = "sossheig"
unit          = "m"

# Read the grid
nc_nemogrid = ncdf.Dataset(NEMO_GRID, 'r')
try:
    # For EAS7 and earlier
    nemo_lat = nc_nemogrid.variables['nav_lat'][:]
    nemo_lon = nc_nemogrid.variables['nav_lon'][:]
    nemo_mask = (nc_nemogrid.variables['tmask'][0, 0, :, :] == 0)
    nemo_mbathy = np.squeeze(nc_nemogrid.variables['mbathy'][0, :, :])
except:
    nemo_lat = nc_nemogrid.variables['y'][:]
    nemo_lon = nc_nemogrid.variables['x'][:]
    nemo_mask = (nc_nemogrid.variables['tmask'][0, 0, :, :] == 0)
    nemo_mbathy = np.squeeze(nc_nemogrid.variables['mbathy'][0, :, :])
nc_nemogrid.close()

# Select time range
infile = []
print ('Infiles:',infiles)
for f in infiles:
    print('file', f)
    parts = f.split("/")
    file_date = parts[7]
    if inidate <= file_date <= enddate:
        infile.append(f)

# Read the NetCDF data
first = 0
for nc2open in infile:
    print('Processing:', nc2open)
    model = nc.Dataset(nc2open, 'r')
    ssh_ts = np.array(model.variables['sossheig'])
    print('ssh_ts shape', ssh_ts.shape)
    if first == 0:
        ssh_ts_all = ssh_ts
        first = 1
    else:
        ssh_ts_all = np.concatenate((ssh_ts_all, ssh_ts))
        print('ssh_ts_all shape', ssh_ts_all.shape)
    model.close()

hours = len(ssh_ts_all[:, 0, 0])
print('Num of hours', hours)

# Time loop
for idx_time in range(0, hours):

    bat1 = np.array(ssh_ts_all[idx_time, :, :])

    nemo_mbathy = nemo_mbathy - 1
    nemo_mbathy = np.where(nemo_mbathy < 0, 0, nemo_mbathy)

    # Surface layer
    wrk = np.squeeze(bat1)
    dims = np.shape(bat1)

    # Mask land
    mask_Med = ~rt_stats_tools.get_med_mask(nemo_lon, nemo_lat, ~nemo_mask)
    wrk = np.ma.array(wrk, mask=mask_Med)

    # Prepare color palette
    try:
        step = (float(maxV) - float(minV)) / 100
        pa = np.arange(float(minV), float(maxV), step)
    except:
        print("WARNING: something wrong with this field.. I am fixing -100,100 as min max palette values")
        step = 25
        pa = np.arange(round(float(-100), 2), round(float(100), 2), step)

    # Create figure with 2 vertically aligned subplots
    fig, axes = plt.subplots(2, 1, figsize=(14, 14))

    # First subplot: full Mediterranean map
    ax1 = axes[0]
    m = Basemap(projection='cyl', llcrnrlat=30.2, urcrnrlat=46.,
                llcrnrlon=-6.0, urcrnrlon=36.3, resolution='h', ax=ax1)
    m.drawparallels(np.arange(30., 46., 5.), labels=[1, 0, 0, 0], fontsize=15)
    m.drawmeridians(np.arange(-18., 36.3, 5.), labels=[0, 0, 0, 1], fontsize=15)

    CS = m.contourf(nemo_lon, nemo_lat, np.squeeze(wrk), pa, cmap=cmap, extend='both')
    if thV != 'nan':
        CSth = m.contour(nemo_lon, nemo_lat, np.squeeze(wrk), [thV], color='green')
    m.drawcoastlines(linewidth=1.0, linestyle='solid', color='k', ax=ax1, zorder=None)
    m.drawmapboundary(fill_color='white')
    m.fillcontinents()
    cbar = plt.colorbar(CS, orientation='horizontal', aspect=40, fraction=0.05, shrink=0.8, ax=ax1)
    cbar.set_label(var_name + ' [' + unit + ']', fontsize=15)
    cbar.ax.tick_params(labelsize=15)
    ax1.set_title('SSH - HH: ' + str(idx_time).zfill(3), fontsize=18, fontweight='bold')

    # Second subplot: SSH along Gibraltar parallel (Lat: 36.0°N)
    ax2 = axes[1]
    gib_lat_idx = 138  # Known index for 36.02°N in this model grid
    ssh_gibraltar = np.squeeze(wrk[gib_lat_idx, :])
    #ax2.plot(nemo_lon[0, :], ssh_gibraltar, color='b', linewidth=2)

    # Plot and fill under the SSH curve
    ax2.set_ylim(-0.15,0.15)
    ax2.set_xlim(-6.0,36.3)
    ax2.plot(nemo_lon[0, :], ssh_gibraltar, color='b', linewidth=2)
    ymin = ax2.get_ylim()[0]  # Get the minimum y-value of the axis
    ax2.fill_between(nemo_lon[0, :], ssh_gibraltar, ymin, color='blue', alpha=0.3)
    ax2.set_title(f"SSH along Lat: 36°N - HH: {str(idx_time).zfill(3)}", fontsize=18, fontweight='bold')
    ax2.set_xlabel('Longitude', fontsize=15)
    ax2.set_ylabel(f'{var_name} [' + unit + ']', fontsize=15)
    ax2.grid(True)
    ax2.tick_params(axis='both', labelsize=13)

    # Save figure
    plt.tight_layout()
    plt.savefig(outfile_temp + str(idx_time).zfill(3) + '_Med_Gibraltar.png')
    plt.close(fig)
