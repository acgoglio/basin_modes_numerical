#!/usr/bin/env python
from netCDF4 import Dataset
import sys
import netCDF4 as ncdf
import matplotlib as mpl # Palettes
from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt
from copy import copy, deepcopy
import glob
#import astropy.time
import dateutil.parser
import rt_stats_tools
import rt_plotbox
import warnings
from matplotlib.colors import ListedColormap
from matplotlib.colors import LogNorm
warnings.filterwarnings("ignore")
from point_ini import *
mpl.use('Agg')
####################################
exp_name  = "Bathymetry"
NEMO_GRID = mesh_mask
var       = "Bathymetry"
infile    = bathy_meter
outfile   = work_dir+"bathy_points.png"
minV      = "-0.01" #"-0.1"
maxV      = "2500" #"0.1"
thV       = "nan"
cmap      = "Blues" #"YlGnBu"
inidate   = "20150101"
hour      = "00"
enddate   = "20150101"
var_name  = "Bathymetry"
unit      = "m"
index_file = coo_file

# Read the grid
nc_nemogrid = ncdf.Dataset(NEMO_GRID,'r')
try:
   # EAS7 and previous 
   nemo_lat = nc_nemogrid.variables['nav_lat'][:]
   nemo_lon = nc_nemogrid.variables['nav_lon'][:]
   nemo_mask = (nc_nemogrid.variables['tmask'][0,0,:,:] == 0)
   nemo_mbathy = np.squeeze(nc_nemogrid.variables['mbathy'][0,:,:])
except:
   nemo_lat = nc_nemogrid.variables['y'][:]
   nemo_lon = nc_nemogrid.variables['x'][:]
   nemo_mask = (nc_nemogrid.variables['tmask'][0,0,:,:] == 0)
   nemo_mbathy = np.squeeze(nc_nemogrid.variables['mbathy'][0,:,:])
nc_nemogrid.close()

# Read the points indexes

grid_indices = np.loadtxt(index_file, dtype={'names': ('i', 'j', 'name'), 'formats': (int, int, 'U20')})

i_list = [row[0] for row in grid_indices]
j_list = [row[1] for row in grid_indices]
names = [row[2] for row in grid_indices]

lat_points = nemo_lat[j_list, i_list]
lon_points = nemo_lon[j_list, i_list]


# Read the field 
print ('Infile',infile)
fh1=ncdf.Dataset(infile,'r')
bat1 = fh1.variables[var][:]
#var_name = fh1.variables[var].standard_name
#unit = fh1.variables[var].units
fh1.close()

nemo_mbathy = nemo_mbathy-1
nemo_mbathy = np.where(nemo_mbathy<0,0,nemo_mbathy)

# Take the sfc and the bottom if 3D
dims=np.shape(bat1)
# 3D fields
if np.size(dims) > 3:
   # SFC
   wrk=np.squeeze(bat1[0,0,:,:])
   # BOTTOM
   wrk_bot=np.zeros(np.shape(wrk)) 
   for idx_x in range (0,len(bat1[0,0,:,0])):
       for idx_y in range (0,len(bat1[0,0,0,:])):
           wrk_bot[idx_x,idx_y]=np.squeeze(bat1[0,nemo_mbathy[idx_x,idx_y],idx_x,idx_y])
else:
   # 2D fields
   wrk=np.squeeze(bat1[:,:])


# Mask the land
mask_Med = ~rt_stats_tools.get_med_mask(nemo_lon,nemo_lat,~nemo_mask)
wrk = np.ma.array(wrk,mask=mask_Med)

#######################
# PLOT 
# Whole Mediterranean Sea

# Build the plot palettes
try:
   step=(float(maxV)-float(minV))/50 #25
   pa = np.arange(float(minV),float(maxV),step)
except:
   print ("WARNING: something wrong with this field.. I am fixing -100,100 as min max palette values")
   step=25
   pa = np.arange(round(float(-100),2),round(float(100),2),step)

vmin = 0.1 #np.nanmin(wrk[wrk > 0])  
vmax = float(maxV) #np.nanmax(wrk)

plt.figure(figsize=(20, 10))
m = Basemap(projection='cyl',llcrnrlat=30.2,urcrnrlat=46.,\
            llcrnrlon=-6.0,urcrnrlon=36.3,resolution='h')
m.drawparallels(np.arange(30.,46.,5.), labels=[1,0,0,0], fontsize=20)
m.drawmeridians(np.arange(-18.,36.3,5.), labels=[0,0,0,1], fontsize=20)

CS   = plt.contourf(nemo_lon,nemo_lat,np.squeeze(wrk), pa, cmap=cmap, extend='max') #, norm=LogNorm(vmin=vmin, vmax=vmax))
if thV != 'nan' :
   CSth = plt.contour(nemo_lon,nemo_lat,np.squeeze(wrk), [thV], color='green')
m.drawcoastlines(linewidth=1.0, linestyle='solid', color='k', antialiased=1, ax=None, zorder=None)
m.drawmapboundary(fill_color='white')
m.fillcontinents()
cbar=plt.colorbar(CS,orientation='horizontal',aspect=40,fraction=0.1,shrink=0.8)
cbar.set_label( var_name+' [' + unit +']',fontsize=20)
cbar.ax.tick_params(labelsize=20)
if thV != 'nan' :
   cbar.add_lines(CSth)
# Add points
for lon, lat, name in zip(lon_points, lat_points, names):
    xp, yp = m(lon, lat)
    plt.scatter(xp, yp, marker='o', color='red', s=200, alpha=1, edgecolors='black', zorder=10) 

#plt.title( var_name+' '+exp_name+' '+inidate[0:8]+' HH:'+hour, fontsize=20, fontweight='bold')
plt.savefig(outfile)

