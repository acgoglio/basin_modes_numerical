# -*- coding: utf-8 -*-
"""
Created on Mon Aug 29 15:41:14 2016

@author: romain
"""

import numpy as np
import matplotlib.colors
import matplotlib.pylab as plt
import matplotlib.cm as cm
from mpl_toolkits.basemap import Basemap
import math
import sys
from decimal import getcontext, Decimal

def rt_plot_2D( lon,lat,var,\
               geo=[-98, -45, 10, 52],clim=[float('nan') ,1],cstep=25,pal=cm.jet,\
               coastcolor=[.9,.9,.9],coastlinecolor=[.5,.5,.5],lakecolor=[1,1,1],
               plotitle='',subpltid=111,fontsize=16,colorbar=True,stepX=10.,stepY=10.,
               ax='none',pcolor=False,coast_resol='l',ccont=[],
               parallel_labels=[True,False,False,True],
               meridian_labels=[False,False,False,True]):
   """
   Plot a map with 2D field
   """    
   getcontext().prec = 2 
   # Initialize
   if (type(var) == np.ndarray) | (type(var) == np.ma.core.MaskedArray):
      if not var.any():
         is_plot=False
      else:
         is_plot=True
   else:
      is_plot=False
   if (math.isnan(clim[0]) and is_plot):
      clim = [np.nanmin(var[:]),np.nanmax(var[:])];
   if not isinstance(cstep,float):
      cstep=(clim[1]-clim[0])/float(cstep)
   # Create figure
   fig=plt.gcf()
   if ax == 'none':
      if isinstance(subpltid,tuple):
         ax  = fig.add_subplot(subpltid[0],subpltid[1],subpltid[2])
      else:
         ax  = fig.add_subplot(subpltid)
   # Create projection
   m = Basemap(projection='cyl',llcrnrlon=geo[0],urcrnrlon=geo[1],\
               llcrnrlat=geo[2],urcrnrlat=geo[3],\
               resolution=coast_resol)
   # Draw coast
   m.drawcoastlines(color=coastlinecolor)
   m.fillcontinents(color=coastcolor,lake_color=lakecolor)
   # Grid
   parallels = np.arange(-90.,90.,stepY)
   m.drawparallels(parallels,labels=parallel_labels,fontsize=fontsize,color=[.6,.6,.6]) # labels = [left,right,top,bottom]
   meridians = np.arange(-180.,180.,stepX)
   m.drawmeridians(meridians,labels=meridian_labels,fontsize=fontsize,color=[.6,.6,.6]) # labels = [left,right,top,bottom]
   # contour filled
   if is_plot:
      if pcolor:
         m.pcolor(lon,lat,var,cmap=pal,vmin=clim[0], vmax=clim[1])
      else:
         norm = matplotlib.colors.Normalize(vmin=clim[0], vmax=clim[1])
         if len(ccont)>0:
            contours = np.array(ccont)
         else:
            contours = np.arange(clim[0],clim[1]+cstep,cstep)
         C = m.contourf(lon,lat,var,contours,cmap=pal,norm=norm,extend='both')
      # colorbar
      if colorbar:
         cbar = plt.colorbar(C,orientation='horizontal',shrink=0.9)
         cbar.ax.tick_params(labelsize=fontsize) 
   # title
   plt.title(plotitle,fontsize=fontsize+2)
   m.ax = ax

   return m

def xaxis_date(fds,step=1,ax=0,date_format='%Y/%m/%d',monthyear_step=1,label_step=1,label_offset=0):
   """
   Put dates on X-axis
   """
   #from matplotlib.ticker import AutoMinorLocator
   import matplotlib.dates as mpl_dates
   try:
      ax.xaxis
   except :
      ax=plt.gca()
   if (step == 'month'):
      datevec = mpl_dates.num2date(fds)
      datetick = []
      for i in xrange(len(datevec)):
         if datevec[i].day == 1:
            if monthyear_step>1:
               if (datevec[i].month % monthyear_step) == 1:
                  datetick.append(datevec[i])
            else:
               datetick.append(datevec[i])
      fdstick = mpl_dates.date2num(datetick)
   elif (step == 'year'):
      datevec = mpl_dates.num2date(fds)
      datetick = []
      for i in xrange(len(datevec)):
         if (datevec[i].day == 1 and datevec[i].month == 1):
            if monthyear_step>1:
               if (datevec[i].year % monthyear_step) == 1:
                  datetick.append(datevec[i])
            else:
               datetick.append(datevec[i])
      fdstick = mpl_dates.date2num(datetick)
   else:
      fdstick = fds[::step]
   ax.xaxis.set_ticks(fdstick)
   hfmt = mpl_dates.DateFormatter(date_format) # matplotlib date format object
   ax.xaxis.set_major_formatter(hfmt)
   # Only plot label every label_step
   if (label_step > 1):
      labels = [item.get_text() for item in ax.xaxis.get_ticklabels()]
      for ilabel in xrange(len(labels)):
         if not (ilabel % label_step) == label_offset:
            labels[ilabel] = ''
      ax.set_xticklabels(labels)
   #minor_locator = AutoMinorLocator(step)
   #ax.xaxis.set_minor_locator(minor_locator)

def nanrms(X,dim):
   """
   Root mean square of an array along one dimension
   """
   siz = np.shape(X)
   my_rms=np.sqrt(np.nansum(np.square(X),axis=dim))/np.sqrt(siz[dim]-1)
   return my_rms


def target_diagram(series_dict, ref_dict, colors=plt.cm.rainbow, markers=['o','*','s','P','X','D'], linestyle='-',markersize=[10,15,10],fontsize=18):
   ''' Plot a target diagram for the series compared to the reference 
       colors can be a dict, a color or a colormap '''

   series_keys = series_dict.keys()
   N_series = len(series_dict)
   N_max = 0
   for mykey in series_keys:
      N_max = max(N_max,series_dict[mykey].shape[1])

   # If colors is colormap, extract colors from it 
   if (type(colors) == matplotlib.colors.ListedColormap) or (type(colors) == matplotlib.colors.LinearSegmentedColormap):
      colors_tmp = get_colors_from_cmap(colors,N_series)
      colors = dict()
      for ikey,mykey in enumerate(series_keys):
         colors[mykey] = colors_tmp[ikey,:]

   # If markers of markersize < size of series, extend them
   if len(markers) < N_max:
      markers = markers + [markers[-1] for i in xrange(N_max-len(markers))]
   if len(markersize) < N_max:
      markersize = markersize + [markersize[-1] for i in xrange(N_max-len(markersize))]

   # Check in series and ref have the same keys
   for mykey in series_keys:
      if mykey not in ref_dict.keys():
         print ('Error! The series and reference need the same keys.')
         return
   # Check in series and colors have the same keys
   if not type(colors) == str:
      for mykey in series_keys:
         if mykey not in colors.keys():
            print ('Error! The series and colors need the same keys.')
            return

   # Compute relevant stats:
   # - Standard deviation
   std_dict = dict()
   stdref_dict = dict()
   for mykey in series_dict.keys():
      stdref_dict[mykey] = np.nanstd(ref_dict[mykey],0)
      std_dict[mykey]    = np.nanstd(series_dict[mykey],0)
   # - Bias
   bias_dict = dict()
   for mykey in series_dict.keys():
      # Special case of 2D field for series
      if (len(series_dict[mykey].shape)>1):
         ref_2D = np.tile(ref_dict[mykey],(series_dict[mykey].shape[1],1)).transpose()
         bias_dict[mykey] = np.nanmean(series_dict[mykey]-ref_2D,0)
      else:
         bias_dict[mykey] = np.nanmean(series_dict[mykey]-ref_dict[mykey],0)
   # - Unbiased RMSD
   rmsd_dict = dict()
   for mykey in series_dict.keys():
      if (len(series_dict[mykey].shape)>1):
         ref_2D  = np.tile(ref_dict[mykey],(series_dict[mykey].shape[1],1)).transpose()
         bias_2D = np.tile(bias_dict[mykey],(series_dict[mykey].shape[0],1))
         rmsd_dict[mykey] = nanrms(series_dict[mykey]-ref_2D-bias_2D,0)
      else:
         rmsd_dict[mykey] = nanrms(series_dict[mykey]-ref_dict[mykey]-bias_dict[mykey],0)

   # Plot data
   fig=plt.gcf()
   for mykey in series_dict.keys():
      for ipoint in xrange(len(rmsd_dict[mykey])):
         if ipoint==0:
            plt.plot(rmsd_dict[mykey][ipoint]*np.sign(std_dict[mykey][ipoint]-stdref_dict[mykey])/stdref_dict[mykey],bias_dict[mykey][ipoint]/stdref_dict[mykey],marker=markers[ipoint],color=colors[mykey],label=mykey,markersize=markersize[ipoint])
         else:
            plt.plot(rmsd_dict[mykey][ipoint]*np.sign(std_dict[mykey][ipoint]-stdref_dict[mykey])/stdref_dict[mykey],bias_dict[mykey][ipoint]/stdref_dict[mykey],marker=markers[ipoint],color=colors[mykey],markersize=markersize[ipoint])
      if not linestyle == 'none':
         plt.plot(rmsd_dict[mykey]*np.sign(std_dict[mykey]-stdref_dict[mykey])/stdref_dict[mykey],bias_dict[mykey]/stdref_dict[mykey],color=colors[mykey])

   # Make the plot nice
   plt.legend(loc='lower right')
   plt.plot([-1,1],[0,0],'k')
   plt.plot([0,0],[-1,1],'k')
   ax = plt.gca()
   ax.set_xlim([-1.4,1.4])
   ax.set_ylim([-1.4,1.4])
   ax.set_xticks(np.arange(-1.3,1.4,0.1));ax.set_xticklabels(['','','','-1','','','','','-.5','','','','','0','','','','','0.5','','','','','1','','',''],fontsize=fontsize)
   ax.set_yticks(np.arange(-1.3,1.4,0.1));ax.set_yticklabels(['','','','-1','','','','','-.5','','','','','0','','','','','0.5','','','','','1','','',''],fontsize=fontsize)
   plt.tick_params(
       axis='both',          # changes apply to the x-axis
       which='both',      # both major and minor ticks are affected
       top='off',         # ticks along the top edge are off
       right='off')       # ticks along the right edge are off

   ax.text(1.02,0.03,'Unbiased',fontsize=fontsize-4)
   ax.text(1.09,-0.07,'RMSD',fontsize=fontsize-4)
   ax.text(0.02,1.15,'Bias',fontsize=fontsize-4)
   ax.spines['left'].set_position(('data', 0.0))
   ax.spines['bottom'].set_position(('data', 0.0))
   ax.spines['right'].set_color('none')
   ax.spines['top'].set_color('none')
   for radi in np.arange(0.1,1.4,0.1):
      circ=plt.Circle((0, 0), radi, fill=False,linewidth=1,linestyle='--',color=[0.8,0.8,0.8]);ax.add_artist(circ)
   circ=plt.Circle((0, 0), 1, fill=False,linewidth=2);ax.add_artist(circ)
   ax.set_aspect('equal')
   
   return ax
   

def rt_getcolormaps(file_cb='/Users/romain/Work/Tools/IMEDEA/Matlab/Matlab-toolbox/Plot/rt_colormaps.txt'):
   rt_colormaps = dict()
   # Open file
   curr = 0
   with open(file_cb,'r') as f:
      for line in f:
         i_tmp=np.mod(curr,4)
         if (i_tmp==0):
            # The first line is the name of the colormap
            if curr>0:
               # If not on the first line of file, save previous cm
               rt_colormaps[name] = cmat2cmpl(val);
               del(val)
            name=line.strip()
         else:
            # get values of colormap
            line=line.strip()
            cols = line.split(',')
            # if Red values (first column), initialize array
            if (i_tmp == 1):
               val=np.zeros([3,len(cols)-1])
            val[i_tmp-1,:]=[float(y) for y in cols[0:-1]]
         curr+=1
   return rt_colormaps


def cmat2cmpl(colormap): 
   """ 
   Convert matlab style colormap to matplotlib style 
   Enter a list non normalized RGB values from 0-255 
   """ 
   r = colormap[0,:]
   g = colormap[1,:]
   b = colormap[2,:]
   cmap = matplotlib.colors.ListedColormap(zip(r,g,b)) 
   return cmap 

def cmpl2cmat(cmap):
   """
   Convert matplotlib style colormap to matlab style
   Enter a matplotlib colormap
   """
   colormap=np.transpose(np.array([list(cmap(i/(cmap.N-1.))) for i in xrange(cmap.N)])[:,0:3])
   return colormap

def inv_cmap(cmap):
   """
   Invert colormap
   """
   return matplotlib.colors.ListedColormap(cmap.colors[::-1])

def set_plt_fontsize(ax,fontsiz):
   """ 
   Change fontsize of given plot
   """ 
   for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
      item.set_fontsize(fontsiz)

def rt_getcolors(file_col='/Users/romain/Work/Tools/IMEDEA/Matlab/Matlab-toolbox/Plot/rt_colors.txt'):
   rt_colors = dict()
   # Open file
   with open(file_col,'r') as f:
      for line in f:
         cols=line.strip().split(',')
         rt_colors[cols[0]] = (float(cols[1]),float(cols[2]),float(cols[3]))
   return rt_colors

def get_colors_from_cmap(cmap,Ncol):
   ''' Get Ncol colors from a colormap cmap '''

   colors = np.zeros((Ncol,3))
   colormap = cmpl2cmat(cmap)
   for i in xrange(3):
      colors[:,i] = np.interp(np.arange(Ncol)*1./(Ncol-1),np.arange(cmap.N)*1./(cmap.N-1),colormap[i,:])
   return colors


def pale_color(color,perc=0.3):
   col_pale=[0,0,0]
   for i in xrange(3):
       inc=(1-color[i])*perc
       col_pale[i]=min(1,color[i]+inc)
   return col_pale


def pale_cmap(cmap,perc=0.3):
   colormap      =cmpl2cmat(cmap)
   colormap_pale = np.empty(colormap.shape)
   for i in xrange(colormap.shape[1]):
       colormap_pale[:,i] = pale_color(colormap[:,i],perc=perc)
   cmap_pale = cmat2cmpl(colormap_pale)
   return cmap_pale 



