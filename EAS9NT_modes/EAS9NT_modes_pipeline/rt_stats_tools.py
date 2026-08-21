
import numpy as np
import matplotlib.pylab as plt

def rms(X,dim):
   """
   Root mean square of an array along one dimension
   """
   siz = np.shape(X)
   my_rms=np.sqrt(np.sum(np.square(X),axis=dim))/np.sqrt(siz[dim]-1)
   return my_rms

def loess(data,fc,t=np.nan,t_final=np.nan,step=1):
   """
   Loess filtering of a time serie
   """
   # Initialize
   if np.isnan(t):
      t       = np.arange(0,len(data)*step,step)
   if np.isnan(t_final):
      t_final = np.arange(0,len(data)*step,step)
   tau = 1./fc   
   data_smooth = np.ones(len(t_final))*np.nan

   # Only compute for the points where t_final is in the range of t
   sx = (t_final >= np.min(t)) & (t_final <= np.max(t))

   # Loop on desired points
   for i in np.arange(len(t_final))[sx]:
      dn = np.abs(t-t_final[i])/tau
      idx_weights = dn<1
      n_pts = idx_weights.sum()
      if (n_pts > 3):
         dn = dn[idx_weights]
         w = (1-dn*dn*dn)
         weights = w*w*w
         X = np.array([np.ones(n_pts),t[idx_weights],t[idx_weights]**2]).transpose()
         W = np.diag(weights)
         B,resid,rank,s = np.linalg.lstsq(np.dot(W,X),np.dot(W,data[idx_weights]))
         data_smooth[i] = B[0]+B[1]*t_final[i]+B[2]*t_final[i]**2
   return data_smooth

def gc_dist(lon1,lat1,lon2,lat2):
   """
   Distance between 2 vector points along a great circle
   """
   lon1 = (np.pi/180.) * lon1
   lon2 = (np.pi/180.) * lon2
   lat1 = (np.pi/180.) * lat1
   lat2 = (np.pi/180.) * lat2

   dlat = lat2-lat1
   dlon = lon2-lon1

   #haversine function
   dang = 2*np.arcsin( np.sqrt( np.sin(dlat/2)**2 + np.cos(lat2)*np.cos(lat1)*np.sin(dlon/2)**2 ) )

   r_earth = 6371315.

   return r_earth*dang

def gc_dist_diff(lon,lat):
   """
   Distance between the points of a vector points along a great circle
   """
   lon1 = lon[0:-1]
   lat1 = lat[0:-1]
   lon2 = lon[1:]
   lat2 = lat[1:]

   return gc_dist(lon1,lat1,lon2,lat2)

def get_med_mask(lon,lat,mask):
   # Region
   geo = [-5.5,36.4,30.2,45.9]
   mask_Med = np.logical_and(np.logical_and(lon>geo[0],lon<geo[1]),np.logical_and(lat>geo[2],lat<geo[3]))
   # Remove Bay of Biscay
   geo = [-7,-.5,42.9,46]
   mask_BoB = np.logical_and(np.logical_and(lon>geo[0],lon<geo[1]),np.logical_and(lat>geo[2],lat<geo[3]))
   mask_Med = np.logical_and(mask_Med,~mask_BoB)
   # Remove Black Sea
   geo = [26.9,44,40.1,48]
   mask_BS = np.logical_and(np.logical_and(lon>geo[0],lon<geo[1]),np.logical_and(lat>geo[2],lat<geo[3]))
   mask_Med = np.logical_and(mask_Med,~mask_BS)
   # Combine masks
   return np.logical_and(mask_Med,mask)


def get_med_regmasks(lon,lat,lon_def,lat_def,mask_regs_def,verbose=False):
   """
   Get masks of the differents area of Med from the NEMO masks
   """
   from matplotlib import path
   pts_grid = np.concatenate((np.reshape(lon,(-1,1)),np.reshape(lat,(-1,1))),axis=1)
   N_regs = np.int(np.max(mask_regs_def))
   lon_ext,lat_ext,mask_regs_ext = add_border_to_matrix(lon_def,lat_def,mask_regs_def,18,1)
   mask_regs_out = np.ones((lon.shape[0],lon.shape[1],N_regs),dtype=bool)
   fig = plt.figure()
   for i_mask in range(0,N_regs):
      if verbose:
         print (" -- Region "+str(i_mask))
      # Get contour
      cont_t = plt.contour(lon_ext,lat_ext,mask_regs_ext==i_mask,[0.5]).allsegs[0]
      if (len(cont_t)>0):
         # Create path
         p = path.Path(cont_t[0])
         # Check if grid points are inside the path
         mask_tmp = p.contains_points(pts_grid)
         # Reshape onto the 2D
         mask_regs_out[:,:,i_mask] = np.reshape(mask_tmp,[lon.shape[0],lon.shape[1]])
   plt.close(fig)
   return mask_regs_out

def get_med_regconts(lon_def,lat_def,mask_regs_def):
   N_regs = np.int(np.max(mask_regs_def))
   lon_ext,lat_ext,mask_regs_ext = add_border_to_matrix(lon_def,lat_def,mask_regs_def,18,1)
   cont_regs = []
   fig = plt.figure()
   for i_mask in range(0,N_regs):
      # Get contour
      cont_t = plt.contour(lon_ext,lat_ext,mask_regs_ext==i_mask,[0.5]).allsegs[0]
      if (len(cont_t)>0):
         cont_regs.append(cont_t[0])
      else:
         cont_regs.append([])
   plt.close(fig)
   return cont_regs


def add_border_to_matrix(lon,lat,mat,val,n_pts,axis=2):
   """
   Add a border with value "val" for "n_pts" around the matrix "mat" with updated coordinates
   """
   n_pts_x = n_pts
   n_pts_y = n_pts
   if (axis == 0):
      n_pts_y = 0
   elif (axis == 1):
      n_pts_x = 0


   # Create the empty arrays
   Nlat,Nlon = mat.shape
   mat_out = np.ones((Nlat+n_pts_y*2,Nlon+n_pts_x*2),dtype=mat.dtype)*val
   lon_out = np.zeros((Nlat+n_pts_y*2,Nlon+n_pts_x*2),dtype=lon.dtype)
   lat_out = np.zeros((Nlat+n_pts_y*2,Nlon+n_pts_x*2),dtype=lat.dtype)
   # Fill the center
   if (n_pts_x == 0):
      mat_out[n_pts_y:-n_pts_y,:] = mat
      lon_out[n_pts_y:-n_pts_y,:] = lon
      lat_out[n_pts_y:-n_pts_y,:] = lat
   elif (n_pts_y == 0):
      mat_out[:,n_pts_x:-n_pts_x] = mat
      lon_out[:,n_pts_x:-n_pts_x] = lon
      lat_out[:,n_pts_x:-n_pts_x] = lat
   else:
      mat_out[n_pts_y:-n_pts_y,n_pts_x:-n_pts_x] = mat
      lon_out[n_pts_y:-n_pts_y,n_pts_x:-n_pts_x] = lon
      lat_out[n_pts_y:-n_pts_y,n_pts_x:-n_pts_x] = lat
   # Extrapolate the coordinates
   # First along the x-axis
   if (n_pts_x != 0):
      dlon1 = lon[:,1]-lon[:,0]
      dlon2 = lon[:,-1]-lon[:,-2]
      dlat1 = lat[:,1]-lat[:,0]
      dlat2 = lat[:,-1]-lat[:,-2]
      for i_pt in range(n_pts_x-1,-1,-1):
         if (n_pts_y != 0):
            lon_out[n_pts_y:-n_pts_y,i_pt]    = lon_out[n_pts_y:-n_pts_y,i_pt+1]-dlon1
            lon_out[n_pts_y:-n_pts_y,-i_pt-1] = lon_out[n_pts_y:-n_pts_y,-i_pt-2]+dlon2
            lat_out[n_pts_y:-n_pts_y,i_pt]    = lat_out[n_pts_y:-n_pts_y,i_pt+1]-dlat1
            lat_out[n_pts_y:-n_pts_y,-i_pt-1] = lat_out[n_pts_y:-n_pts_y,-i_pt-2]+dlat2
         else:
            lon_out[:,i_pt]    = lon_out[:,i_pt+1]-dlon1
            lon_out[:,-i_pt-1] = lon_out[:,i_pt-2]+dlon2
            lat_out[:,i_pt]    = lat_out[:,i_pt+1]-dlat1
            lat_out[:,-i_pt-1] = lat_out[:,i_pt-2]+dlat2
   # Then along the y-axis
   if (n_pts_y != 0):
      dlon1 = lon_out[n_pts_y+1,:]-lon_out[n_pts_y,:]
      dlon2 = lon_out[-n_pts_y-1,:]-lon_out[-n_pts_y-2,:]
      dlat1 = lat_out[n_pts_y+1,:]-lat_out[n_pts_y,:]
      dlat2 = lat_out[-n_pts_y-1,:]-lat_out[-n_pts_y-2,:]
      for i_pt in range(n_pts_y-1,-1,-1):
         lon_out[i_pt,:]    = lon_out[i_pt+1,:]-dlon1
         lon_out[-i_pt-1,:] = lon_out[-i_pt-2,:]+dlon2
         lat_out[i_pt,:]    = lat_out[i_pt+1,:]-dlat1
         lat_out[-i_pt-1,:] = lat_out[-i_pt-2,:]+dlat2
 
   return lon_out,lat_out,mat_out




