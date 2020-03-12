from pylab import * 
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
#plt.rc('text', usetex=True)
plt.rc('font', family='serif')
#from IPython.display import set_matplotlib_formats
#set_matplotlib_formats('png', 'pdf')
import glob
import os
import cartopy
import cartopy.crs as ccrs
from mpl_toolkits.axes_grid1 import AxesGrid
from cartopy.mpl.geoaxes import GeoAxes
import matplotlib.path as mpath

### load datasets ###

os.chdir("/nfs3m/archive/sfa_cache09/users/g00/imitevsk/E2.1_CO2_runs/pytropd/output")
i_05 = xr.open_dataset('temp_0.5_surface.nc').ice
i_1 = xr.open_dataset('temp_1_surface.nc').ice
i_2 = xr.open_dataset('temp_2_surface.nc').ice
i_3 = xr.open_dataset('temp_3_surface.nc').ice
i_4 = xr.open_dataset('temp_4_surface.nc').ice
i_5 = xr.open_dataset('temp_5_surface.nc').ice
i_6 = xr.open_dataset('temp_6_surface.nc').ice
i_8 = xr.open_dataset('temp_8_surface.nc').ice
os.chdir("/nfs3m/archive/sfa_cache09/users/g00/imitevsk/E2.1_CO2_runs/pytropd/plots_figures")

def ice_nh_all_year_pd():
	'''
	Plots Figure 3 in the paper
	'''	
	# cbar variables
	cmap_var = 'seismic'
	cbar_start = -100
	cbar_end = 105
	cbar_step = 5
		
	lat = i_1.lat 
	lon = i_1.lon
	dlon = lon[1] - lon[0]
	lon = np.concatenate((lon, lon[-1:] + dlon)) # needed to avoid the empty "line" at 180 lon.
	
	months = ['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']
	
	di_05 = np.array((i_05-i_1)[-50:,:,:,:].sel(month = months).mean(dim=['month','year']))
	di_2 = np.array((i_2-i_1)[-50:,:,:,:].sel(month = months).mean(dim=['month','year']))
	di_3 = np.array((i_3-i_1)[-50:,:,:,:].sel(month = months).mean(dim=['month','year']))
	di_4 = np.array((i_4-i_1)[-50:,:,:,:].sel(month = months).mean(dim=['month','year']))
	di_6 = np.array((i_6-i_1)[-50:,:,:,:].sel(month = months).mean(dim=['month','year']))
	di_8 = np.array((i_8-i_1)[-50:,:,:,:].sel(month = months).mean(dim=['month','year']))
	
	# needed to avoid the empty "line" at 180 lon.
	di_05 = np.concatenate((di_05, di_05[:, 0:1]), axis=1)
	di_2 = np.concatenate((di_2, di_2[:, 0:1]), axis=1)
	di_3 = np.concatenate((di_3, di_3[:, 0:1]), axis=1)
	di_4 = np.concatenate((di_4, di_4[:, 0:1]), axis=1)
	di_6 = np.concatenate((di_6, di_6[:, 0:1]), axis=1)
	di_8 = np.concatenate((di_8, di_8[:, 0:1]), axis=1)
	
	fig = plt.figure(figsize=(4, 4))
	fig.set_figwidth(fig.get_figwidth() * 3)
	fig.set_figheight(fig.get_figheight() * 2)

	# Set projection
	proj = ccrs.NorthPolarStereo(globe=None)
	
	# create an axes class using cartopy.
	axes_class = (GeoAxes, dict(map_projection=proj))
	
	# Make an Axes grid using matplotlib
	axgr = AxesGrid(fig, 111, axes_class=axes_class,
	                nrows_ncols=(2, 3),
	                axes_pad=0.4,
	                cbar_location='right',
	                cbar_mode='single',
	                cbar_pad=0.4,
	                cbar_size='4%',
	                label_mode='')
	
	axgr[0].set_title('A) 0.5xCO$_2$ - PI',fontsize=20, usetex=True)
	axgr[1].set_title('B) 2xCO$_2$ - PI',fontsize=20, usetex=True)
	axgr[2].set_title('C) 3xCO$_2$ - PI',fontsize=20, usetex=True)
	axgr[3].set_title('D) 4xCO$_2$ - PI',fontsize=20, usetex=True)
	axgr[4].set_title('E) 6xCO$_2$ - PI',fontsize=20, usetex=True)
	axgr[5].set_title('F) 8xCO$_2$ - PI',fontsize=20, usetex=True)
	
	axgr[0].set_extent([-180,180,55,90], ccrs.PlateCarree())
	axgr[1].set_extent([-180,180,55,90], ccrs.PlateCarree())
	axgr[2].set_extent([-180,180,55,90], ccrs.PlateCarree())
	axgr[3].set_extent([-180,180,55,90], ccrs.PlateCarree())
	axgr[4].set_extent([-180,180,55,90], ccrs.PlateCarree())
	axgr[5].set_extent([-180,180,55,90], ccrs.PlateCarree())

	theta = np.linspace(0, 2*np.pi, 100)
	center = np.array([0.5,0.5])
	radius = 0.5
	verts = np.vstack([np.sin(theta), np.cos(theta)]).T
	circle = mpath.Path(verts * radius + center)
	
	# adding LAND on top of the data
	axgr[0].add_feature(cartopy.feature.LAND, zorder=1, color='gray')
	axgr[1].add_feature(cartopy.feature.LAND, zorder=1, color='gray')
	axgr[2].add_feature(cartopy.feature.LAND, zorder=1, color='gray')
	axgr[3].add_feature(cartopy.feature.LAND, zorder=1, color='gray')
	axgr[4].add_feature(cartopy.feature.LAND, zorder=1, color='gray')
	axgr[5].add_feature(cartopy.feature.LAND, zorder=1, color='gray')
	
	axgr[0].coastlines()
	axgr[1].coastlines()
	axgr[2].coastlines()
	axgr[3].coastlines()
	axgr[4].coastlines()
	axgr[5].coastlines()
	
	axgr[0].gridlines()
	axgr[1].gridlines()
	axgr[2].gridlines()
	axgr[3].gridlines()
	axgr[4].gridlines()
	axgr[5].gridlines()
		
	axgr[0].set_boundary(circle, transform=axgr[0].transAxes)
	axgr[1].set_boundary(circle, transform=axgr[1].transAxes)
	axgr[2].set_boundary(circle, transform=axgr[2].transAxes)
	axgr[3].set_boundary(circle, transform=axgr[3].transAxes)
	axgr[4].set_boundary(circle, transform=axgr[4].transAxes)
	axgr[5].set_boundary(circle, transform=axgr[5].transAxes)
	
	# use matplotlib.pyplot.contourf to great spatial plot on cartopy map
	axgr[0].contourf(lon,lat,di_05,transform=ccrs.PlateCarree(), levels=np.arange(cbar_start, cbar_end, cbar_step), cmap = cmap_var)#extent = img_extent,	
	axgr[1].contourf(lon,lat,di_2,transform=ccrs.PlateCarree(), levels=np.arange(cbar_start, cbar_end, cbar_step),cmap = cmap_var)
	axgr[2].contourf(lon,lat,di_3,transform=ccrs.PlateCarree(), levels=np.arange(cbar_start, cbar_end, cbar_step),cmap = cmap_var)
	axgr[3].contourf(lon,lat,di_4,transform=ccrs.PlateCarree(), levels=np.arange(cbar_start, cbar_end, cbar_step),cmap = cmap_var)
	axgr[4].contourf(lon,lat,di_6,transform=ccrs.PlateCarree(), levels=np.arange(cbar_start, cbar_end, cbar_step),cmap = cmap_var)
	axgr[5].contourf(lon,lat,di_8,transform=ccrs.PlateCarree(), levels=np.arange(cbar_start, cbar_end, cbar_step),cmap = cmap_var)	
	cf = axgr[5].contourf(lon,lat,di_8,transform=ccrs.PlateCarree(), levels=np.arange(cbar_start, cbar_end, cbar_step),cmap = cmap_var)
		
	# we only want one colorbar so we're going to call axgr.cbar_axes[0].colorbar(cf) where cf is one of the plots
	cb = axgr.cbar_axes[0].colorbar(cf, ticks=[-100,-50,0,50,100])
		
	# Change colorbar label, font size, etc
	axgr.cbar_axes[0].axis[axgr.cbar_axes[0].orientation].label.set_text('% Sea-Ice')
	axgr.cbar_axes[0].axis[axgr.cbar_axes[0].orientation].label.set_size(20)
	cb.ax.tick_params(labelsize=20)	
			
	plt.tight_layout()	
	plt.savefig('ice_nh_all_year_pd.pdf', bbox_inches='tight')	
	plt.show()

ice_nh_all_year_pd()
