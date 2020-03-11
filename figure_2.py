from pylab import * 
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
import matplotlib
#from IPython.display import set_matplotlib_formats
#set_matplotlib_formats('png', 'pdf')
import glob
import os
import cartopy.crs as ccrs
from mpl_toolkits.axes_grid1 import AxesGrid
from cartopy.mpl.geoaxes import GeoAxes

os.chdir("/nfs3m/archive/sfa_cache09/users/g00/imitevsk/E2.1_CO2_runs/pytropd/output")

t_05 = xr.open_dataset('temp_0.5_surface.nc').temp
t_1 = xr.open_dataset('temp_1_surface.nc').temp
t_15 = xr.open_dataset('temp_1.5_surface.nc').temp
t_2 = xr.open_dataset('temp_2_surface.nc').temp
t_3 = xr.open_dataset('temp_3_surface.nc').temp
t_4 = xr.open_dataset('temp_4_surface.nc').temp
t_5 = xr.open_dataset('temp_5_surface.nc').temp
t_6 = xr.open_dataset('temp_6_surface.nc').temp
t_8 = xr.open_dataset('temp_8_surface.nc').temp

os.chdir("/nfs3m/archive/sfa_cache09/users/g00/imitevsk/E2.1_CO2_runs/pytropd/plots_figures")

def t_surf_globe_pd2_r1():
	'''	
	Figure 1 in Paper
	'''

	# defining cbar
	cmap_var = 'RdBu_r'
	cbar_begin = -12
	cbar_end = 13
	cbar_step = 1
	
	lat = t_1.lat 
	lon = t_1.lon
	dt_05 = (t_05-t_1)[-50:,:,:,:].mean(dim='month').mean(dim='year')
	dt_2 = (t_2-t_1)[-50:,:,:,:].mean(dim='month').mean(dim='year')
	dt_3 = (t_3-t_1)[-50:,:,:,:].mean(dim='month').mean(dim='year')
	dt_4 = (t_4-t_1)[-50:,:,:,:].mean(dim='month').mean(dim='year')
	dt_6 = (t_6-t_1)[-50:,:,:,:].mean(dim='month').mean(dim='year')
	dt_8 = (t_8-t_1)[-50:,:,:,:].mean(dim=['month','year'])
	
	fig = plt.figure(figsize=(5.91496607667, 3))
	fig.set_figwidth(fig.get_figwidth() * 3)
	fig.set_figheight(fig.get_figheight() * 2)

	# Set projection
	proj = ccrs.Robinson()
	
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
	
	axgr[0].coastlines()
	axgr[1].coastlines()
	axgr[2].coastlines()
	axgr[3].coastlines()
	axgr[4].coastlines()
	axgr[5].coastlines()
		
	axgr[0].contourf(lon,lat,dt_05,transform=ccrs.PlateCarree(),levels=np.arange(cbar_begin,cbar_end,cbar_step),cmap = cmap_var)	
	axgr[1].contourf(lon,lat,dt_2,transform=ccrs.PlateCarree(), levels=np.arange(cbar_begin,cbar_end,cbar_step),cmap = cmap_var)
	axgr[2].contourf(lon,lat,dt_3,transform=ccrs.PlateCarree(),levels=np.arange(cbar_begin,cbar_end,cbar_step),cmap = cmap_var)
	axgr[3].contourf(lon,lat,dt_4,transform=ccrs.PlateCarree(),levels=np.arange(cbar_begin,cbar_end,cbar_step),cmap = cmap_var)
	axgr[4].contourf(lon,lat,dt_6,transform=ccrs.PlateCarree(),levels=np.arange(cbar_begin,cbar_end,cbar_step),cmap = cmap_var, extend='max')
	axgr[5].contourf(lon,lat,dt_8,transform=ccrs.PlateCarree(),levels=np.arange(cbar_begin,cbar_end,cbar_step),cmap = cmap_var, extend='max')	
		
	cf = axgr[1].contourf(lon,lat,dt_2,transform=ccrs.PlateCarree(), levels=np.arange(cbar_begin,cbar_end,cbar_step),cmap = cmap_var)
	cb = axgr.cbar_axes[0].colorbar(cf, ticks=[-12,-6,0,6,12])
		
	axgr.cbar_axes[0].axis[axgr.cbar_axes[0].orientation].label.set_text(r'K')
	axgr.cbar_axes[0].axis[axgr.cbar_axes[0].orientation].label.set_size(20) 
	cb.ax.tick_params(labelsize=20) 
	
	axgr[0].set_title('A) 0.5xCO$_2$ - PI',fontsize=25)
	axgr[1].set_title('B) 2xCO$_2$ - PI',fontsize=25)
	axgr[2].set_title('C) 3xCO$_2$ - PI',fontsize=25)
	axgr[3].set_title('D) 4xCO$_2$ - PI',fontsize=25)
	axgr[4].set_title('E) 6xCO$_2$ - PI',fontsize=25)
	axgr[5].set_title('F) 8xCO$_2$ - PI',fontsize=25)
	
	plt.tight_layout(pad=1.5, w_pad=1.0, h_pad=0.5) 	
	plt.savefig('t_surf_globe_pd2_r1.pdf')	
	plt.show()	

t_surf_globe_pd2_r1()
