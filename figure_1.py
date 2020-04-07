from pylab import * 
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
#from IPython.display import set_matplotlib_formats
#set_matplotlib_formats('png', 'pdf')
import glob
import os
import cartopy.crs as ccrs

### load datasets ###
os.chdir("/nfs3m/archive/sfa_cache09/users/g00/imitevsk/E2.1_CO2_runs/pytropd/output")
t_05 = xr.open_dataset('temp_0.5_surface.nc')
t_1 = xr.open_dataset('temp_1_surface.nc')
t_15 = xr.open_dataset('temp_1.5_surface.nc')
t_2 = xr.open_dataset('temp_2_surface.nc')
t_3 = xr.open_dataset('temp_3_surface.nc')
t_4 = xr.open_dataset('temp_4_surface.nc')
t_5 = xr.open_dataset('temp_5_surface.nc')
t_6 = xr.open_dataset('temp_6_surface.nc')
t_7 = xr.open_dataset('temp_7_surface.nc')
t_8 = xr.open_dataset('temp_8_surface.nc')
t_1_qf = xr.open_dataset('temp_1_1901_1960_surface_QFLUX.nc')
t_2_qf = xr.open_dataset('temp_2_1901_1960_surface_QFLUX.nc')
t_3_qf = xr.open_dataset('temp_3_1901_1960_surface_QFLUX.nc')
t_4_qf = xr.open_dataset('temp_4_1901_1960_surface_QFLUX.nc')
os.chdir("/nfs3m/archive/sfa_cache09/users/g00/imitevsk/E2.1_CO2_runs/pytropd/plots_figures")

def t_surf_timeseries_pd():
	''' 
	Plots Figure 1 in the paper
	'''

	def y(l_s, l_e):
		"""
		Surface temperature deviations from PI control for fully coupled ocean runs
		
		:Input:
		 - *l_s* (int) - starting latitude
		 - *l_e* (int) - ending latitude
		:Output:
		 - *y* (ndarray) - array with values for 0.5,1,1.5,2,3,4,5,6,7,8xCO2
		"""
		MEAN = ((t_1.temp * t_1.axyp) / t_1.axyp.mean()).sel(year = slice(1950, 2000), lat = slice(l_s, l_e)).mean()
		y = np.array([
			((t_05.temp * t_05.axyp) / t_05.axyp.mean()).sel(year = slice(1950, 2000), lat = slice(l_s, l_e)).mean() - MEAN,	
			((t_1.temp * t_1.axyp) / t_1.axyp.mean()).sel(year = slice(1950, 2000), lat = slice(l_s, l_e)).mean() - MEAN,
			((t_15.temp * t_15.axyp) / t_15.axyp.mean()).sel(year = slice(1950, 2000), lat = slice(l_s, l_e)).mean() - MEAN,
			((t_2.temp * t_2.axyp) / t_2.axyp.mean()).sel(year = slice(1950, 2000), lat = slice(l_s, l_e)).mean() - MEAN,
			((t_3.temp * t_3.axyp) / t_3.axyp.mean()).sel(year = slice(1950, 2000), lat = slice(l_s, l_e)).mean() - MEAN,
			((t_4.temp * t_4.axyp) / t_4.axyp.mean()).sel(year = slice(1950, 2000), lat = slice(l_s, l_e)).mean() - MEAN,
			((t_5.temp * t_5.axyp) / t_5.axyp.mean()).sel(year = slice(1950, 2000), lat = slice(l_s, l_e)).mean() - MEAN,
			((t_6.temp * t_6.axyp) / t_6.axyp.mean()).sel(year = slice(1950, 2000), lat = slice(l_s, l_e)).mean() - MEAN,
			((t_7.temp * t_7.axyp) / t_7.axyp.mean()).sel(year = slice(1950, 2000), lat = slice(l_s, l_e)).mean() - MEAN,
			((t_8.temp * t_8.axyp) / t_8.axyp.mean()).sel(year = slice(1950, 2000), lat = slice(l_s, l_e)).mean() - MEAN ])
		return y
			
	def y_QFLUX(y_s, y_e, l_s, l_e):
		"""
		Surface temperature deviations from PI control for QFLUX runs 
		
		:Input:
		 - *y_s* (int) - starting year 
		 - *y_e* (int) - ending year
		 - *l_s* (int) - starting latitude
		 - *l_e* (int) - ending latitude
		:Output:
		 - *y* (ndarray) - array with values for 0.5,1,1.5,2,3,4,5,6,7,8xCO2
		"""
		MEAN = ((t_1_qf.temp * t_1_qf.axyp) / t_1_qf.axyp.mean()).sel(year = slice(y_s, y_e), lat = slice(l_s, l_e)).mean()
		y = np.array([
			((t_1_qf.temp * t_1_qf.axyp) / t_1_qf.axyp.mean()).sel(year = slice(y_s, y_e), lat = slice(l_s, l_e)).mean() - MEAN,
			((t_2_qf.temp * t_2_qf.axyp) / t_2_qf.axyp.mean()).sel(year = slice(y_s, y_e), lat = slice(l_s, l_e)).mean() - MEAN,
			((t_3_qf.temp * t_3_qf.axyp) / t_3_qf.axyp.mean()).sel(year = slice(y_s, y_e), lat = slice(l_s, l_e)).mean() - MEAN,
			((t_4_qf.temp * t_4_qf.axyp) / t_4_qf.axyp.mean()).sel(year = slice(y_s, y_e), lat = slice(l_s, l_e)).mean() - MEAN ])
		return y
	
	### masking NAWH ###

	td_150 = (t_15.temp-t_1.temp).sel(year=slice(1950,2000)).mean(dim=['year','month']) > 0
	td_150[:45,:] = True
	
	td_20 = (t_2.temp-t_1.temp).sel(year=slice(1950,2000)).mean(dim=['year','month']) > 0
	td_20[:45,:] = True
		
	td_30 = (t_3.temp-t_1.temp).sel(year=slice(1950,2000)).mean(dim=['year','month']) > 0
	td_30[:45,:] = True
	
	td_40 = (t_4.temp-t_1.temp).sel(year=slice(1950,2000)).mean(dim=['year','month']) > 0
	td_40[:45,:] = True
	
	td_50 = (t_5.temp-t_1.temp).sel(year=slice(1950,2000)).mean(dim=['year','month']) > 0
	td_50[:45,:] = True
	
	td_60 = (t_6.temp-t_1.temp).sel(year=slice(1950,2000)).mean(dim=['year','month']) > 0 
	td_60[:45,:] = True
	
	td_70 = (t_7.temp-t_1.temp).sel(year=slice(1950,2000)).mean(dim=['year','month']) > 0 
	td_70[:45,:] = True
	
	td_80 = (t_8.temp-t_1.temp).sel(year=slice(1950,2000)).mean(dim=['year','month']) > 0 
	td_80[:45,:] = True
	
	def y_mask(l_s, l_e):
		"""
		Surface temperature deviations from PI control with NAWH mask applied
		
		:Input:
		 - *l_s* (int) - starting latitude
		 - *l_e* (int) - ending latitude
		:Output:
		 - *y* (ndarray) - array with values for 0.5,1,1.5,2,3,4,5,6,7,8xCO2
		"""	
		MEAN = ((t_1.temp * t_1.axyp) / t_1.axyp.mean()).sel(year = slice(1950, 2000), lat = slice(l_s, l_e)).mean()
		y = np.array([	
			((t_05.temp * t_05.axyp) / t_05.axyp.mean()).sel(year = slice(1950, 2000), lat = slice(l_s, l_e)).mean() - MEAN,
			((t_1.temp * t_1.axyp) / t_1.axyp.mean()).sel(year = slice(1950, 2000), lat = slice(l_s, l_e)).mean() - MEAN,
			((t_15.temp * t_15.axyp) / t_15.axyp.mean()).sel(year = slice(1950, 2000), lat = slice(l_s, l_e)).mean(dim=['month','year']).where(td_150 == True).mean() - MEAN,
			((t_2.temp * t_2.axyp) / t_2.axyp.mean()).sel(year = slice(1950, 2000), lat = slice(l_s, l_e)).mean(dim=['month','year']).where(td_20 == True).mean() - MEAN,
			((t_3.temp * t_3.axyp) / t_3.axyp.mean()).sel(year = slice(1950, 2000), lat = slice(l_s, l_e)).mean(dim=['month','year']).where(td_30 == True).mean() - MEAN,
			((t_4.temp * t_4.axyp) / t_4.axyp.mean()).sel(year = slice(1950, 2000), lat = slice(l_s, l_e)).mean(dim=['month','year']).where(td_40 == True).mean() - MEAN,
			((t_5.temp * t_5.axyp) / t_5.axyp.mean()).sel(year = slice(1950, 2000), lat = slice(l_s, l_e)).mean(dim=['month','year']).where(td_50 == True).mean() - MEAN,
			((t_6.temp * t_6.axyp) / t_6.axyp.mean()).sel(year = slice(1950, 2000), lat = slice(l_s, l_e)).mean(dim=['month','year']).where(td_60 == True).mean() - MEAN,
			((t_7.temp * t_7.axyp) / t_7.axyp.mean()).sel(year = slice(1950, 2000), lat = slice(l_s, l_e)).mean(dim=['month','year']).where(td_70 == True).mean() - MEAN,
			((t_8.temp * t_8.axyp) / t_8.axyp.mean()).sel(year = slice(1950, 2000), lat = slice(l_s, l_e)).mean(dim=['month','year']).where(td_80 == True).mean() - MEAN ])
		return y
	
	
	x = np.array([0.5,1,1.5,2,3,4,5,6,7,8])	
	x = 5.35*np.log(x)	
	x_QFLUX = np.array([1,2,3,4])	
	x_QFLUX = 5.35*np.log(x_QFLUX)	
			
	fig = plt.figure()
	fig.set_figwidth(fig.get_figwidth() * 1.8)		
	axes = fig.add_subplot(1,2,1)	
	axes.errorbar(x, y(-90, 90), marker='o', markersize = 8, linewidth=3, color = 'green', label='global')	
	axes.errorbar(x, y(-90, 0), marker='o', markersize = 8, linewidth=3, color = 'blue', label='SH')	
	axes.errorbar(x, y(0, 90), marker='o', markersize = 8, linewidth=3, color = 'red', label='NH')
	plt.xlabel('Radiative Forcing (W/m$^2$)', fontsize = 20)
	plt.ylabel('$\Delta T_{s}$ (K)', fontsize = 20)
	axes.tick_params(axis='both', labelsize=15)	
	plt.title('A) Surface Temperature', fontsize = 20)
	axes.legend(loc = 2, fontsize = 12)
	plt.yticks([0,4,8])
	plt.xticks([-4,0,4,8,12])
	plt.ylim([-3.1,8.1])
	plt.xlim([-4.2,12])
	
	axes = fig.add_subplot(1,2,2)
	axes.errorbar(x, y(0, 90), marker='o', markersize = 8, linewidth=3, color = 'red', label='NH')
	axes.errorbar(x, y_mask(0,90), marker='o', markersize = 8, linewidth=3, color = 'black', label='NH w/o NAWH')
	axes.errorbar(x_QFLUX, y_QFLUX(1950, 1960, 0, 90), marker='o', markersize = 8, linewidth=3, color = 'cyan', label='NH, QFLUX')
	plt.legend(loc=2)
	plt.xlabel('Radiative Forcing (W/m$^2$)', fontsize = 20)	
	axes.tick_params(axis='both', labelsize=15)
	plt.title('B) Surface Temperature w/o NAWH', fontsize = 20)
	axes.legend(loc = 2, fontsize = 10)
	plt.yticks([0,4,8])
	plt.xticks([-4,0,4,8,12])
	plt.ylim([-3.1,8.1])
	plt.xlim([-4.2,12])
	fig.tight_layout()
	plt.savefig('figure_1_pd3.pdf')
	plt.show()

t_surf_timeseries_pd()
