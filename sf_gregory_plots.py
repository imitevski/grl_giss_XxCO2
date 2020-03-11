from pylab import * 
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
#from IPython.display import set_matplotlib_formats
#set_matplotlib_formats('png', 'pdf')
import glob
import os
import cartopy.crs as ccrs
import scipy 
import scipy.stats
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
SMALL_SIZE = 15
MEDIUM_SIZE = 20
BIGGER_SIZE = 25
plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=8)    	   	# legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title



def gregory_plot():

	''' loading the data'''	
	
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
	
	t_1_qflux = xr.open_dataset('temp_1_1901_1960_surface_QFLUX.nc')
	t_2_qflux = xr.open_dataset('temp_2_1901_1960_surface_QFLUX.nc')
	t_3_qflux = xr.open_dataset('temp_3_1901_1960_surface_QFLUX.nc')
	t_4_qflux = xr.open_dataset('temp_4_1901_1960_surface_QFLUX.nc')
		
	rh_05 = xr.open_dataset('net_rad_planet_hemis_avg_0.5_xCO2_1850_2000.nc').rad
	rh_1 = xr.open_dataset('net_rad_planet_hemis_avg_1_xCO2_1850_2000.nc').rad
	rh_15 = xr.open_dataset('net_rad_planet_hemis_avg_1.5_xCO2_1850_2000.nc').rad
	rh_2 = xr.open_dataset('net_rad_planet_hemis_avg_2_xCO2_1850_2000.nc').rad
	rh_3 = xr.open_dataset('net_rad_planet_hemis_avg_3_xCO2_1850_2000.nc').rad
	rh_4 = xr.open_dataset('net_rad_planet_hemis_avg_4_xCO2_1850_2000.nc').rad
	rh_5 = xr.open_dataset('net_rad_planet_hemis_avg_5_xCO2_1850_2000.nc').rad
	rh_6 = xr.open_dataset('net_rad_planet_hemis_avg_6_xCO2_1850_2000.nc').rad
	rh_7 = xr.open_dataset('net_rad_planet_hemis_avg_7_xCO2_1850_2000.nc').rad
	rh_8 = xr.open_dataset('net_rad_planet_hemis_avg_8_xCO2_1850_2000.nc').rad

	rh_1_qflux = xr.open_dataset('net_rad_planet_hemis_avg_1_xCO2_1901_1961_QFLUX.nc').rad	
	rh_2_qflux = xr.open_dataset('net_rad_planet_hemis_avg_2_xCO2_1901_1961_QFLUX.nc').rad
	rh_3_qflux = xr.open_dataset('net_rad_planet_hemis_avg_3_xCO2_1901_1961_QFLUX.nc').rad
	rh_4_qflux = xr.open_dataset('net_rad_planet_hemis_avg_4_xCO2_1901_1961_QFLUX.nc').rad
	
	area = xr.open_dataset('temp_1_surface.nc').axyp.mean(dim=['month','year'])
	
	os.chdir("/nfs3m/archive/sfa_cache09/users/g00/imitevsk/E2.1_CO2_runs/pytropd/plots_figures")
			
	''' taking yearly and global weighted averages of temperature and net radiation of planet '''

	temp_05 = (t_05.temp * t_05.axyp / t_05.axyp.mean()).mean('month').mean(dim=['lat','lon'])
	temp_1 = (t_1.temp * t_1.axyp / t_1.axyp.mean()).mean('month').mean(dim=['lat','lon'])
	temp_15 = (t_15.temp * t_15.axyp / t_15.axyp.mean()).mean('month').mean(dim=['lat','lon'])
	temp_2 = (t_2.temp * t_2.axyp / t_2.axyp.mean()).mean('month').mean(dim=['lat','lon'])
	temp_3 = (t_3.temp * t_3.axyp / t_3.axyp.mean()).mean('month').mean(dim=['lat','lon'])
	temp_4 = (t_4.temp * t_4.axyp / t_4.axyp.mean()).mean('month').mean(dim=['lat','lon'])
	temp_5 = (t_5.temp * t_5.axyp / t_5.axyp.mean()).mean('month').mean(dim=['lat','lon'])
	temp_6 = (t_6.temp * t_6.axyp / t_6.axyp.mean()).mean('month').mean(dim=['lat','lon'])
	temp_7 = (t_7.temp * t_7.axyp / t_7.axyp.mean()).mean('month').mean(dim=['lat','lon'])
	temp_8 = (t_8.temp * t_8.axyp / t_8.axyp.mean()).mean('month').mean(dim=['lat','lon'])
		
	temp_1_qflux = (t_1_qflux.temp * t_1_qflux.axyp / t_1_qflux.axyp.mean()).mean('month').mean(dim=['lat','lon'])
	temp_2_qflux = (t_2_qflux.temp * t_2_qflux.axyp / t_2_qflux.axyp.mean()).mean('month').mean(dim=['lat','lon'])
	temp_3_qflux = (t_3_qflux.temp * t_3_qflux.axyp / t_3_qflux.axyp.mean()).mean('month').mean(dim=['lat','lon'])
	temp_4_qflux = (t_4_qflux.temp * t_4_qflux.axyp / t_4_qflux.axyp.mean()).mean('month').mean(dim=['lat','lon'])
	
	rad_h_05 = rh_05.sel(hemis='global').groupby('time.year').mean('time')	
	rad_h_1 = rh_1.sel(hemis='global').groupby('time.year').mean('time')
	rad_h_15 = rh_15.sel(hemis='global').groupby('time.year').mean('time')
	rad_h_2 = rh_2.sel(hemis='global').groupby('time.year').mean('time')
	rad_h_3 = rh_3.sel(hemis='global').groupby('time.year').mean('time')
	rad_h_4 = rh_4.sel(hemis='global').groupby('time.year').mean('time')
	rad_h_5 = rh_5.sel(hemis='global').groupby('time.year').mean('time')
	rad_h_6 = rh_6.sel(hemis='global').groupby('time.year').mean('time')
	rad_h_7 = rh_7.sel(hemis='global').groupby('time.year').mean('time')
	rad_h_8 = rh_8.sel(hemis='global').groupby('time.year').mean('time')	
	
	rad_h_1_qflux = rh_1_qflux.sel(hemis='global').groupby('time.year').mean('time')
	rad_h_2_qflux = rh_2_qflux.sel(hemis='global').groupby('time.year').mean('time')
	rad_h_3_qflux = rh_3_qflux.sel(hemis='global').groupby('time.year').mean('time')
	rad_h_4_qflux = rh_4_qflux.sel(hemis='global').groupby('time.year').mean('time')
		
	''' 149 years '''

	X = np.zeros((149, 9))
	Y = np.zeros((149, 9))

	X[:,0] = np.delete(np.array(temp_05 - temp_1), 42)
	X[:,1] = np.delete(np.array(temp_15 - temp_1), 42)
	X[:,2] = np.delete(np.array(temp_2 - temp_1), 42)
	X[:,3] = np.delete(np.array(temp_3 - temp_1), 42)
	X[:,4] = np.delete(np.array(temp_4 - temp_1), 42)
	X[:,5] = np.delete(np.array(temp_5 - temp_1), 42)
	X[:,6] = np.delete(np.array(temp_6 - temp_1), 42)
	X[:,7] = np.delete(np.array(temp_7 - temp_1), 42)
	X[:,8] = np.delete(np.array(temp_8 - temp_1), 42)
	
	Y[:,0] = np.delete(np.array(rad_h_05 - rad_h_1), 42)	
	Y[:,1] = np.delete(np.array(rad_h_15 - rad_h_1), 42)
	Y[:,2] = np.delete(np.array(rad_h_2 - rad_h_1), 42)
	Y[:,3] = np.delete(np.array(rad_h_3 - rad_h_1), 42)
	Y[:,4] = np.delete(np.array(rad_h_4 - rad_h_1), 42)
	Y[:,5] = np.delete(np.array(rad_h_5 - rad_h_1), 42)
	Y[:,6] = np.delete(np.array(rad_h_6 - rad_h_1), 42)
	Y[:,7] = np.delete(np.array(rad_h_7 - rad_h_1), 42)
	Y[:,8] = np.delete(np.array(rad_h_8 - rad_h_1), 42)	
	
	X_qflux = np.zeros((60, 3))
	Y_qflux = np.zeros((60, 3))
	
	X_qflux[:,0] = np.array(temp_2_qflux - temp_1_qflux)
	X_qflux[:,1] = np.array(temp_3_qflux - temp_1_qflux)
	X_qflux[:,2] = np.array(temp_4_qflux - temp_1_qflux)
	
	Y_qflux[:,0] = np.array(rad_h_2_qflux - rad_h_1_qflux)	
	Y_qflux[:,1] = np.array(rad_h_3_qflux - rad_h_1_qflux)
	Y_qflux[:,2] = np.array(rad_h_4_qflux - rad_h_1_qflux)
	
	title = ['0.5xCO$_2$', '1.5xCO$_2$', '2xCO$_2$', '3xCO$_2$', '4xCO$_2$', '5xCO$_2$', '6xCO$_2$', '7xCO$_2$', '8xCO$_2$']
	title_qflux = ['2xCO$_2$', '3xCO$_2$', '4xCO$_2$']
	colors = ['blue','green','orange','lime','red','magenta','deepskyblue','brown','navy']	
	
	fig = plt.figure()#figsize=[6,6])
	fig.set_figwidth(fig.get_figwidth()*1.8)
	#fig.set_figheight(fig.get_figheight())

	axes = fig.add_subplot(1,2,1)
	plt.axvline(x=0, color = 'black', linestyle='-', label = None)		
	plt.axhline(y=0, color = 'black', linestyle='-', label = None)
	for i in arange(1,10,1):	
		x = X[:,i-1]
		y = Y[:,i-1]
		c = scipy.stats.linregress(x,y)
		f = lambda x: c.intercept + c.slope * x
		axes.scatter(x, y, color = colors[i-1], label = None, s = 40, alpha=0.3)
		#label = title[i-1] + ', y = ' + str(np.around(c.intercept, 2)) + str(np.around(c.slope, 2)) + 'x, ' + 'x = ' + str(np.around(-c.intercept/c.slope,2))	
		label = title[i-1]	
		axes.plot([0, -c.intercept/c.slope], [c.intercept, 0], linewidth = 3, label = str(label), color = colors[i-1])
		axes.scatter([0], [c.intercept], color = 'black', s=50)
		axes.scatter([-c.intercept/c.slope], [0], color = 'black', s=50) 
	axes.set_title('A) Fully coupled')
	axes.set_ylabel('$\Delta R$ (W/m$^2$)')
	axes.set_xlabel('$\Delta T$ (K)')
	axes.legend(loc=0, ncol = 2)
	axes.set_xlim([-3,10])
	axes.set_ylim([-4,12])
	
	axes = fig.add_subplot(1,2,2)
	plt.axvline(x=0, color = 'black', linestyle='-', label = None)		
	plt.axhline(y=0, color = 'black', linestyle='-', label = None)
	for i in arange(1,4,1):
		x = X_qflux[:,i-1]
		y = Y_qflux[:,i-1]
		c = scipy.stats.linregress(x,y)
		f = lambda x: c.intercept + c.slope * x
		axes.scatter(x, y, color = colors[i-1+2], label = None, s = 40, alpha=0.3)
		#label = title_qflux[i-1] + ', y = ' + str(np.around(c.intercept, 2)) + str(np.around(c.slope, 2)) + 'x, ' + 'x = ' + str(np.around(-c.intercept/c.slope,2))	
		label = title_qflux[i-1]
		axes.plot([0, -c.intercept/c.slope], [c.intercept, 0], linewidth = 3, label = str(label), color = colors[i-1+2])
		axes.scatter([0], [c.intercept], color = 'black', s=50)
		axes.scatter([-c.intercept/c.slope], [0], color = 'black', s=50) 
	axes.set_title('B) Slab ocean', usetex=True)
	#axes.set_ylabel('$\Delta R$ (W/m$^2$)')
	axes.set_xlabel('$\Delta T$ (K)')
	axes.legend(loc=0)
	axes.set_xlim([-3,10])
	axes.set_ylim([-4,12])
	
	fig.tight_layout()
	plt.savefig('gregory_plot_supplemental_figure.pdf')
	plt.show()

gregory_plot()

