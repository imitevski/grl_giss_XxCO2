from pylab import * 
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
from matplotlib.ticker import StrMethodFormatter, NullFormatter
#from IPython.display import set_matplotlib_formats
#set_matplotlib_formats('png', 'pdf')
import glob
import os

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

os.chdir("/nfs3m/archive/sfa_cache09/users/g00/imitevsk/E2.1_CO2_runs/pytropd/plots_figures")

def supplemental_fig_ice_extent():
	x = np.array([0.5, 1, 1.5, 2, 3, 4, 5, 6, 7, 8])
	x = np.log(x) * 5.5
	
	#### Ice calculations ### 
	y_ice = np.array([	
		(t_05.ice * t_05.axyp)[-50:,:,:,:].mean() / t_05.axyp[-50:,:,:,:].mean(),
		(t_1.ice * t_1.axyp)[-50:,:,:,:].mean() / t_1.axyp[-50:,:,:,:].mean(),
		(t_15.ice * t_15.axyp)[-50:,:,:,:].mean() / t_15.axyp[-50:,:,:,:].mean(),
	   	(t_2.ice * t_2.axyp)[-50:,:,:,:].mean() / t_2.axyp[-50:,:,:,:].mean(),
		(t_3.ice * t_3.axyp)[-50:,:,:,:].mean() / t_3.axyp[-50:,:,:,:].mean(),
		(t_4.ice * t_4.axyp)[-50:,:,:,:].mean() / t_4.axyp[-50:,:,:,:].mean(),
		(t_5.ice * t_5.axyp)[-50:,:,:,:].mean() / t_5.axyp[-50:,:,:,:].mean(),
		(t_6.ice * t_6.axyp)[-50:,:,:,:].mean() / t_6.axyp[-50:,:,:,:].mean(),
		(t_7.ice * t_7.axyp)[-50:,:,:,:].mean() / t_7.axyp[-50:,:,:,:].mean(),
		(t_8.ice * t_8.axyp)[-50:,:,:,:].mean() / t_8.axyp[-50:,:,:,:].mean()])
	
	sie_05 = t_05.ice.where(t_05.ice > 15) / 100
	sie_1 = t_1.ice.where(t_1.ice > 15) / 100
	sie_15 = t_15.ice.where(t_15.ice > 15) / 100
	sie_2 = t_2.ice.where(t_2.ice > 15) / 100
	sie_3 = t_3.ice.where(t_3.ice > 15) / 100
	sie_4 = t_4.ice.where(t_4.ice > 15) / 100
	sie_5 = t_5.ice.where(t_5.ice > 15) / 100
	sie_6 = t_6.ice.where(t_6.ice > 15) / 100
	sie_7 = t_7.ice.where(t_7.ice > 15) / 100
	sie_8 = t_8.ice.where(t_8.ice > 15) / 100

	y_sie_nh = np.array([	
		(sie_05 * t_05.axyp).sel(year=slice(1950,2000), lat=slice(0,90)).sum(dim=['lat','lon']).mean('month').mean('year') / 1e12 ,
		(sie_1 * t_1.axyp).sel(year=slice(1950,2000), lat=slice(0,90)).sum(dim=['lat','lon']).mean('month').mean('year') / 1e12 ,
		(sie_15 * t_15.axyp).sel(year=slice(1950,2000), lat=slice(0,90)).sum(dim=['lat','lon']).mean('month').mean('year') / 1e12 ,
		(sie_2 * t_2.axyp).sel(year=slice(1950,2000), lat=slice(0,90)).sum(dim=['lat','lon']).mean('month').mean('year') / 1e12 ,
		(sie_3 * t_3.axyp).sel(year=slice(1950,2000), lat=slice(0,90)).sum(dim=['lat','lon']).mean('month').mean('year') / 1e12 ,
		(sie_4 * t_4.axyp).sel(year=slice(1950,2000), lat=slice(0,90)).sum(dim=['lat','lon']).mean('month').mean('year') / 1e12 ,
		(sie_5 * t_5.axyp).sel(year=slice(1950,2000), lat=slice(0,90)).sum(dim=['lat','lon']).mean('month').mean('year') / 1e12 ,
		(sie_6 * t_6.axyp).sel(year=slice(1950,2000), lat=slice(0,90)).sum(dim=['lat','lon']).mean('month').mean('year') / 1e12 ,
		(sie_7 * t_7.axyp).sel(year=slice(1950,2000), lat=slice(0,90)).sum(dim=['lat','lon']).mean('month').mean('year') / 1e12 ,
		(sie_8 * t_8.axyp).sel(year=slice(1950,2000), lat=slice(0,90)).sum(dim=['lat','lon']).mean('month').mean('year') / 1e12 ])
		
	y_sie_nh_error = np.array([	
		((sie_05 * t_05.axyp).sel(year=slice(1950,2000), lat=slice(0,90)).sum(dim=['lat','lon']).mean('month') / 1e12).std() ,
		((sie_1 * t_1.axyp).sel(year=slice(1950,2000), lat=slice(0,90)).sum(dim=['lat','lon']).mean('month') / 1e12).std() ,
		((sie_15 * t_15.axyp).sel(year=slice(1950,2000), lat=slice(0,90)).sum(dim=['lat','lon']).mean('month') / 1e12).std() ,
		((sie_2 * t_2.axyp).sel(year=slice(1950,2000), lat=slice(0,90)).sum(dim=['lat','lon']).mean('month') / 1e12).std() ,
		((sie_3 * t_3.axyp).sel(year=slice(1950,2000), lat=slice(0,90)).sum(dim=['lat','lon']).mean('month') / 1e12).std() ,
		((sie_4 * t_4.axyp).sel(year=slice(1950,2000), lat=slice(0,90)).sum(dim=['lat','lon']).mean('month') / 1e12).std() ,
		((sie_5 * t_5.axyp).sel(year=slice(1950,2000), lat=slice(0,90)).sum(dim=['lat','lon']).mean('month') / 1e12).std() ,
		((sie_6 * t_6.axyp).sel(year=slice(1950,2000), lat=slice(0,90)).sum(dim=['lat','lon']).mean('month') / 1e12).std() ,
		((sie_7 * t_7.axyp).sel(year=slice(1950,2000), lat=slice(0,90)).sum(dim=['lat','lon']).mean('month') / 1e12).std() ,
		((sie_8 * t_8.axyp).sel(year=slice(1950,2000), lat=slice(0,90)).sum(dim=['lat','lon']).mean('month') / 1e12).std() ])

	y_sie_sh = np.array([	
		(sie_05 * t_05.axyp).sel(year=slice(1950,2000), lat=slice(-90,0)).sum(dim=['lat','lon']).mean('month').mean('year') / 1e12 ,
		(sie_1 * t_1.axyp).sel(year=slice(1950,2000), lat=slice(-90,0)).sum(dim=['lat','lon']).mean('month').mean('year') / 1e12 ,
		(sie_15 * t_15.axyp).sel(year=slice(1950,2000), lat=slice(-90,0)).sum(dim=['lat','lon']).mean('month').mean('year') / 1e12 ,
		(sie_2 * t_2.axyp).sel(year=slice(1950,2000), lat=slice(-90,0)).sum(dim=['lat','lon']).mean('month').mean('year') / 1e12 ,
		(sie_3 * t_3.axyp).sel(year=slice(1950,2000), lat=slice(-90,0)).sum(dim=['lat','lon']).mean('month').mean('year') / 1e12 ,
		(sie_4 * t_4.axyp).sel(year=slice(1950,2000), lat=slice(-90,0)).sum(dim=['lat','lon']).mean('month').mean('year') / 1e12 ,
		(sie_5 * t_5.axyp).sel(year=slice(1950,2000), lat=slice(-90,0)).sum(dim=['lat','lon']).mean('month').mean('year') / 1e12 ,
		(sie_6 * t_6.axyp).sel(year=slice(1950,2000), lat=slice(-90,0)).sum(dim=['lat','lon']).mean('month').mean('year') / 1e12 ,
		(sie_7 * t_7.axyp).sel(year=slice(1950,2000), lat=slice(-90,0)).sum(dim=['lat','lon']).mean('month').mean('year') / 1e12 ,
		(sie_8 * t_8.axyp).sel(year=slice(1950,2000), lat=slice(-90,0)).sum(dim=['lat','lon']).mean('month').mean('year') / 1e12 ])

	y_sie_sh_error = np.array([	
		((sie_05 * t_05.axyp).sel(year=slice(1950,2000), lat=slice(-90,0)).sum(dim=['lat','lon']).mean('month') / 1e12).std() ,
		((sie_1 * t_1.axyp).sel(year=slice(1950,2000), lat=slice(-90,0)).sum(dim=['lat','lon']).mean('month') / 1e12).std() ,
		((sie_15 * t_15.axyp).sel(year=slice(1950,2000), lat=slice(-90,0)).sum(dim=['lat','lon']).mean('month') / 1e12).std() ,
		((sie_2 * t_2.axyp).sel(year=slice(1950,2000), lat=slice(-90,0)).sum(dim=['lat','lon']).mean('month') / 1e12).std() ,
		((sie_3 * t_3.axyp).sel(year=slice(1950,2000), lat=slice(-90,0)).sum(dim=['lat','lon']).mean('month') / 1e12).std() ,
		((sie_4 * t_4.axyp).sel(year=slice(1950,2000), lat=slice(-90,0)).sum(dim=['lat','lon']).mean('month') / 1e12).std() ,
		((sie_5 * t_5.axyp).sel(year=slice(1950,2000), lat=slice(-90,0)).sum(dim=['lat','lon']).mean('month') / 1e12).std() ,
		((sie_6 * t_6.axyp).sel(year=slice(1950,2000), lat=slice(-90,0)).sum(dim=['lat','lon']).mean('month') / 1e12).std() ,
		((sie_7 * t_7.axyp).sel(year=slice(1950,2000), lat=slice(-90,0)).sum(dim=['lat','lon']).mean('month') / 1e12).std() ,
		((sie_8 * t_8.axyp).sel(year=slice(1950,2000), lat=slice(-90,0)).sum(dim=['lat','lon']).mean('month') / 1e12).std() ])		
		
	fig = plt.figure()
	fig.set_figwidth(fig.get_figwidth() * 1)
	fig.set_figheight(fig.get_figheight() * 1)
			
	axes = fig.add_subplot(1,1,1)
	axes.errorbar(x, y_sie_nh, y_sie_nh_error, linewidth = 3, marker='o', markersize=12, capsize = 8, color = 'red', label='NH')
	axes.errorbar(x, y_sie_sh, y_sie_sh_error, linewidth = 3, marker='o', markersize=12, capsize = 8, color = 'blue', label='SH')	
	plt.xlabel('Radiative Forcing (W/m$^2$)', fontsize=20)
	plt.ylabel('10$^6$ km$^2$', fontsize=18)
	plt.title("Sea-ice extent", fontsize = 20)
	axes.legend(loc = 0)
	plt.xticks([-4,0,4,8,12], fontsize = 15)
	plt.yticks([2.5,7.5,12.5,17.5], fontsize = 15)
		
	plt.tight_layout()
	plt.savefig('supplemental_fig_ice_extent.pdf')
	plt.show()
	
supplemental_fig_ice_extent() 
