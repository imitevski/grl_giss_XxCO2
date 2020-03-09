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
import scipy 
import scipy.stats


	
def amoc_paper_figure(l_s, l_e, d_s, d_e, mov_avg):
		
	'''
	Plots AMOC in simulations from A) abrupt XxCO2 runs in GISS E2.1 and B) abrupt 4xCO2 CMIP6 models
 
	Arguments:
	- l_s: latitude start 
	- l_e: latitude end
	- d_s: depth start 
	- d_e: depth end
	- mov_avg: moving average in years 
	'''

	os.chdir("/nfs3m/archive/sfa_cache09/users/g00/imitevsk/E2.1_CO2_runs/pytropd/output")

	psi_05 = xr.open_dataset('psi_sf_Atl_0.5.nc').sf_Atl
	psi_1 = xr.open_dataset('psi_sf_Atl_1.nc').sf_Atl
	psi_15 = xr.open_dataset('psi_sf_Atl_1.5.nc').sf_Atl
	psi_2 = xr.open_dataset('psi_sf_Atl_2.nc').sf_Atl
	psi_3 = xr.open_dataset('psi_sf_Atl_3.nc').sf_Atl
	psi_4 = xr.open_dataset('psi_sf_Atl_4.nc').sf_Atl
	psi_5 = xr.open_dataset('psi_sf_Atl_5.nc').sf_Atl
	psi_6 = xr.open_dataset('psi_sf_Atl_6.nc').sf_Atl
	psi_8 = xr.open_dataset('psi_sf_Atl_8.nc').sf_Atl

	c6_CanESM5_1 = np.array(xr.open_dataset('/gpfsm/dnb02/projects/p54/users/imitevsk/cmip6/abrupt-4xCO2/CanESM5/msftmz_Omon_CanESM5_abrupt-4xCO2_r1i1p1f1_gn_185001-200012.nc')\
			.msftmz.sel(lev=slice(d_s,d_e),lat=slice(l_s,l_e)).max(dim=['lev','lat']).groupby('time.year').mean('time').rolling(year = mov_avg, center=True).mean()[:150,0])
	c6_CanESM5_2 = np.array(xr.open_dataset('/gpfsm/dnb02/projects/p54/users/imitevsk/cmip6/abrupt-4xCO2/CanESM5/msftmz_Omon_CanESM5_abrupt-4xCO2_r1i1p2f1_gn_185001-200012.nc')\
			.msftmz.sel(lev=slice(d_s,d_e),lat=slice(l_s,l_e)).max(dim=['lev','lat']).groupby('time.year').mean('time').rolling(year = mov_avg, center=True).mean()[:150,0])
	c6_CESM2 = np.array(xr.open_dataset('/gpfsm/dnb02/projects/p54/users/imitevsk/cmip6/abrupt-4xCO2/CESM2/merged_msftmz_Omon_CESM2_abrupt-4xCO2_r1i1p1f1_gn.nc')\
			.msftmz.sel(lev=slice(d_s*100,d_e*100),lat=slice(l_s,l_e)).max(dim=['lev','lat']).groupby('time.year').mean('time').rolling(year = mov_avg, center=True).mean()[:150,0])
	c6_CESM2_WACCM = np.array(xr.open_dataset('/gpfsm/dnb02/projects/p54/users/imitevsk/cmip6/abrupt-4xCO2/CESM2-WACCM/merged_msftmz_Omon_CESM2-WACCM_abrupt-4xCO2_r1i1p1f1_gn.nc')\
			.msftmz.sel(lev=slice(d_s*100,d_e*100),lat=slice(l_s,l_e)).max(dim=['lev','lat']).groupby('time.year').mean('time').rolling(year = mov_avg, center=True).mean()[:150,0])
	c6_E3SM_1_0 = np.array(xr.open_dataset('/gpfsm/dnb02/projects/p54/users/imitevsk/cmip6/abrupt-4xCO2/E3SM-1-0/merged_msftmz_Omon_E3SM-1-0_abrupt-4xCO2_r1i1p1f1_gr.nc')\
			.msftmz.sel(lev=slice(d_s,d_e),lat=slice(l_s,l_e)).max(dim=['lev','lat']).groupby('time.year').mean('time').rolling(year = mov_avg, center=True).mean()[:150,0])
	c6_EC_Earth3 = np.array(xr.open_dataset('/gpfsm/dnb02/projects/p54/users/imitevsk/cmip6/abrupt-4xCO2/EC-Earth3/merged_msftyz_Omon_EC-Earth3_1pctCO2_r3i1p1f1_gn.nc')\
			.msftyz.sel(lev=slice(d_s,d_e),rlat=slice(l_s,l_e)).max(dim=['lev','rlat']).groupby('time.year').mean('time').rolling(year = mov_avg, center=True).mean()[:150,1]) # difference
	c6_EC_Earth3_Veg = np.array(xr.open_dataset('/gpfsm/dnb02/projects/p54/users/imitevsk/cmip6/abrupt-4xCO2/EC-Earth3-Veg/merged_msftyz_Omon_EC-Earth3-Veg_1pctCO2_r1i1p1f1_gn.nc')\
			.msftyz.sel(lev=slice(d_s,d_e),rlat=slice(l_s,l_e)).max(dim=['lev','rlat']).groupby('time.year').mean('time').rolling(year = mov_avg, center=True).mean()[:150,1]) # difference
	c6_FGOALS_f3_L_1 = np.array(xr.open_dataset('/gpfsm/dnb02/projects/p54/users/imitevsk/cmip6/abrupt-4xCO2/FGOALS-f3-L/msftmz_Omon_FGOALS-f3-L_abrupt-4xCO2_r1i1p1f1_gn_185001-200912.nc')\
			 .msftmz.sel(lev=slice(d_s,d_e),lat=slice(l_s,l_e)).max(dim=['lev','lat']).groupby('time.year').mean('time').rolling(year = mov_avg, center=True).mean()[:150,0])
	c6_FGOALS_f3_L_2 = np.array(xr.open_dataset('/gpfsm/dnb02/projects/p54/users/imitevsk/cmip6/abrupt-4xCO2/FGOALS-f3-L/msftmz_Omon_FGOALS-f3-L_abrupt-4xCO2_r2i1p1f1_gn_185001-200912.nc')\
			 .msftmz.sel(lev=slice(d_s,d_e),lat=slice(l_s,l_e)).max(dim=['lev','lat']).groupby('time.year').mean('time').rolling(year = mov_avg, center=True).mean()[:150,0])
	c6_FGOALS_f3_L_3 = np.array(xr.open_dataset('/gpfsm/dnb02/projects/p54/users/imitevsk/cmip6/abrupt-4xCO2/FGOALS-f3-L/msftmz_Omon_FGOALS-f3-L_abrupt-4xCO2_r3i1p1f1_gn_185001-200912.nc')\
			 .msftmz.sel(lev=slice(d_s,d_e),lat=slice(l_s,l_e)).max(dim=['lev','lat']).groupby('time.year').mean('time').rolling(year = mov_avg, center=True).mean()[:150,0])
	#c6_GISS_E2_1_G_1 = np.array(xr.open_dataset('/gpfsm/dnb02/projects/p54/users/imitevsk/cmip6/abrupt-4xCO2/GISS-E2-1-G/merged_msftmz_Omon_GISS-E2-1-G_abrupt-4xCO2_r102i1p1f1_gn.nc')\
	#		.msftmz.sel(lev=slice(d_s,d_e),lat=slice(l_s,l_e)).max(dim=['lev','lat']).groupby('time.year').mean('time').rolling(year = mov_avg, center=True).mean()[:150,0])
	c6_GISS_E2_1_G_2 = np.array(xr.open_dataset('/gpfsm/dnb02/projects/p54/users/imitevsk/cmip6/abrupt-4xCO2/GISS-E2-1-G/merged_msftmz_Omon_GISS-E2-1-G_abrupt-4xCO2_r1i1p3f1_gn.nc')\
			.msftmz.sel(lev=slice(d_s,d_e),lat=slice(l_s,l_e)).max(dim=['lev','lat']).groupby('time.year').mean('time').rolling(year = mov_avg, center=True).mean()[:150,0])	
	c6_INM_CM4_8 = np.array(xr.open_dataset('/gpfsm/dnb02/projects/p54/users/imitevsk/cmip6/abrupt-4xCO2/INM-CM4-8/msftmz_Omon_INM-CM4-8_abrupt-4xCO2_r1i1p1f1_gr1_185001-199912.nc')\
			.msftmz.sel(lev=slice(d_s,d_e),lat=slice(l_s,l_e)).max(dim=['lev','lat']).groupby('time.year').mean('time').rolling(year = mov_avg, center=True).mean()[:150,0])
	c6_INM_CM5_0 = np.array(xr.open_dataset('/gpfsm/dnb02/projects/p54/users/imitevsk/cmip6/abrupt-4xCO2/INM-CM5-0/msftmz_Omon_INM-CM5-0_abrupt-4xCO2_r1i1p1f1_gr1_185001-199912.nc')\
			.msftmz.sel(lev=slice(d_s,d_e),lat=slice(l_s,l_e)).max(dim=['lev','lat']).groupby('time.year').mean('time').rolling(year = mov_avg, center=True).mean()[:150,0])
	c6_MPI_ESM1_2_HR = np.array(xr.open_dataset('/gpfsm/dnb02/projects/p54/users/imitevsk/cmip6/abrupt-4xCO2/MPI-ESM1-2-HR/merged_msftmz_Omon_MPI-ESM1-2-HR_abrupt-4xCO2_r1i1p1f1_gn.nc')\
			.msftmz.sel(lev=slice(d_s,d_e),lat=slice(l_s,l_e)).max(dim=['lev','lat']).groupby('time.year').mean('time').rolling(year = mov_avg, center=True).mean()[:150,1]) # difference
	c6_MRI_ESM2_0_1 = np.array(xr.open_dataset('/gpfsm/dnb02/projects/p54/users/imitevsk/cmip6/abrupt-4xCO2/MRI-ESM2-0/msftmz_Omon_MRI-ESM2-0_abrupt-4xCO2_r1i1p1f1_gr2z_185001-200012.nc')\
			.msftmz.sel(lev=slice(d_s,d_e),lat=slice(l_s,l_e)).max(dim=['lev','lat']).groupby('time.year').mean('time').rolling(year = mov_avg, center=True).mean()[:150,0]) 
	c6_MRI_ESM2_0_2 = np.array(xr.open_dataset('/gpfsm/dnb02/projects/p54/users/imitevsk/cmip6/abrupt-4xCO2/MRI-ESM2-0/msftmz_Omon_MRI-ESM2-0_abrupt-4xCO2_r4i1p1f1_gr2z_185001-200012.nc')\
			.msftmz.sel(lev=slice(d_s,d_e),lat=slice(l_s,l_e)).max(dim=['lev','lat']).groupby('time.year').mean('time').rolling(year = mov_avg, center=True).mean()[:150,0])	
	c6_MRI_ESM2_0_3 = np.array(xr.open_dataset('/gpfsm/dnb02/projects/p54/users/imitevsk/cmip6/abrupt-4xCO2/MRI-ESM2-0/msftmz_Omon_MRI-ESM2-0_abrupt-4xCO2_r7i1p1f1_gr2z_185001-200012.nc')\
			.msftmz.sel(lev=slice(d_s,d_e),lat=slice(l_s,l_e)).max(dim=['lev','lat']).groupby('time.year').mean('time').rolling(year = mov_avg, center=True).mean()[:150,0])
	c6_MRI_ESM2_0_4 = np.array(xr.open_dataset('/gpfsm/dnb02/projects/p54/users/imitevsk/cmip6/abrupt-4xCO2/MRI-ESM2-0/msftmz_Omon_MRI-ESM2-0_abrupt-4xCO2_r10i1p1f1_gr2z_185001-200012.nc')\
			.msftmz.sel(lev=slice(d_s,d_e),lat=slice(l_s,l_e)).max(dim=['lev','lat']).groupby('time.year').mean('time').rolling(year = mov_avg, center=True).mean()[:150,0])
	c6_MRI_ESM2_0_5 = np.array(xr.open_dataset('/gpfsm/dnb02/projects/p54/users/imitevsk/cmip6/abrupt-4xCO2/MRI-ESM2-0/msftmz_Omon_MRI-ESM2-0_abrupt-4xCO2_r13i1p1f1_gr2z_185001-200012.nc')\
			.msftmz.sel(lev=slice(d_s,d_e),lat=slice(l_s,l_e)).max(dim=['lev','lat']).groupby('time.year').mean('time').rolling(year = mov_avg, center=True).mean()[:,0])
	c6_NorCPM1 = np.array(xr.open_dataset('/gpfsm/dnb02/projects/p54/users/imitevsk/cmip6/abrupt-4xCO2/NorCPM1/msftmz_Omon_NorCPM1_abrupt-4xCO2_r1i1p1f1_grz_000101-008012.nc')\
			.msftmz.sel(lev=slice(d_s,d_e),lat=slice(l_s,l_e)).max(dim=['lev','lat']).groupby('time.year').mean('time').rolling(year = mov_avg, center=True).mean()[:150,0])
	c6_NorESM2_LM = np.array(xr.open_dataset('/gpfsm/dnb02/projects/p54/users/imitevsk/cmip6/abrupt-4xCO2/NorESM2-LM/merged_msftmz_Omon_NorESM2-LM_abrupt-4xCO2_r1i1p1f1_grz.nc')\
			.msftmz.sel(lev=slice(d_s,d_e),lat=slice(l_s,l_e)).max(dim=['lev','lat']).groupby('time.year').mean('time').rolling(year = mov_avg, center=True).mean()[:150,0])
	c6_SAM0_UNICON = np.array(xr.open_dataset('/gpfsm/dnb02/projects/p54/users/imitevsk/cmip6/abrupt-4xCO2/SAM0-UNICON/merged_msftmz_Omon_SAM0-UNICON_abrupt-4xCO2_r1i1p1f1_gn.nc')\
			.msftmz.sel(lev=slice(d_s,d_e),lat=slice(l_s,l_e)).max(dim=['lev','lat']).groupby('time.year').mean('time').rolling(year = mov_avg, center=True).mean()[:150,0])
	
	os.chdir("/nfs3m/archive/sfa_cache09/users/g00/imitevsk/E2.1_CO2_runs/pytropd/plots_figures")

	fig = plt.figure()
	fig.set_figwidth(fig.get_figwidth() * 2)
	#fig.set_figheight(fig.get_figheight() * 3)
	
	axes = fig.add_subplot(1,2,1)
	psi_05.sel(lato2=slice(l_s,l_e), zoce=slice(d_s,d_e)).max(dim=['zoce','lato2']).groupby('time.year').mean('time').rolling(year = mov_avg, center = True).mean().plot(label = '0.5xCO$_2$', color = 'blue')
	psi_1.sel(lato2=slice(l_s,l_e), zoce=slice(d_s,d_e)).max(dim=['zoce','lato2']).groupby('time.year').mean('time').rolling(year = mov_avg, center = True).mean().plot(label = '1xCO$_2$', color = 'black')
	psi_15.sel(lato2=slice(l_s,l_e), zoce=slice(d_s,d_e)).max(dim=['zoce','lato2']).groupby('time.year').mean('time').rolling(year = mov_avg, center = True).mean().plot(label = '1.5xCO$_2$', color = 'green')
	psi_2.sel(lato2=slice(l_s,l_e), zoce=slice(d_s,d_e)).max(dim=['zoce','lato2']).groupby('time.year').mean('time').rolling(year = mov_avg, center = True).mean().plot(label = '2xCO$_2$', color = 'orange')
	psi_3.sel(lato2=slice(l_s,l_e), zoce=slice(d_s,d_e)).max(dim=['zoce','lato2']).groupby('time.year').mean('time').rolling(year = mov_avg, center = True).mean().plot(label = '3xCO$_2$', color = 'lime')	
	psi_5.sel(lato2=slice(l_s,l_e), zoce=slice(d_s,d_e)).max(dim=['zoce','lato2']).groupby('time.year').mean('time').rolling(year = mov_avg, center = True).mean().plot(label = '5xCO$_2$', color = 'magenta')
	psi_6.sel(lato2=slice(l_s,l_e), zoce=slice(d_s,d_e)).max(dim=['zoce','lato2']).groupby('time.year').mean('time').rolling(year = mov_avg, center = True).mean().plot(label = '6xCO$_2$', color = 'deepskyblue')
	psi_8.sel(lato2=slice(l_s,l_e), zoce=slice(d_s,d_e)).max(dim=['zoce','lato2']).groupby('time.year').mean('time').rolling(year = mov_avg, center = True).mean().plot(label = '8xCO$_2$', color = 'brown')	
	psi_4.sel(lato2=slice(l_s,l_e), zoce=slice(d_s,d_e)).max(dim=['zoce','lato2']).groupby('time.year').mean('time').rolling(year = mov_avg, center = True).mean().plot(label = '4xCO$_2$', color = 'red', linewidth = 3)
	plt.title('A) Abrupt XxCO$_2$ from GISS Model E2.1', fontsize = 17, usetex=True)
	plt.xlabel('Years', fontsize = 15, usetex=True)
	plt.xlim([1850,2000])
	plt.ylim([0,37])
	plt.xticks(np.arange(1850,2025,25), np.arange(0,175,25), fontsize = 15, usetex=True)
	plt.yticks([0, 5, 10, 15, 20, 25, 35], fontsize = 15, usetex=True)
	plt.ylabel('AMOC (Sv)', fontsize = 15, usetex=True)
	plt.rc('text', usetex=True)
	plt.legend(loc=0, ncol = 2)
	
	axes = fig.add_subplot(1,2,2)	
	plt.plot(c6_CanESM5_1, label = 'CanESM5', color = 'pink')# label = 'CanESM5, r1i1p1f1'
	plt.plot(c6_CanESM5_2, color = 'pink')#, label = 'CanESM5, r1i1p2f1')
	#plt.plot(c6_E3SM_1_0, label = 'E3SM-1-0, r1i1p1f1', color = 'lightskyblue') # it shows zero but it appears to have some dimension issue (10^3 difference)
	plt.plot(c6_EC_Earth3, label = 'EC\_Earth3', color = 'black')# label = 'EC_Earth3, r3i1p1f1'
	plt.plot(c6_EC_Earth3_Veg, label = 'EC\_Earth3\_Veg', color = 'black', linestyle = '--')# label = 'EC_Earth3_Veg, r1i1p1f1'
	plt.plot(c6_FGOALS_f3_L_1, color = 'gray', label = 'FGOALS-f3-L')# label = 'FGOALS-f3-L, r1i1p1f1'
	plt.plot(c6_FGOALS_f3_L_2, color = 'gray')#, label = 'FGOALS-f3-L, r2i1p1f1')
	plt.plot(c6_FGOALS_f3_L_3, color = 'gray')#, label = 'FGOALS-f3-L, r3i1p1f1')	
	plt.plot(c6_INM_CM5_0, label = 'INM\_CM5\_0', color = 'blue')# label = 'INM_CM5_0, r1i1p1f1'
	plt.plot(c6_INM_CM4_8, label = 'INM\_CM4\_8', color = 'blue', linestyle='--')# label = 'INM_CM4_8, r1i1p1f1'
	plt.plot(c6_MPI_ESM1_2_HR, label = 'MPI-ESM1-2-HR', color = 'cyan')# label = 'MPI-ESM1-2-HR, r1i1p1f1'
	plt.plot(c6_MRI_ESM2_0_1, color = 'yellow', label = 'MRI-ESM2-0')# label = 'MRI-ESM2-0, r1i1p1f1'
	plt.plot(c6_MRI_ESM2_0_2, color = 'yellow')#, label = 'MRI-ESM2-0, r4i1p1f1')
	plt.plot(c6_MRI_ESM2_0_3, color = 'yellow')#, label = 'MRI-ESM2-0, r7i1p1f1')
	plt.plot(c6_MRI_ESM2_0_4, color = 'yellow')#, label = 'MRI-ESM2-0, r10i1p1f1')
	plt.plot(c6_MRI_ESM2_0_5, color = 'yellow')#, label = 'MRI-ESM2-0, r13i1p1f1')
	#plt.plot(c6_NorCPM1, label = 'NorCPM1, r1i1p1f1', color = 'lime', linestyle='--') it only runs for 80 years
	plt.plot(c6_NorESM2_LM, label = 'NorESM2-LM', color = 'lime')# label = 'NorESM2-LM, r1i1p1f1'
	plt.plot(c6_SAM0_UNICON, label = 'SAM0-UNICON, r1i1p1f1', color = 'chocolate')
	plt.plot(c6_CESM2, label = 'CESM2', color = 'purple')# label = 'CESM2, r1i1p1f1'
	plt.plot(c6_CESM2_WACCM, label = 'CESM2-WACCM', color = 'purple', linestyle='--')# label = 'CESM2-WACCM, r1i1p1f1'
	#plt.plot(c6_GISS_E2_1_G_1, label = 'GISS-E2-1-G, r102i1p1f1', color = 'red', linestyle = '--')
	plt.plot(c6_GISS_E2_1_G_2, label = 'GISS-E2-1-G', color = 'red', linewidth = 3)# label = 'GISS-E2-1-G, r1i1p3f1'
	plt.title('B) Abrupt 4xCO$_2$ runs from CMIP6 models', fontsize = 17, usetex=True)
	plt.xlabel('Years', fontsize = 15, usetex=True)
	plt.xticks([0,25,50,75,100,125,150], fontsize = 15, usetex=True)
	plt.xlim([0,150])
	plt.yticks([0,0.5e10, 1.0e10, 1.5e10, 2.0e10, 2.5e10, 3.0e10, 3.5e10], [0, 5, 10, 15, 20, 25, 30, 35], fontsize = 15, usetex=True)
	plt.ylim([0,3.7e10])
	plt.ylabel('AMOC (Sv)', fontsize = 15, usetex=True)
	plt.legend(loc=0, ncol = 2)#, fontsize = 4, ncol = 3)	
	
	
	
	fig.tight_layout()
	name_file = 'amoc_paper_figure_mov_avg_' + str(mov_avg) + '_years.pdf'
	plt.savefig(name_file)
	plt.show()

amoc_paper_figure(l_s = 30, l_e = 55, d_s = 800, d_e = 2000, mov_avg = 5)

