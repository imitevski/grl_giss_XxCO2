'''
This is similar to Figure 4 (figure_4_pd3.pdf) except here we use QFLUX data
'''

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


### load datasets ###

os.chdir("/nfs3m/archive/sfa_cache09/users/g00/imitevsk/E2.1_CO2_runs/pytropd/output")

pe_1 = xr.open_dataset('u_pe_1_all_QFLUX.nc').pe
pe_2 = xr.open_dataset('u_pe_2_all_QFLUX.nc').pe
pe_3 = xr.open_dataset('u_pe_3_all_QFLUX.nc').pe
pe_4 = xr.open_dataset('u_pe_4_all_QFLUX.nc').pe

t_1 = xr.open_dataset('temp_1_1901_1960_surface_QFLUX.nc')
t_2 = xr.open_dataset('temp_2_1901_1960_surface_QFLUX.nc')
t_3 = xr.open_dataset('temp_3_1901_1960_surface_QFLUX.nc')
t_4 = xr.open_dataset('temp_4_1901_1960_surface_QFLUX.nc')

u_1 = xr.open_dataset('u_psi_1_all_QFLUX.nc').psi
u_2 = xr.open_dataset('u_psi_2_all_QFLUX.nc').psi
u_3 = xr.open_dataset('u_psi_3_all_QFLUX.nc').psi
u_4 = xr.open_dataset('u_psi_4_all_QFLUX.nc').psi

psi_max = xr.open_dataset('psi_max_abs_500_mbar_all_xCO2_QFLUX_rev2.nc').psi_max

os.chdir("/nfs3m/archive/sfa_cache09/users/g00/imitevsk/E2.1_CO2_runs/pytropd/plots_figures")


def figure_4_pd3_QFLUX(t_s, t_e):
	"""
	Plots supplemental figure similar to Figure 4 in paper except with QFLUX runs data
	
	:Input:
	 - *t_s* (int) - start equilibrium year (range 1901-1960)
	 - *t_e* (int) - end equilibrium year (range 1901-1960)
	:Output:
	 - *figure* (.pdf) - figure with 4 panels showing A) Precipitation, B) P-E=0, C) HC width, D) HC strength
	"""

	#### A) Precipitation ###
	
	def y_prec(l_s, l_e):
		"""
		Calculates  precipitation at certain latitudinal band at equilibrium (last 10 years)
		
		:Input:		
		 - *l_s* (int) - starting latitude 
		 - *l_e* (int) - ending latitude
		:Output:
		 - *y* (ndarray) - precipitation averaged over last 50 years for 1,2,3,4xCO2
		"""
		y = np.array([
			(t_1.prec * t_1.axyp).sel(year = slice(t_s, t_e), lat=slice(l_s, l_e)).mean() / t_1.axyp.sel(year = slice(t_s, t_e), lat=slice(l_s, l_e)).mean(),
			(t_2.prec * t_2.axyp).sel(year = slice(t_s, t_e), lat=slice(l_s, l_e)).mean() / t_2.axyp.sel(year = slice(t_s, t_e), lat=slice(l_s, l_e)).mean(),
			(t_3.prec * t_3.axyp).sel(year = slice(t_s, t_e), lat=slice(l_s, l_e)).mean() / t_3.axyp.sel(year = slice(t_s, t_e), lat=slice(l_s, l_e)).mean(),
			(t_4.prec * t_4.axyp).sel(year = slice(t_s, t_e), lat=slice(l_s, l_e)).mean() / t_4.axyp.sel(year = slice(t_s, t_e), lat=slice(l_s, l_e)).mean() ])
		return y
			
	def e_prec(l_s, l_e):
		'''
		Calculates 1 \sigma of yearly variability of precipitation at certain latitudinal band at equilibrium (last 50 years) 
		
		:Input:		
		 - *l_s* (int) - starting latitude 
		 - *l_e* (int) - ending latitude
		:Output:
		 - *e* (ndarray) - 1 \sigma of yearly variations of precipitation for 1,2,3,4xCO2
		'''
		e = np.array([
			((t_1.prec * t_1.axyp).sel(year = slice(t_s, t_e), lat=slice(l_s, l_e)).mean(dim=['lat','lon','month']) / t_1.axyp.sel(year = slice(t_s, t_e), lat=slice(l_s, l_e)).mean(dim=['lat','lon','month'])).std(),
			((t_2.prec * t_2.axyp).sel(year = slice(t_s, t_e), lat=slice(l_s, l_e)).mean(dim=['lat','lon','month']) / t_2.axyp.sel(year = slice(t_s, t_e), lat=slice(l_s, l_e)).mean(dim=['lat','lon','month'])).std(),
			((t_3.prec * t_3.axyp).sel(year = slice(t_s, t_e), lat=slice(l_s, l_e)).mean(dim=['lat','lon','month']) / t_3.axyp.sel(year = slice(t_s, t_e), lat=slice(l_s, l_e)).mean(dim=['lat','lon','month'])).std(),
			((t_4.prec * t_4.axyp).sel(year = slice(t_s, t_e), lat=slice(l_s, l_e)).mean(dim=['lat','lon','month']) / t_4.axyp.sel(year = slice(t_s, t_e), lat=slice(l_s, l_e)).mean(dim=['lat','lon','month'])).std() ])
		return e	

		
	### B) P-E=0 ###	

	def y_pe(h):
		"""
		Brings latitude of P-E=0 for each hemisphere for 1,2,3,4xCO2
		
		:Input:		
		 - *h* (str) - hemisphere in either 'SH' or 'NH'
		:Output:
		 - *y* (ndarray) - latitude of P-E = 0 for 1,2,3,4xCO2
		"""
		Mean = pe_1.sel(year = slice(t_s, t_e), season='all_year', pe_location = h).mean()
		y = np.array([	
			(pe_1 - pe_1).sel(year = slice(t_s, t_e), season='all_year', pe_location = h).mean(),
			(pe_2 - pe_1).sel(year = slice(t_s, t_e), season='all_year', pe_location = h).mean(),
			(pe_3 - pe_1).sel(year = slice(t_s, t_e), season='all_year', pe_location = h).mean(),
			(pe_4 - pe_1).sel(year = slice(t_s, t_e), season='all_year', pe_location = h).mean() ])
		return y

	def e_pe(h):
		"""
		Brings 1 \sigma of latitude of P-E=0 for each hemisphere for 1,2,3,4xCO2
		
		:Input:		
		 - *h* (str) - hemisphere in either 'SH' or 'NH'
		:Output:
		 - *e* (ndarray) - 1 \sigma of yearly variations of P-E=0 latitude at a hemisphere for 1,2,3,4xCO2
		"""
		e = np.array([	
			pe_1.sel(year = slice(t_s, t_e), season='all_year', pe_location = h).std(),
			pe_2.sel(year = slice(t_s, t_e), season='all_year', pe_location = h).std(),
			pe_3.sel(year = slice(t_s, t_e), season='all_year', pe_location = h).std(),
			pe_4.sel(year = slice(t_s, t_e), season='all_year', pe_location = h).std() ])
		return e
	

	### C) HC Width ###	

	def y_psi(h):
		"""
		Calculates Hadley Cell Width for 1,2,3,4xCO2 runs
		
		:Inputs:		
		 - *h* (str) - hemisphere in either 'SH' or 'NH'
		:Output:
		 - *y* (ndarray) - HC Width 50 year average for 1,2,3,4xCO2
		"""
		y = np.array([
			(u_1 - u_1).sel(year = slice(t_s, t_e), season = 'all_year', psi_location = h).mean(),
			(u_2 - u_1).sel(year = slice(t_s, t_e), season = 'all_year', psi_location = h).mean(),
			(u_3 - u_1).sel(year = slice(t_s, t_e), season = 'all_year', psi_location = h).mean(),
			(u_4 - u_1).sel(year = slice(t_s, t_e), season = 'all_year', psi_location = h).mean()])
		return y

	def e_psi(h):
		"""
		Calculates 1 \sigma of Hadley Cell width for 1,2,3,4xCO2 
		
		:Input:		
		 - *hemisphere* (str) - hemisphere in either 'SH' or 'NH'	
		:Output:
		 - *e* (ndarray) - 1 \sigma of yearly variations of HC width for 1,2,3,4xCO2
		"""
		e = np.array([
			u_1.sel(year = slice(t_s, t_e), season = 'all_year', psi_location = h).std(),
			u_2.sel(year = slice(t_s, t_e), season = 'all_year', psi_location = h).std(),
			u_3.sel(year = slice(t_s, t_e), season = 'all_year', psi_location = h).std(),
			u_4.sel(year = slice(t_s, t_e), season = 'all_year', psi_location = h).std() ])
		return e
	
		
	### D) HC Strength ### 

	def y_max_psi(ms):
		"""
		Calculates maximum streamfunction over last 50 years of specified months for 1,2,3,4xCO2
		
		:Input:
		 - *ms* (list of strings) -  months specified (e.g. ['JAN','FEB','DEC'])
		:Output:
		 - *y* (ndarray) - HC Strength as 50 year average for 1,2,3,4xCO2
		"""
		y = np.array([
			abs(psi_max.sel(year = slice(t_s, t_e), month = ms, CO2 = 1).mean(dim=['month'])).max('lat').mean('year'),
			abs(psi_max.sel(year = slice(t_s, t_e), month = ms, CO2 = 2).mean(dim=['month'])).max('lat').mean('year'),
			abs(psi_max.sel(year = slice(t_s, t_e), month = ms, CO2 = 3).mean(dim=['month'])).max('lat').mean('year'),
			abs(psi_max.sel(year = slice(t_s, t_e), month = ms, CO2 = 4).mean(dim=['month'])).max('lat').mean('year') ]) 
		return y

	def e_max_psi(ms):
		"""
		Calculates 1 \sigma of  maximum streamfunction over last 50 years of specified months for 1,2,3,4xCO2
		
		:Inputs:
		 - *ms* (list of strings): months specified (e.g. ['JAN','FEB','DEC'])
		:Output:
		 - *e* (ndarray) - 1 \sigma of yearly variations of HC Strength for 1,2,3,4xCO2
		"""
		e = np.array([
			abs(psi_max.sel(year = slice(t_s, t_e), month = ms, CO2 = 1).mean(dim=['month'])).max('lat').std(),
			abs(psi_max.sel(year = slice(t_s, t_e), month = ms, CO2 = 2).mean(dim=['month'])).max('lat').std(),
			abs(psi_max.sel(year = slice(t_s, t_e), month = ms, CO2 = 3).mean(dim=['month'])).max('lat').std(),
			abs(psi_max.sel(year = slice(t_s, t_e), month = ms, CO2 = 4).mean(dim=['month'])).max('lat').std() ]) 
		return e


	### plot ###

	fig = plt.figure()
	fig.set_figwidth(fig.get_figwidth() * 2)
	fig.set_figheight(fig.get_figheight() * 2)
	
	x = 5.35 * np.log(np.array([1, 2, 3, 4])) 
	
	axes = fig.add_subplot(2,2,1)
	axes.errorbar(x, y_prec(0,90), e_prec(0,90), linewidth = 3, marker='o', markersize=12, capsize = 8, color = 'red', label='NH')
	axes.errorbar(x, y_prec(-90,0), e_prec(-90,0), linewidth = 3, marker='o', markersize=12, capsize = 8, color = 'blue', label='SH')
	plt.ylabel('mm/day', fontsize=18)
	plt.title("A) Precipitation", fontsize = 20)
	axes.legend(loc = 0)
	plt.xticks([-4,0,4,8,12], fontsize = 15)
	plt.xlim([-4.3,12])
	plt.yticks([2.7,3.0,3.3], fontsize = 15)
				
	axes = fig.add_subplot(2,2,2)
	axes.errorbar(x, y_pe('NH'), e_pe('NH'), linewidth=3, marker='o', markersize=12, capsize = 8, color = 'red', label='NH')
	axes.errorbar(x, y_pe('SH'), e_pe('SH'), linewidth=3, marker='o', markersize=12, capsize = 8, color = 'blue', label='SH')	
	axes.plot(np.linspace(-4.3,12),np.zeros((50)), color = 'black')
	axes.set_ylabel("latitude (deg.)", fontsize=18)
	axes.set_title("B) P-E=0", fontsize = 20)
	plt.xticks([-4,0,4,8,12], fontsize = 15)
	plt.yticks([-2,0,2], fontsize = 15)
	plt.ylim(-3.3,3.3)	
	plt.xlim([-4.3,12])
	plt.legend(loc=0, fontsize = 10)
	
	axes = fig.add_subplot(2,2,3)
	axes.errorbar(x, y_psi('SH'), e_psi('SH'), linewidth=3, marker='o', markersize=12, capsize = 8, color = 'blue', label='SH')
	axes.errorbar(x, y_psi('NH'), e_psi('NH'), linewidth=3, marker='o', markersize=12, capsize = 8, color = 'red', label='NH')
	axes.plot(np.linspace(-4.3,12),np.zeros((50)), color = 'black')
	axes.set_ylim([-3.2, 3.2])
	axes.set_xlim([-4.3,12])
	axes.set_xlabel("Radiative Forcing (W/m$^2$)", fontsize = 20)
	axes.set_ylabel("latitude (deg.)", fontsize = 20)
	plt.title("C) HC width", fontsize = 20)
	plt.xticks([-4,0,4,8,12], fontsize = 15)
	plt.yticks([-2, 0, 2], fontsize = 15)	
	plt.legend(loc=0, fontsize = 10)
		
	axes = fig.add_subplot(2,2,4)
	axes.errorbar(x, y_max_psi(ms = ['JUN','JUL','AUG']), e_max_psi(ms = ['JUN','JUL','AUG']), linewidth=3, marker='o', markersize=12, capsize = 8, color = 'blue', label='JJA')
	axes.errorbar(x, y_max_psi(ms = ['JAN','FEB','DEC']), e_max_psi(ms = ['JAN','FEB','DEC']), linewidth=3, marker='o', markersize=12, capsize = 8, color = 'red', label='DJF')
	axes.set_xlabel("Radiative Forcing (W/m$^2$)", fontsize = 20)
	axes.set_ylabel("10$^{10}$ kg/s", fontsize = 20)
	axes.set_title("D) HC strength", fontsize = 20)	
	plt.xticks([-4,0,4,8,12], fontsize = 15)
	plt.yticks([16e10, 20e10, 24e10, 28e10], [16, 20, 24, 28], fontsize = 15)
	plt.xlim([-4.3,12])
	axes.legend(loc = 0, fontsize = 10)
	
	plt.tight_layout()
	plt.savefig('sf_figure_4_pd3_QFLUX.pdf')
	plt.show()

figure_4_pd3_QFLUX(t_s = 1930, t_e = 1960)
 
