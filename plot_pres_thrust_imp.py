import matplotlib.pyplot as plt
from scipy import stats
import math
import pandas as pd
import numpy as np
import sys
import seaborn as sns
import matplotlib as mpl

## ==================================================================================
## ---- LOAD DATA AND PARAMETERS ----------------------------------------------------
## ==================================================================================
all_data = pd.read_csv('simulation_data.csv')
params = pd.read_csv('simulation_parameters.csv')

dpi = 300


## ==================================================================================
## ---- THIS BS NEEDS TO BE DONE BETTER ---------------------------------------------
## ==================================================================================
data = 	{	'P_t': 			[0, all_data['time'], 	[x*10_6 for x in all_data['P_t']], 		'Plenum Total Pres, $Pa$', 		'-'],
			'thrust': 		[1, all_data['time'], 	[x*1000 for x in all_data['thrust']], 	'Thrust, $mN$', 				'-'],
			'avg_thrust': 	[1, all_data['time'], 	all_data['avg_thrust'], 				'Time Average Thrust, $mN$',	'--'],
			'cum_impulse': 	[2, all_data['time'], 	all_data['cum_impulse'], 				'Impulse, $mN-s$', 				'-']
		}



## ==================================================================================
## ---- PLOT FORMATTER --------------------------------------------------------------
## ==================================================================================
class ScalarFormatterForceFormat(mpl.ticker.ScalarFormatter):
		def _set_format(self):  # Override function that finds format to use.
			self.format = "%1.1f"  # Give format here
sns.axes_style("white")
sns.set_style("whitegrid", {"xtick.major.size": 0, "ytick.major.size": 0, 'grid.linestyle': '--'})
sns.set_context("paper", font_scale = 1, rc={"grid.linewidth": .5})
yfmt = ScalarFormatterForceFormat()



## ==================================================================================
## ---- JUST DO IT ------------------------------------------------------------------
## ==================================================================================
for i, gas in enumerate(params['gas_type']):
	# -------------------------------------------------------------------------------
	# Setup Common Things
	T_t_init = params['T_t_init'][i].round(2)
	gas_label = params['Propellant'][i]
	plenum_vol = params['vol'][i]*10**6
	d_star = params['d_star'][i]*1000
	expansion_ratio = params['expansion_ratio'][i]

	if gas == 'CO2':
		suptitle_var = 'On-Ground'
	else:
		suptitle_var = 'In-Space'


	# -------------------------------------------------------------------------------
	# Pressure, Thrust, and Net Impulse 
	fig1, axs = plt.subplots(3,1, figsize=(6,4), dpi=dpi, sharex=True)
	fig1.suptitle('{} Single Plenum Discharge Performance'.format(suptitle_var), y=0.98)
	axs[0].set_title(r'({} at $T_0$={} K, $V_{{p}}=${} cm$^3$, Nozzle $\varnothing =${} mm, $\lambda$={})'.format(gas_label, T_t_init, round(plenum_vol, 1), d_star, expansion_ratio), fontsize=8)
	fig1.canvas.set_window_title('Sim Thrust and Impulse {}'.format(gas))

	axs[0].plot(all_data[all_data['gas_type'] == gas]['time'], [x*10**-3 for x in all_data[all_data['gas_type'] == gas]['P_t']], 											linestyle='-', color='#1f77b4')
	axs[1].plot(all_data[all_data['gas_type'] == gas]['time'], [x*1000 for x in all_data[all_data['gas_type'] == gas]['thrust']], 		label='Instantaneous Thrust',		linestyle='-', color='#1f77b4')
	axs[1].plot(all_data[all_data['gas_type'] == gas]['time'], [x*1000 for x in all_data[all_data['gas_type'] == gas]['avg_thrust']], 	label='Cumulative Average Thrust',	linestyle='--', color='#1f77b4')
	axs[2].plot('time', 'cum_impulse', 	data=all_data[all_data['gas_type'] == gas], linestyle='-', color='#1f77b4')

	ylabels = {	'P_t': 			'Pressure, $kPa$',
				'thrust': 		'Thrust, $mN$',
				'cum_impulse': 	'Impulse, $mN-s$'
			}

	for i, key in enumerate(ylabels.keys()):
		axs[i].set(ylabel=ylabels[key])
		if gas == 'CO2':
			axs[i].set_xlim([0, 1.4])
		else:
			axs[i].set_xlim([0, 28])

		axs[i].yaxis.set_major_formatter(yfmt)
		axs[i].yaxis.offsetText.set_fontsize(7)
		axs[i].tick_params(axis='x', labelsize=7, pad=0)
		axs[i].tick_params(axis='y', labelsize=7, pad=0)
		axs[i].grid(which='major', axis='both', linestyle='--')
		axs[i].xaxis.label.set_size(8)
		axs[i].yaxis.label.set_size(8)

	axs[1].legend(loc='upper right', fontsize=8, framealpha=0.9)
	axs[0].set_ylim(bottom=-30)
	axs[-1].set(xlabel=r'Time $(sec)$')

	plt.tight_layout()
	plt.subplots_adjust(top=0.88, hspace=0.115)
	fig1.align_ylabels()
	# axs[1].set_ylim(bottom=-0.01)



	# -------------------------------------------------------------------------------
	# Specific Impulse
	
	fig2, axs = plt.subplots(1,1, figsize=(6, 2), dpi=dpi, sharex=True)
	fig2.suptitle('{} Specific Impulse'.format(suptitle_var), y=0.98)
	axs.set_title(r'({} at $T_0$={} K, $V_{{p}}=${} cm$^3$, Nozzle $\varnothing =${} mm, $\lambda$={})'.format(gas_label, T_t_init, round(plenum_vol, 1), d_star, expansion_ratio), fontsize=8)
	fig2.canvas.set_window_title('Sim ISP {}'.format(gas))

	axs.plot(all_data[all_data['gas_type'] == gas]['time'], all_data[all_data['gas_type'] == gas]['ISP'], linestyle='-', color='#1f77b4')

	if gas == 'CO2':
		axs.set_xlim([0, 1.4])
	else:
		axs.set_xlim([0, 28])
		
	axs.yaxis.set_major_formatter(yfmt)
	axs.yaxis.offsetText.set_fontsize(7)
	axs.tick_params(axis='x', labelsize=7, pad=0)
	axs.tick_params(axis='y', labelsize=7, pad=0)
	axs.xaxis.label.set_size(8)
	axs.yaxis.label.set_size(8)
	axs.set(xlabel=r'Time $(sec)$')
	axs.set(ylabel=r'Specific Impulse, $s$')

	axs.set_ylim(top=50, bottom=-1)

	plt.tight_layout()
	plt.subplots_adjust(top=0.8)
	fig2.align_ylabels()

plt.show()