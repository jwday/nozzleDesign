import matplotlib.pyplot as plt
from scipy import stats
import math
import pandas as pd
import numpy as np
import sys
import seaborn as sns
import matplotlib as mpl


all_data = pd.read_csv('simulation_data.csv')
params = pd.read_csv('simulation_parameters.csv')
figsize = (5,3.5)
dpi = 300

data = 	{	'P_t': 			[0, all_data['time'], 	[x*10_6 for x in all_data['P_t']], 		'Plenum Total Pres, $Pa$', 		'-'],
			'thrust': 		[1, all_data['time'], 	[x*1000 for x in all_data['thrust']], 	'Thrust, $mN$', 				'-'],
			'avg_thrust': 	[1, all_data['time'], 	all_data['avg_thrust'], 				'Time Average Thrust, $mN$',	'--'],
			'cum_impulse': 	[2, all_data['time'], 	all_data['cum_impulse'], 				'Impulse, $mN-s$', 				'-']
		}

ylabels = {	'P_t': 			'Pressure, $kPa$',
			'thrust': 		'Thrust, $mN$',
			'cum_impulse': 	'Impulse, $mN-s$'
		}


class ScalarFormatterForceFormat(mpl.ticker.ScalarFormatter):
		def _set_format(self):  # Override function that finds format to use.
			self.format = "%1.1f"  # Give format here
sns.axes_style("white")
sns.set_style("whitegrid", {"xtick.major.size": 0, "ytick.major.size": 0, 'grid.linestyle': '--'})
sns.set_context("paper", font_scale = 1, rc={"grid.linewidth": .5})


for i, gas in enumerate(params['gas_type']):
	gas_label = params['gas_type'][i]
	T_t_init = params['T_t_init'][i].round(2)
	gas_label = params['Propellant'][i]
	plenum_vol = params['vol'][i]*10**6
	d_star = params['d_star'][i]*1000
	expansion_ratio = params['expansion_ratio'][i]

	fig, axs = plt.subplots(3,1, figsize=figsize, dpi=dpi, sharex=True)
	if gas == 'CO2':
		fig.suptitle('{} Single Plenum Discharge Performance'.format('On-Ground'), y=0.98)
	else:
		fig.suptitle('{} Single Plenum Discharge Performance'.format('In-Space'), y=0.98)
	axs[0].set_title(r'({} at $T_0$={} K, $V_{{p}}=${} cm$^3$, Nozzle $\varnothing${} mm, $\lambda$={})'.format(gas_label, T_t_init, round(plenum_vol, 1), d_star, expansion_ratio), fontsize=7)
	fig.canvas.set_window_title('Performance Specs')

	axs[0].plot(all_data[all_data['gas_type'] == gas]['time'], [x*10**-3 for x in all_data[all_data['gas_type'] == gas]['P_t']], 											linestyle='-', color='#1f77b4')
	axs[1].plot(all_data[all_data['gas_type'] == gas]['time'], [x*1000 for x in all_data[all_data['gas_type'] == gas]['thrust']], 		label='Instantaneous Thrust',		linestyle='-', color='#1f77b4')
	axs[1].plot(all_data[all_data['gas_type'] == gas]['time'], [x*1000 for x in all_data[all_data['gas_type'] == gas]['avg_thrust']], 	label='Cumulative Average Thrust',	linestyle='--', color='#1f77b4')
	axs[2].plot('time', 'cum_impulse', 	data=all_data[all_data['gas_type'] == gas], linestyle='-', color='#1f77b4')

	for i, key in enumerate(ylabels.keys()):
		axs[i].set_ylabel(ylabels[key], color='#413839', fontsize=8)
		if gas == 'CO2':
			axs[i].set_xlim([0, 1.4])
		else:
			axs[i].set_xlim([0, 28])
		# axs[i].legend().texts[0].set_text('Flow Regimes')

		yfmt = ScalarFormatterForceFormat()
		# yfmt.set_powerlimits((0,0))
		axs[i].yaxis.set_major_formatter(yfmt)
		# axs[i].ticklabel_format(axis='y', style='sci', scilimits=(0,0), useMathText=True)
		axs[i].tick_params(axis='y', labelsize=6, pad=0)
		axs[i].yaxis.offsetText.set_fontsize(6)

		axs[i].tick_params(axis='x', labelsize=6, pad=0)

	axs[1].legend(loc='upper right', fontsize=6, framealpha=0.9)
	axs[2].xaxis.label.set_size(8)
	axs[2].set(xlabel=r'Time $(sec)$')

	plt.tight_layout()
	plt.subplots_adjust(top=0.88, hspace=0.115)
	fig.align_ylabels()
	# axs[0].set_ylim([-0.01, 2.6])
	# axs[1].set_ylim(bottom=-0.01)
	plt.show()