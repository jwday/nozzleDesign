# plot_ideal_vs_exp.py
# Comparison of ideal (simulated) results with experimental data for single plenum discharge
# Requires user to first run simulation to generate 'simulation_data.csv' (if not already present, or if parameters are changed)
# Will plot 2x1 subplots to compare pressure and thrust to experimental data collected on UC Davis HRVIP thrust test stand, 11-30-2019

from nozzle_code_3 import nozzle
import matplotlib.pyplot as plt
from scipy import stats
import math
import pandas as pd
import numpy as np
import sys
import seaborn as sns
from data_handling_funcs import all_data as all_exp_data
import matplotlib as mpl

all_data_heat_xfer = pd.read_csv('simulation_data.csv')
all_data_isen = pd.read_csv('simulation_data.csv')
d_star = 0.6 / 1000
expansion_ratio = 1.34
gas_label = all_data_heat_xfer[all_data_heat_xfer['gas_type'] == 'CO2']['gas_type'][0]
T_t_init = all_data_heat_xfer[all_data_heat_xfer['gas_type'] == 'CO2']['T_t'][0].round(2)

class ScalarFormatterForceFormat(mpl.ticker.ScalarFormatter):
		def _set_format(self):  # Override function that finds format to use.
			self.format = "%1.1f"  # Give format here
sns.axes_style("white")
sns.set_style("whitegrid", {"xtick.major.size": 0, "ytick.major.size": 0, 'grid.linestyle': '--'})
sns.set_context("paper", font_scale = 1, rc={"grid.linewidth": .5})
yfmt = ScalarFormatterForceFormat()

## ==================================================================================
## ---- PRESSURE & THRUST VS TIME ---------------------------------------------------
## ==================================================================================

## ---- Steady-state thrust tests on load cell --------------------------------------

# test_nos = ['20191205_191138', # 114.7 psia, 0.6 mm nozzle, raw data in g (multiply by 9.81)
# 			'20191205_191210', 
# 			'20191205_191332', 
# 			'20191205_191402', 
# 			'20191205_191433'
# 			]  
# steady_state = True
# mult_by_g = True

# test_nos = [ 
# 			'20191204_143449',  # 3x trials @ 5x pressures, 0.6 mm nozzle, raw data in g (multiply by 9.81)
# 			'20191204_143538',
# 			'20191204_143636',

# 			'20191204_143735',
# 			'20191204_143817',
# 			'20191204_143931',

# 			'20191204_144142',
# 			'20191204_144207',
# 			'20191204_144256',

# 			'20191204_144715',
# 			'20191204_144754',
# 			'20191204_144840',

# 			'20191204_145400',
# 			'20191204_145419',
# 			'20191204_145442'
# 			]  
# steady_state = True
# mult_by_g = True

# test_nos = [ 
# 			'20191205_191138',  # 5x trials @ 114.7 psia, 0.6 mm nozzle, raw data in g (multiply by 9.81)
# 			'20191205_191210',
# 			'20191205_191332',
# 			'20191205_191402',
# 			'20191205_191433'
# 			]  
# steady_state = True
# mult_by_g = True

# test_nos = [ 
# 			'20191219_205802',  # Trash
# 			'20191219_205915',  # Trash
# 			'20191219_205943',  # Trash
# 			]  
# steady_state = False
# mult_by_g = True




# ## ---- Single plenum discharge tests -----------------------------------------------

test_nos = [
			'20191130_131419', # 90 psig, 0.6 mm nozzle, raw data in g (multiply by 9.81)
			'20191130_131515',
			'20191130_131607',
			'20191130_131624',
			'20191130_131644'
			]  
steady_state = False
mult_by_g = True


# test_nos = [
# 			'20191223_183658', # 114.7 psia, 0.4 mm nozzle, raw data in mN (do not multiply by 9.81)
# 			'20191223_183725',
# 			'20191223_183832',
# 			'20191223_183908',
# 			'20191223_183945'
# 			]  
# steady_state = False
# mult_by_g = False

## ----------------------------------------------------------------------------------

linewidth = 1.5
fontsize = 10
if steady_state:
	data_marker = 'None'
else:
	data_marker = 'o'

fig1, axs = plt.subplots(2, 1, figsize=(6,4), dpi=300, sharex='col')
if steady_state:
	fig1.suptitle('Steady-State Pressure & Thrust Measurements'.format(len(test_nos), d_star*1000, expansion_ratio))
else:
	fig1.suptitle('Single Plenum Discharge Pressure & Thrust Measurements'.format(len(test_nos), d_star*1000, expansion_ratio), y=0.98)
	# fig1.suptitle('Simulation Results Convergence for Varying Time Step', y=0.95)
	# axs[0].set_title(r'({}, $V_{{p}}=${} cm$^3$, Nozzle $\varnothing${} mm, $\lambda$={})'.format(gas_label, vol*10**6, d_star*1000, expansion_ratio), fontsize=7)
	# axs[0].set_title(r'CO$_2$, Nozzle $\varnothing$0.6 mm, $\lambda$=1.34 (5$x$ Trials/Set Point)', fontsize=7, color='dimgray', y=1.03)
	axs[0].set_title(r'({} at $T_0$={} K, $V_{{p}}=${} cm$^3$, Nozzle $\varnothing =${} mm, $\lambda$={})'.format(gas_label, T_t_init, 30, d_star*1000, expansion_ratio), fontsize=8)
fig1.canvas.set_window_title('Single Plenum Discharge ExpVsIdeal')


td1 = []  # Pressure
for trial, test_no in reversed(list(enumerate(test_nos))):
	test_data = all_exp_data('./testdata/' + test_no, mult_by_g)[0]  # All the Float Pressure data
	test_data['Trial'] = trial  # Used to label the data for showing individually
	test_data['Setpoint'] = '{} kPa'.format(int(test_data['Float Pressure (kPa)'][0].round(-1)))  # Used to label the data for showing individually (converted to kPa)
	td1.append(test_data)
td1 = pd.concat(td1)
td1['Time (s)'] = td1['Time (s)'].round(1)
sns.lineplot(ax=axs[0],
			 x='Time (s)',
			 y='Float Pressure (kPa)',
			 data=td1,
			#  hue='Setpoint',  	# Show each trial individually instead of an aggregate
			#  style='Setpoint',  	# Show each trial individually instead of an aggregate
			#  estimator=np.mean,
			 marker=data_marker,
			 label='Experimental')
			 
if not steady_state:
	# sns.lineplot(ax=axs[0],
	# 			 data=all_data_heat_xfer,
	# 			 x='time',
	# 			 y='P_t',
	# 			 hue='Time Step'
	# )
	axs[0].plot([x+0.68 for x in all_data_heat_xfer[all_data_heat_xfer['gas_type']=='CO2']['time']], [x/1000 for x in all_data_heat_xfer[all_data_heat_xfer['gas_type']=='CO2']['P_t']], color='#ff7f0e', label='Isentropic Theory', linestyle='-', linewidth=linewidth)
	# axs[0].plot([x+0.68 for x in all_data_isen[all_data_isen['gas_type']=='CO2']['time']], [x/1000 for x in all_data_isen[all_data_isen['gas_type']=='CO2']['P_t']], color='#ff7f0e', label='Isentropic', linestyle='--', linewidth=linewidth)
	axs[0].set_xlim(left=0.5, right=2.25)
	axs[0].set_ylim(bottom=0)



td2 = []  # Thrust
for trial, test_no in reversed(list(enumerate(test_nos))):
	test_data = all_exp_data('./testdata/' + test_no, mult_by_g)[2]  # All the Thrust data
	test_data['Trial'] = trial  # Used to label the data for showing individually
	test_data['Setpoint'] = '{} kPa'.format(int(all_exp_data('./testdata/' + test_no, mult_by_g)[0]['Float Pressure (kPa)'][0].round(-1)))  # Used to label the data for showing individually
	td2.append(test_data)
td2 = pd.concat(td2)
td2['Time (s)'] = td2['Time (s)'].round(1)
sns.lineplot(ax=axs[1],
			 x='Time (s)',
			 y='Thrust Corrected (mN)',
			 data=td2,
			#  hue='Setpoint',  	# Show each trial individually instead of an aggregate
			#  style='Setpoint',  	# Show each trial individually instead of an aggregate
			#  estimator=np.mean,  
			 marker=data_marker,
			 label='Experimental')

if not steady_state:
	# sns.lineplot(ax=axs[1],
	# 			x='time',
	# 			y='thrust',
	# 			data=all_data_heat_xfer,
	# 			hue='Time Step'
	# )
	axs[1].plot([x+0.68 for x in all_data_heat_xfer[all_data_heat_xfer['gas_type']=='CO2']['time']], all_data_heat_xfer[all_data_heat_xfer['gas_type']=='CO2']['thrust'].multiply(1000), color='#ff7f0e', label='Isentropic Theory', linestyle='-', linewidth=linewidth)
	# axs[1].plot([x+0.68 for x in all_data_isen[all_data_isen['gas_type']=='CO2']['time']], all_data_isen[all_data_isen['gas_type']=='CO2']['thrust'].multiply(1000), color='#ff7f0e', label='Isentropic', linestyle='--', linewidth=linewidth)
	axs[1].set_xlim(left=0.5, right=2.25)
	axs[1].set_ylim(bottom=0)


for ax in axs:
	ax.yaxis.set_major_formatter(yfmt)
	ax.yaxis.offsetText.set_fontsize(7)
	ax.tick_params(axis='x', colors='#413839', labelsize=7, pad=0)
	ax.tick_params(axis='y', labelsize=7, pad=0)
	ax.grid(which='major', axis='both', linestyle='--')
	ax.xaxis.label.set_size(8)
	ax.yaxis.label.set_size(8)
	ax.legend(loc=1, fontsize=8)

axs[0].set(ylabel=r'Pressure, $kPa$')
axs[0].set_ylim(bottom=-30)

axs[1].set(ylabel=r'Thrust, $mN$')
axs[1].set_ylim(bottom=-10)

axs[-1].set(xlabel=r'Time $(sec)$')

# Net impulse
# if not steady_state:
# 	sns.lineplot(ax=axs[2],
# 				x='time',
# 				y='net impulse',
# 				data=all_data_heat_xfer,
# 				hue='Time Step'
# 	)
# 	# axs[1].plot([x+0.57 for x in all_data_heat_xfer[str(time_step)]['time']], [x *1000 for x in all_data_heat_xfer[str(time_step)]['thrust']], color='#ff7f0e', label='thrust', linestyle='-', linewidth=linewidth)
# 	axs[2].set_xlim(left=0, right=right_limit)
# 	axs[2].set_ylim(bottom=0)

# axs[2].set_xlabel('Time, s', color='#413839', fontsize=fontsize)
# axs[2].set_ylabel('Impulse, $mN-s$', color='#413839', fontsize=fontsize)
# axs[2].tick_params(colors='#413839')
# axs[2].grid(which='major', axis='both', linestyle='--')
# axs[2].legend(loc=1, fontsize=9)


# box1 = axs[2].get_position()
# axs[1].set_position([box1.x0 + box1.width * 0.05, box1.y0 + box1.height * 0.05, box1.width, box1.height])




## ==================================================================================
## ---- THRUST VS PRESSURE ----------------------------------------------------------
## ==================================================================================
# test_nos = [ 
# 			'20191204_143449',  # 3x trials @ 5x pressures, 0.6 mm nozzle, raw data in g (multiply by 9.81)
# 			'20191204_143538',
# 			'20191204_143636',

# 			'20191204_143735',
# 			'20191204_143817',
# 			'20191204_143931',

# 			'20191204_144142',
# 			'20191204_144207',
# 			'20191204_144256',

# 			'20191204_144715',
# 			'20191204_144754',
# 			'20191204_144840',

# 			'20191204_145400',
# 			'20191204_145419',
# 			'20191204_145442'
# 			]  
# steady_state = True
# mult_by_g = True

# linewidth = 2
# fontsize = 12
# if steady_state:
# 	data_marker = 'None'
# else:
# 	data_marker = 'o'

# td1 = []  # Pressure
# for trial, test_no in reversed(list(enumerate(test_nos))):
# 	test_data = all_data_heat_xfer('Test Data/' + test_no, mult_by_g)[0]  # All the Float Pressure data
# 	test_data['trial'] = trial  # Used to label the data for showing individually
# 	test_data['Setpoint'] = '{} psig'.format(int(test_data['Float Pressure (psig)'][0].round(-1)))  # Used to label the data for showing individually
# 	td1.append(test_data)
# td1 = pd.concat(td1)
# pressure_data = pd.concat([td1['Float Pressure Resampled (psig)'][12], td1['Float Pressure Resampled (psig)'][22]]) # psig, at 1.2 and 2.2 seconds

# td2 = []  # Thrust
# for trial, test_no in reversed(list(enumerate(test_nos))):
# 	test_data = all_data_heat_xfer('Test Data/' + test_no, mult_by_g)[2]  # All the Thrust data
# 	test_data['trial'] = trial  # Used to label the data for showing individually
# 	test_data['Setpoint'] = '{} psig'.format(int(all_data_heat_xfer('Test Data/' + test_no, mult_by_g)[0]['Float Pressure (psig)'][0].round(-1)))  # Used to label the data for showing individually
# 	td2.append(test_data)
# td2 = pd.concat(td2)
# thrust_data = pd.concat([td2['Thrust Resampled (mN)'][12], td2['Thrust Resampled (mN)'][22]])  # mN, at 1.2 and 2.2 seconds


# fig8, ax1 = plt.subplots(figsize=figsize, dpi=dpi)
# 	# Blue: #1f77b4 (Inviscid)
# 	# Green: #2ca02c
# 	# Orange: #ff7f0e
# 	# Gray: #808080
# ax1.set_xlabel('Pressure (psig)', color='#413839')
# ax1.set_ylabel('Thrust (mN)', color='#413839')
# colors = ['#1f77b4', '#ff7f0e']
# linestyles = ['-', ':']

# for i, T in enumerate(list_of_Tts):
# 	label = 'Inviscid (' + str(T-273) + ' C)'
# 	pres = [(x/6894.76)-14.7 for x in list_of_P_ts]
# 	thrust = list(all_data_heat_xfer[str(T)]['thrust']*1000)
# 	ax1.plot(pres, thrust, color=colors[i], label=label, linestyle=linestyles[i])

# ax1.plot([x for x in pressure_data], thrust_data, color='#2ca02c', linestyle='none', label='Experimental', marker='x')

# ax1.tick_params(colors='#413839')
# ax1.grid(which='major', axis='both', linestyle='--')
# ax1.set_xlim(left=0)
# ax1.set_ylim(bottom=0)

# box = ax1.get_position()
# ax1.set_position([box.x0, box.y0 + box.height*0.1, box.width, box.height*0.9])
# fig8.legend(loc='center', bbox_to_anchor=(0.5, 0.03), ncol=3, frameon=False )

# # fig8.legend(['Inviscid', 'Experimental'], loc='center', bbox_to_anchor=(0.5, 0.03), ncol=3, frameon=False )
# plt.title('Ideal & Measured Thrust vs. Pressure\n($\\varnothing${0} mm, $\\lambda$={1})'.format(d_star*1000, expansion_ratio), y=1.0, color='#413839')



# ---- Plot Everything
plt.tight_layout()
plt.subplots_adjust(top=0.885)
fig1.align_ylabels()

plt.show()
