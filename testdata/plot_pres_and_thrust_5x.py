from nozzle_code_2 import nozzle
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Ellipse
from matplotlib import lines
from scipy import stats
from scipy import signal
import math
import pandas as pd
import numpy as np
import sys
import seaborn as sns
import os

# sns.set()
# sns.axes_style("white")
# sns.set_style("whitegrid", {"xtick.major.size": 0, "ytick.major.size": 0, 'grid.linestyle': '--'})
# sns.set_context("paper", font_scale = 1, rc={"grid.linewidth": .5})
# sns.set_palette("colorblind")
dpi=300

# constrct color map for matplotlib because god damnit this is so annyoying that i can't use seaborn to print the mean and stdev of a few goddamn points
# my_cmap = ListedColormap(sns.color_palette("colorblind"))

class ScalarFormatterForceFormat(mpl.ticker.ScalarFormatter):
		def _set_format(self):  # Override function that finds format to use.
			self.format = "%1.1f"  # Give format here
sns.axes_style("white")
sns.set_style("whitegrid", {"xtick.major.size": 0, "ytick.major.size": 0, 'grid.linestyle': '--'})
sns.set_context("paper", font_scale = 1, rc={"grid.linewidth": .5})
yfmt = ScalarFormatterForceFormat()


# ---=== Steady-state Thrust, mdot vs. Pressure Simulation ===---


## ==================================================================================
## ---- USER OPTIONS ----------------------------------------------------------------
## ==================================================================================

P_amb = 14.7 * 6894.76  			# Ambient Pressure, units of Pa (psia * 6894.76)
d_star = 0.6 / 1000  			# Nozzle throat diameter, units of m (mm / 1000)
expansion_ratio = 1.17			# Nozzle expansion ratio (Exit Area / Throat Area)
half_angle = 10  				# (Conical) Nozzle expansion angle (degrees)
gas_type = 'CO2'				# Gas Choices: R236fa, R134a, N2, CO2, H2, air




## ==================================================================================
## ---- INIT DATA LISTS -------------------------------------------------------------
## ==================================================================================

list_of_P_ts = [(x+14.7)*6894.76 for x in list(range(101))]
list_of_T_ts = [220, 293]
sim_results = pd.DataFrame()

for T_t in list_of_T_ts:
	list_of_thrusts = []
	list_of_mdots = []

	for P_t in list_of_P_ts:
		m_dot, M_crit_sub, M_crit_sup, PR_crit_sub, PR_crit_sup, PR_exit_shock, M_exit_behindshock, M_exit, P_exit, v_exit, F, P_star, T_star, rho_star, Re_star, T_exit, rho_exit = nozzle(P_t, T_t, P_amb, d_star, expansion_ratio, half_angle, gas_type)
		
		list_of_mdots.append(m_dot*1000)
		list_of_thrusts.append(F)

	# Sim results
	df = pd.DataFrame(data={'Pressure (kPa)': [x/1000 for x in list_of_P_ts], 'Temperature (K)':[T_t for x in list_of_P_ts], 'Thrust (mN)': [x*1000 for x in list_of_thrusts], 'mdot (g/s)': list_of_mdots})
	sim_results = sim_results.append(df)






# ---=== Steady-state Thrust & Pressure Data ===---

test_nos = {
	'100 psig': ["20191204_145400", "20191204_145419", "20191204_145442"],
	 '80 psig': ["20191204_144715", "20191204_144754", "20191204_144840"],
	 '60 psig': ["20191204_144142", "20191204_144207", "20191204_144256"],
	 '40 psig': ["20191204_143735", "20191204_143817", "20191204_143931"],
	 '20 psig': ["20191204_143449", "20191204_143538", "20191204_143636"]}

all_float_data = pd.DataFrame()
all_prop_data = pd.DataFrame()
all_thrust_data = pd.DataFrame()
all_data = pd.DataFrame()

for i, setpoint in enumerate(test_nos):
	for j, trial_no in enumerate(test_nos[setpoint]):
		setpoint_kPa = str( round((int(setpoint.split(' ')[0])+14.7)*6.895)) + ' kPa'
		prefix = trial_no
		data_float_psi = pd.read_csv('/home/josh/nozzleDesign/testdata/' + str(prefix + '_float_data.csv'), header=1)																# Load trial float pressure data
		resampled_n = int(round(data_float_psi['Time (s)'][-1:])/0.1)																										# Determine how many points are necessary to have even spacing at 0.1 sec
		resampled_time = np.linspace(0, resampled_n/10, resampled_n+1)																										# Make a new series of times based on the number of samples
		resampled_data = signal.resample(data_float_psi["Float Pressure (psig)"], resampled_n+1)																			# Resample the data
		data_float_psi_resampled = pd.DataFrame({'Time (s)': resampled_time, 'Pressure (kPa)': [(x+14.7)*6.895 for x in resampled_data]})															# Make a DataFrame object out of the resampled data
		data_float_psi_resampled = pd.concat([data_float_psi_resampled, pd.DataFrame([j+1 for x in data_float_psi_resampled.index], columns=['Trial'])], axis=1)			# Concatenate with column of Trial nos.
		# data_float_psi_resampled = pd.concat([data_float_psi_resampled, pd.DataFrame([setpoint for x in data_float_psi_resampled.index], columns=['Set Point (psig)'])], axis=1)	# Concatenate with column of  Ps
		# data_float_all = pd.concat([data_float_psi_resampled, pd.DataFrame([(x+14.7)*6.895 for x in data_float_psi_resampled['Pressure (psig)']], columns=['Pressure (kPa)'])], axis=1)
		data_float_all = pd.concat([data_float_psi_resampled, pd.DataFrame([setpoint_kPa for x in data_float_psi_resampled.index], columns=['Set Point'])], axis=1)

		data_prop_psi = pd.read_csv('/home/josh/nozzleDesign/testdata/' + str(prefix + '_prop_data.csv'), header=1)
		resampled_n = int(round(data_prop_psi['Time (s)'][-1:])/0.1)																										# Determine how many points are necessary to have even spacing at 0.1 sec
		resampled_time = np.linspace(0, resampled_n/10, resampled_n+1)																										# Make a new series of times based on the number of samples
		resampled_data = signal.resample(data_prop_psi["Prop Pressure (psig)"], resampled_n+1)																				# Resample the data
		data_prop_psi_resampled = pd.DataFrame({'Time (s)': resampled_time, 'Pressure (kPa)': [(x+14.7)*6.895 for x in resampled_data]})																# Make a DataFrame object out of the resampled data
		data_prop_psi_resampled = pd.concat([data_prop_psi_resampled, pd.DataFrame([j+1 for x in data_prop_psi_resampled.index], columns=['Trial'])], axis=1)				# Concatenate with column of Trial nos.
		# data_prop_psi_resampled = pd.concat([data_prop_psi_resampled, pd.DataFrame([setpoint for x in data_prop_psi_resampled.index], columns=['Set Point (psig)'])], axis=1)		# Concatenate with column of Setpoints
		# data_prop_all = pd.concat([data_prop_psi_resampled, pd.DataFrame([(x+14.7)*6.895 for x in data_prop_psi_resampled['Pressure (psig)']], columns=['Pressure (kPa)'])], axis=1)
		data_prop_all = pd.concat([data_prop_psi_resampled, pd.DataFrame([setpoint_kPa for x in data_prop_psi_resampled.index], columns=['Set Point'])], axis=1)

		data_thrust = pd.read_csv('/home/josh/nozzleDesign/testdata/' + str(prefix + '_loadcell_data.csv'), header=1)
		data_thrust.insert(2, "Thrust (mN)", [x*9.81 for x in data_thrust["Weight (?)"]])
		resampled_n = int(round(data_thrust['Time (s)'][-1:])/0.1)																											# Determine how many points are necessary to have even spacing at 0.1 sec
		resampled_time = np.linspace(0, resampled_n/10, resampled_n+1)																										# Make a new series of times based on the number of samples
		resampled_data = signal.resample(data_thrust["Thrust (mN)"], resampled_n+1)																							# Resample the data
		data_thrust_resampled = pd.DataFrame({'Time (s)': resampled_time, 'Thrust (mN)': resampled_data})																	# Make a DataFrame object out of the resampled data
		data_thrust_resampled = pd.concat([data_thrust_resampled, pd.DataFrame([j+1 for x in data_thrust_resampled.index], columns=['Trial'])], axis=1)						# Concatenate with column of Trial nos.
		data_thrust_resampled = pd.concat([data_thrust_resampled, pd.DataFrame([setpoint_kPa for x in data_thrust_resampled.index], columns=['Set Point'])], axis=1)				# Concatenate with column of Setpoints

		all_float_data = all_float_data.append(data_float_all)
		all_prop_data  = all_prop_data.append(data_prop_all)
		all_thrust_data = all_thrust_data.append(data_thrust_resampled)


df_pres_at_2sec = all_float_data.loc[all_float_data['Time (s)'] == 2.0].reset_index(drop=True)
df_thrust_at_2sec = all_thrust_data.loc[all_thrust_data['Time (s)'] == 2.0].reset_index(drop=True)

avg_of_data = df_pres_at_2sec.merge(df_thrust_at_2sec, on=['Set Point', 'Trial', 'Time (s)'])
avg_of_data.drop(['Time (s)'], axis=1, inplace=True)

thrust_means = avg_of_data.groupby('Set Point').mean().drop(['Trial'], axis=1).sort_values(['Pressure (kPa)'], ascending=False).reset_index()
thrust_stdevs = avg_of_data.groupby('Set Point').std().drop(['Trial'], axis=1).sort_values(['Pressure (kPa)'], ascending=False).reset_index()

# Now we have one DataFrame with mean + stdev thrust & pressure for each setpoint
# I want the same for mdot now





# ---=== Steady-state Mass Flow Rate & Pressure Data ===---

data1 = pd.read_csv('/home/josh/nozzleDesign/testdata/11132019_test11.csv')
data1.insert(3, "P init (psia)", [x+14.7 for x in data1["P init (psig)"]], True)
data1.insert(5, "P final (psia)", [x+14.7 for x in data1["P final (psig)"]], True)
data1.insert(6, "P avg (psia)", [x for x in (data1["P init (psia)"] + data1["P final (psia)"])/2], True)
data1.insert(7, "P avg (psig)", [x-14.7 for x in data1["P avg (psia)"]], True)
data1.insert(8, "P act_est (psig)", [x for x in data1["P init (psig)"]], True)
data1.insert(9, "P act_est (kPa)", [(x+14.7)*6.895 for x in data1["P act_est (psig)"]], True)
data1.insert(12, "dM", [round(x,2) for x in data1['M final (g)']-data1['M init (g)']])
data1.insert(13, "dM err", [round(np.sqrt(0.02),2) for x in data1["dM"]])
# Adding or subtracting uncertain data means the uncertanties add *in quadrature* (sqaure-root of the sum of the squares)
# This assumes the uncertainy provided by the manufacturer is given in terms of standard deviation
# Thus, the error is: +/- sqrt(0.02) grams for all data points

data1.insert(14, "dM/dt", [-round(x,3) for x in data1["dM"]/data1["Time (s)"]])
data1.insert(15, "dM/dt err", [round(x,3) for x in np.sqrt((np.sqrt(0.02)/data1["dM"])**2 + ((4E-6)/data1["Time (s)"])**2)])
# Multiplying or dividing uncertain data is a little more complicated
# All measurements must be put in the form of *fractional* uncertainties, then added in quadrature
# Thus the error is: +/- sqrt((dm_err/dm)**2 + (dt_err/dt)**2) grams/sec for each data point

mdot_means = data1.groupby('Set Point').mean().drop(['Trial'], axis=1).sort_values(['Set Point'], ascending=False).reset_index()
mdot_stdevs = data1.groupby('Set Point').std().drop(['Trial'], axis=1).sort_values(['Set Point'], ascending=False).reset_index()





# ---=== Plot it! ===---
# We're going to make a two 2x1 subplot plots
# First one is (1) Thrust vs. Time at 5 pressure setpoints, and (2) Pressure vs. Time at 5 pressure setpoints
# Second one (1) Thrust vs. Pressure, and (2) mdot vs. Pressure


# First one: (1) Thrust vs. Time, (2) Pressure vs. Time
# This is ALL experimental data ONLY!!!!
fig1, axs = plt.subplots(2, sharex=True, dpi=dpi, figsize=(6, 4))
fig1.suptitle('Pressure & Thrust vs. Time for Steady-State Conditions', y=0.98)
axs[0].set_title(r'(CO$_2$, Nozzle $\varnothing$0.6 mm, $\lambda$=1.34 [3$x$ Trials/Set Point])', fontsize=8)
fig1.canvas.set_window_title('SteadyState Pressure and Thrust UpdatedStyle')

sns.lineplot(ax=axs[0], x="Time (s)", y="Pressure (kPa)", hue="Set Point", style="Set Point", data=all_float_data, ci='sd')
sns.lineplot(ax=axs[1], x="Time (s)", y="Thrust (mN)", hue="Set Point", style="Set Point", data=all_thrust_data, ci='sd')

for ax in axs:
	ax.yaxis.set_major_formatter(yfmt)
	ax.yaxis.offsetText.set_fontsize(7)
	ax.tick_params(axis='x', colors='#413839', labelsize=7, pad=0)
	ax.tick_params(axis='y', labelsize=7, pad=0)
	ax.grid(which='major', axis='both', linestyle='--')
	ax.xaxis.label.set_size(8)
	ax.yaxis.label.set_size(8)
	ax.legend(loc=1, fontsize=8, framealpha=0.9)

axs[0].set(ylabel=r'Pressure, $kPa$')
axs[0].set_ylim(bottom=-30)

axs[1].set(ylabel=r'Thrust, $mN$')
axs[1].set_ylim(bottom=-10)

axs[-1].set(xlabel=r'Time, $sec$')

plt.tight_layout()
plt.subplots_adjust(top=0.875)
fig1.align_ylabels()

plt.xlim(left=0, right=16)


# Second one: (1) Thrust vs. Pressure (Isentropic + Experimental), (2) mdot vs. Pressure (Isentropic + Experimental)
fig2, axs = plt.subplots(2, 1, sharex=True, dpi=dpi, figsize=(6,4))
fig2.suptitle('Thrust & Mass Flow Rate vs. Pressure for Steady-State Conditions', y=0.98)
axs[0].set_title(r'(CO$_2$, Nozzle $\varnothing$0.6 mm, $\lambda$=1.34 [3$x$ Trials/Set Point)]', fontsize=8)
fig2.canvas.set_window_title('SteadyState Thrust and MassFlowRate UpdatedStyle')

# ====== Plot (1): Thrust vs Pressure ============
axs[0].legend()
for i, T_t in enumerate(list_of_T_ts):
	label = 'Isentropic (' + str(T_t) + ' K)'
	axs[0].plot("Pressure (kPa)", "Thrust (mN)", data=sim_results.loc[sim_results['Temperature (K)'] == T_t], label=label, linestyle=list(lines.lineStyles.keys())[i])
sns.scatterplot(ax=axs[0], x="Pressure (kPa)", y="Thrust (mN)", data=thrust_means, s=6**2, marker='o', label='Experimental')


# ====== Plot (2): mdot vs Pressure ==============
axs[1].legend()
for i, T_t in enumerate(list_of_T_ts):
	label = 'Isentropic (' + str(T_t) + ' K)'
	axs[1].plot("Pressure (kPa)", "mdot (g/s)", data=sim_results.loc[sim_results['Temperature (K)'] == T_t], label=label, linestyle=list(lines.lineStyles.keys())[i])
sns.scatterplot(ax=axs[1], x="P act_est (kPa)", y="dM/dt", data=mdot_means, s=6**2, marker='o', label='Experimental')

for ax in axs:
	ax.yaxis.set_major_formatter(yfmt)
	ax.yaxis.offsetText.set_fontsize(7)
	ax.tick_params(axis='x', colors='#413839', labelsize=7, pad=0)
	ax.tick_params(axis='y', labelsize=7, pad=0)
	ax.grid(which='major', axis='both', linestyle='--')
	ax.xaxis.label.set_size(8)
	ax.yaxis.label.set_size(8)
	ax.legend(loc=4, fontsize=8, framealpha=0.9)

axs[0].set(ylabel=r'Thrust, $mN$')
axs[0].set_ylim(bottom=0, top=300)

axs[1].set(ylabel=r'$\dot{m}$, $g/s$')
axs[1].set_ylim(bottom=0, top=0.8)

axs[-1].set(xlabel=r'Pressure, $kPa$')
axs[-1].set_xlim(left=100, right=800)

plt.tight_layout()
plt.subplots_adjust(top=0.875)
fig2.align_ylabels()

plt.show()