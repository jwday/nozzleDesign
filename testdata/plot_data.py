import matplotlib.pyplot as plt
from scipy import stats
from scipy import signal
import math
import pandas as pd
import numpy as np
import sys
import seaborn as sns
import os


def all_data(prefix):
	data_float_psi = pd.read_csv('/home/josh/remoteProp/data/' + str(prefix + '_float_data.csv'), header=1)
	data_float_psi.insert(2, "Float Pressure (psia)", [x+14.7 for x in data_float_psi["Float Pressure (psig)"]])
	b, a = signal.butter(5, 0.5)
	data_float_psi_filtered = signal.filtfilt(b, a, data_float_psi['Float Pressure (psig)'])
	data_float_psi = pd.concat([data_float_psi, pd.DataFrame(data_float_psi_filtered, columns=['Float Pressure Filtered (psig)'])], axis=1)
	tmax_float_psi = round(data_float_psi['Time (s)'].iloc[-1], 1)
	# n_float_psi_resampled = int(tmax_float_psi*10 + 1)
	n_all_resampled =int(tmax_float_psi*10 + 1)
	tnew_float_psi = np.linspace(0, tmax_float_psi, n_all_resampled)
	data_float_psi_resampled = signal.resample(data_float_psi['Float Pressure Filtered (psig)'], n_all_resampled)
	# data_float_psi_resampled = signal.resample_poly(data_float_psi['Float Pressure (psig)'], 10, 9)
	tnew_float_psi = np.linspace(0, tmax_float_psi, data_float_psi_resampled.size)
	data_float_psi = pd.concat([data_float_psi, pd.DataFrame(np.transpose(np.array([tnew_float_psi, data_float_psi_resampled])), columns=['Time Resampled (s)', 'Float Pressure Resampled (psig)'])], axis=1)


	data_prop_psi = pd.read_csv('/home/josh/remoteProp/data/' + str(prefix + '_prop_data.csv'), header=1)
	data_prop_psi.insert(2, "Prop Pressure (psia)", [x+14.7 for x in data_prop_psi["Prop Pressure (psig)"]])
	data_prop_psi_filtered = signal.filtfilt(b, a, data_prop_psi['Prop Pressure (psig)'])
	data_prop_psi = pd.concat([data_prop_psi, pd.DataFrame(data_prop_psi_filtered, columns=['Prop Pressure Filtered (psig)'])], axis=1)
	tmax_prop_psi = round(data_prop_psi['Time (s)'].iloc[-1], 1)
	# n_prop_psi_resampled = int(tmax_prop_psi*10 + 1)
	tnew_prop_psi = np.linspace(0, tmax_prop_psi, n_all_resampled)
	data_prop_psi_resampled = signal.resample(data_prop_psi['Prop Pressure Filtered (psig)'], n_all_resampled)
	# data_prop_psi_resampled = signal.resample_poly(data_prop_psi['Prop Pressure (psig)'], 10, 9)
	tnew_prop_psi = np.linspace(0, tmax_prop_psi, data_prop_psi_resampled.size)
	data_prop_psi = pd.concat([data_prop_psi, pd.DataFrame(np.transpose(np.array([tnew_prop_psi, data_prop_psi_resampled])), columns=['Time Resampled (s)', 'Prop Pressure Resampled (psig)'])], axis=1)


	data_weight = pd.read_csv('/home/josh/remoteProp/data/' + str(prefix + '_loadcell_data.csv'), header=1)
	data_weight.insert(2, "Thrust (mN)", [x*9.81 for x in data_weight["Weight (?)"]])
	b, a = signal.butter(3, 0.85)
	data_weight_filtered = signal.filtfilt(b, a, data_weight['Thrust (mN)'])
	data_weight = pd.concat([data_weight, pd.DataFrame(data_weight_filtered, columns=['Thrust Filtered (mN)'])], axis=1)
	tmax_weight = round(data_weight['Time (s)'].iloc[-1], 1)
	# n_weight_resampled = int(tmax_weight*10 + 1)
	tnew_weight = np.linspace(0, tmax_weight, n_all_resampled)
	data_weight_resampled = signal.resample(data_weight['Thrust Filtered (mN)'], n_all_resampled)
	# data_weight_resampled = signal.resample_poly(data_weight['Thrust (mN)'], 10, 9)
	tnew_weight = np.linspace(0, tmax_weight, data_weight_resampled.size)
	data_weight = pd.concat([data_weight, pd.DataFrame(np.transpose(np.array([tnew_weight, data_weight_resampled])), columns=['Time Resampled (s)', 'Thrust Resampled (mN)'])], axis=1)

	thrust_corrected = []
	for i, x in enumerate(data_weight['Thrust Resampled (mN)']):
		thrust_corrected.append(x + (100 - data_float_psi['Float Pressure Filtered (psig)'][i])*(37.6/100))
	data_weight = pd.concat([data_weight, pd.DataFrame(thrust_corrected, columns=['Thrust Corrected (mN)'])], axis=1)


	# data_temp = pd.read_csv('/home/josh/remoteProp/data/' + str(prefix + '_temp_data.csv'), header=1)
	# data_temp.insert(2, "Temperature (K)", [x+273.15 for x in data_temp["Exit Temperature (Celsius)"]])
	# data_temp_filtered = signal.filtfilt(b, a, data_temp['Temperature (K)'])
	# data_temp = pd.concat([data_temp, pd.DataFrame(data_temp_filtered, columns=['Temperature Filtered (K)'])], axis=1)
	# tmax_temp = round(data_temp['Time (s)'].iloc[-1], 1)
	# # n_temp_resampled = int(tmax_temp*10 + 1)
	# tnew_temp = np.linspace(0, tmax_temp, n_all_resampled)
	# data_temp_resampled = signal.resample(data_temp['Temperature Filtered (K)'], n_all_resampled)
	# # data_temp_resampled = signal.resample_poly(data_temp['Temperature (K)'], 10, 9)
	# tnew_temp = np.linspace(0, tmax_temp, data_temp_resampled.size)
	# data_temp = pd.concat([data_temp, pd.DataFrame(np.transpose(np.array([tnew_temp, data_temp_resampled])), columns=['Time Resampled (s)', 'Temperature Resampled (K)'])], axis=1)


	# return [data_float_psi, data_prop_psi, data_weight, data_temp]

	return [data_float_psi, data_prop_psi, data_weight]

# Calculate water hammer pressure spike
dia_inner = 2.8/1000  # mm/1000 = m
area = math.pi*(dia_inner**2)/4  # m^2

temp_0 = -40  # C
temp_1 = -40.0385
m_dot = 0.76  # g/s
psia_0 = 113.2  # psia
psia_1 = 103  # psia

# ---------------------------------------------------------------
k = 1.298
R = 8.314/0.044041  # J/kg-K
temp_0 += 273.15
temp_1 += 273.15
m_dot /= 1000
pa_0 = psia_0*6894.76  # psia*6894.76 = Pascals (N/m^2)
pa_1 = psia_1*6894.76  # psia*6894.76 = Pascals (N/m^2)
rho_0 = pa_0/(R*temp_0)
rho_1 = pa_1/(R*temp_1)

a = math.sqrt(k*R*temp_0)

# Bernoulli says...
delta_v_bern = math.sqrt( 2*( (k*R/(k-1))*(temp_0 - temp_1) ) )
delta_v = m_dot/(rho_1*area)
hammer_pres = rho_0*a*delta_v  # Pascals (N/m^2)

print('')
print('dV Bernoulli: ' + str(delta_v_bern) + ' m/s')
print('dV from mdot: ' + str(delta_v) + ' m/s')
print('Speed of sound 0: ' + str(a) + ' m/s')
print('')
print("Hammer pressure (delta-P): " + str(hammer_pres/6894.76) + " psid")
print("Max pressure (psia): " + str((hammer_pres + pa_0)/6894.76) + " psia")
print('')
print('rho_0: ' + str(rho_0) + ' kg/m^3')
print('rho_1: ' + str(rho_1) + ' kg/m^3')
print('')

# Plot Data
# test_nos = ["20191130_120249"]  # Plenum discharge, bad data
# test_nos = ["20191130_120842"]  # Plenum discharge, bad data
# test_nos = ["20191130_122101"]  # Plenum discharge, bad data
# test_nos = ["20191130_124508", "20191130_124745"]  # Plenum discharge, bad data

# ---=== Plenum Discharge, 115 psia ===---
# test_nos = ["20191130_130935", "20191130_130954", "20191130_131419", "20191130_131515", "20191130_131607", "20191130_131624", "20191130_131644"]

# ---=== Steady-state Thrust Tests ===---
test_nos = ["20191204_143449", "20191204_143538", "20191204_143636"]  # 35 psia steady-state
test_nos = ["20191204_143735", "20191204_143817", "20191204_143931"]  # 55 psia steady-state
test_nos = ["20191204_144142", "20191204_144207", "20191204_144256"]  # 75 psia steady-state
test_nos = ["20191204_144715", "20191204_144754", "20191204_144840"]  # 95 psia steady-state
test_nos = ["20191204_145400", "20191204_145419", "20191204_145442"]  # 115 psia steady-state
test_nos = ["20191205_191138",  "20191205_191210",  "20191205_191332",  "20191205_191402",  "20191205_191433"]  # 115 psia steady-state

# test_nos = ["20191223_183908", "20191223_183945", "20191223_183832", "20191223_183725", "20191223_183658"]  # Plenum discharge
# test_nos = ["20191223_183908"]  # Plenum discharge

colors = []
fig1, ax1 = plt.subplots(figsize=(6.5, 4), dpi=90)
ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

for trial, test in enumerate(test_nos):
	test_data = all_data(test)

	ax1.plot(test_data[0]['Time (s)'], test_data[0]['Float Pressure (psig)'], 
		# color='#2ca02c', 
		label='Trial #' + str(trial+1) + ' Pres Raw.', 
		# marker='+', 
		markersize='10', 
		linestyle='--', 
		linewidth='1')
	ax1.plot(test_data[0]['Time Resampled (s)'], test_data[0]['Float Pressure Resampled (psig)'], 
		# color='#1f77b4', 
		label='Trial #' + str(trial+1) + ' Pres. Resampled', 
		# marker='x', 
		markersize='10', 
		linestyle='-', 
		linewidth='1')

	# ax2.plot(test_data[2]['Time (s)'], test_data[2]['Thrust (mN)'], 
	# 	# color='#ff7f0e', 
	# 	label='Trial #' + str(trial+1) + ' Thrust', 
	# 	marker='o',
	# 	fillstyle='none', 
	# 	linestyle='-', 
	# 	linewidth='1')		
	# ax2.plot(test_data[2]['Time Resampled (s)'], test_data[2]['Thrust Resampled (mN)'], 
	# 	# color='#ff7f0e', 
	# 	label='Trial #' + str(trial+1) + ' Thrust', 
	# 	marker='o',
	# 	fillstyle='none', 
	# 	linestyle='-', 
	# 	linewidth='1')

	# Blue: #1f77b4
	# Orange: #ff7f0e
	# Green: #2ca02c
ax1.set_xlabel('Time (s)', color='#413839')
ax1.set_ylabel('Pressure (psig)', color='#413839')
ax2.set_ylabel('Thrust (mN)', color='#413839')

ax1.set_xlim([0, 11])

ax1.tick_params(colors='#413839')
ax1.grid(which='major', axis='both', linestyle='--')
box = ax1.get_position()
ax1.set_position([box.x0, box.y0 + box.height*0.1, box.width, box.height*0.9])

fig1.legend(loc='center', bbox_to_anchor=(0.5, 0.03), ncol=10, frameon=False )
# plt.title('Single Plenum Discharge Pressure and Thrust ({1} mm Nozzle)'.format(test_nos,0.6), y=1.03, color='#413839')

# plt.show()

plt.close('all')
linewidth = 2
fontsize = 12
fig2, ax1 = plt.subplots(figsize=(5.75, 4), dpi=200)

td = []

for trial, test in enumerate(test_nos):
	test_data = all_data(test)[0]  # All the Float Pressure data
	test_data["Trial"] = trial
	td.append(test_data)

td = pd.concat(td)
td['Time (s)'] = td['Time (s)'].round(1)
sns.lineplot(ax=ax1, x=td['Time (s)'], y=td['Float Pressure (psia)'], data=td, label='Pressure (psia)', linewidth=linewidth)

# plt.legend(loc='upper left', bbox_to_anchor=(0.68, 1), ncol=1, frameon=False )
plt.legend(loc='center', bbox_to_anchor=(0.2, -0.25), ncol=1, frameon=False )

plt.tick_params(colors='#413839')
plt.grid(which='major', axis='both', linestyle='--')


ax2 = plt.twinx()  # instantiate a second axes that shares the same x-axis
td = []

for trial, test in enumerate(test_nos):
	test_data = all_data(test)[2]  # All the Thrust data
	test_data["Trial"] = trial
	td.append(test_data)

td = pd.concat(td)
td['Time (s)'] = td['Time (s)'].round(1)
sns.lineplot(x=td['Time (s)'], y=td['Thrust (mN)'], color='Orange', data=td, ax=ax2, label='Thrust (mN)', linewidth=linewidth)


box = ax2.get_position()
ax1.set_position([box.x0 + box.width * 0.03, box.y0 + box.height * 0.15, box.width * 0.95, box.height * 0.95])
ax1.set_xlabel('Time, s', color='#413839', fontsize=fontsize)
ax1.set_ylim
ax1.set_ylabel('Pressure, psia', color='#413839', fontsize=fontsize)
ax2.set_ylabel('Thrust, mN', color='#413839', fontsize=fontsize)

# plt.legend(loc='upper left', bbox_to_anchor=(0.68, 0.95), ncol=1, frameon=False )
plt.legend(loc='center', bbox_to_anchor=(0.7, -0.25), ncol=1, frameon=False )

# plt.title('Supply Pressure and Thrust (0.6 mm Nozzle)')
# plt.savefig('/mnt/d/OneDrive - UC Davis/HRVIP/Writing/AIAA SciTech 2019 Paper/Images/Sim Results/image.png')

plt.show()