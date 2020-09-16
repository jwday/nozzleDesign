from scipy import stats
from scipy import signal
import math
import pandas as pd
import numpy as np
import sys
import seaborn as sns
import os

def all_data(prefix, mult_by_g=False):
	data_float_psi = pd.read_csv(str(prefix + '_float_data.csv'), header=1)  # Init the dataframe
	data_float_psi.insert(2, "Float Pressure (psia)", [x+14.7 for x in data_float_psi["Float Pressure (psig)"]])  # Insert psia column
	b, a = signal.butter(5, 0.5)  # Create a filter
	data_float_psi_filtered = signal.filtfilt(b, a, data_float_psi['Float Pressure (psig)'])  # Apply the filter to the psig data
	data_float_psi = pd.concat([data_float_psi, pd.DataFrame(data_float_psi_filtered, columns=['Float Pressure Filtered (psig)'])], axis=1)  # Insert filtered psig column
	tmax_float_psi = round(data_float_psi['Time (s)'].iloc[-1], 1)  # Determine a max time for the next step
	n_all_resampled =int(tmax_float_psi*10 + 1)  # Use the max time to determine number of resampled data points to create
	tnew_float_psi = np.linspace(0, tmax_float_psi, n_all_resampled)  # Create a new list of equally spaced times based on the number of desired resampled data points
	data_float_psi_resampled = signal.resample(data_float_psi['Float Pressure Filtered (psig)'], n_all_resampled)  # Resample the filtered psig data
	tnew_float_psi = np.linspace(0, tmax_float_psi, data_float_psi_resampled.size)  # Shouldn't this be the same? Yeah it should be
	data_float_psi = pd.concat([data_float_psi, pd.DataFrame(np.transpose(np.array([tnew_float_psi, data_float_psi_resampled])), columns=['Time Resampled (s)', 'Float Pressure Resampled (psig)'])], axis=1)
	data_float_psi.insert(3, 'Float Pressure (kPa)', data_float_psi['Float Pressure (psia)'].multiply(6.8948))


	data_prop_psi = pd.read_csv(str(prefix + '_prop_data.csv'), header=1)
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
	data_prop_psi.insert(3, 'Prop Pressure (kPa)', data_prop_psi['Prop Pressure (psia)'].multiply(6.8948))


	if mult_by_g: 
		g = 9.81
	else:
		g = 1
	data_weight = pd.read_csv(str(prefix + '_loadcell_data.csv'), header=1)
	data_weight.insert(2, "Thrust (mN)", [x*g for x in data_weight["Weight (?)"]])
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
	thrust_correction_factor = data_weight['Thrust Resampled (mN)'].dropna().tail().mean()
	for i, x in enumerate(data_weight['Thrust Resampled (mN)'].dropna()):
		thrust_corrected.append(x - (100 - data_float_psi['Float Pressure Filtered (psig)'][i])*(thrust_correction_factor/100))
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

	
# def all_data(prefix):
# 	data_float_psi = pd.read_csv(str(prefix + '_float_data.csv'))
# 	data_float_psi.insert(2, "Float Pressure (psia)", [x+14.7 for x in data_float_psi["Float Pressure (psig)"]])

# 	data_prop_psi = pd.read_csv(str(prefix + '_prop_data.csv'))
# 	data_prop_psi.insert(2, "Prop Pressure (psia)", [x+14.7 for x in data_prop_psi["Prop Pressure (psig)"]])

# 	data_weight = pd.read_csv(str(prefix + '_loadcell_data.csv'))
# 	data_weight.insert(2, "Thrust (mN)", [x*9.81 for x in data_weight["Weight (?)"]])

# 	data_temp = pd.read_csv(str(prefix + '_temp_data.csv'))
# 	data_temp.insert(2, "Temperature (K)", [x+273.15 for x in data_temp["Exit Temperature (Celsius)"]])

# 	return [data_float_psi, data_prop_psi, data_weight, data_temp]



def massflow_data_test11(filename, sys_mass, sys_psia):
	data = pd.read_csv(str(filename))
	data.insert(2, "P init (psia)", [x+14.7 for x in data["P init (psig)"]], True)
	
	sys_char_max = (sys_mass + 0.03)/(sys_psia - 2)  # sys_char = (m1 / P1)
	sys_char_min = (sys_mass - 0.03)/(sys_psia + 2)

	data.insert(7, "sys mass max (g)", [ round(x,2) for x in (data["P init (psia)"]+2)*sys_char_max ])
	data.insert(8, "sys mass min (g)", [ round(x,2) for x in (data["P init (psia)"]-2)*sys_char_min ])

	data.insert(9, "dM max", [ round(x,2) for x in (data["M init (g)"]+0.03) - (data["M final (g)"]-0.03) - data["sys mass min (g)"] ])
	data.insert(10, "dM min", [ round(x,2) for x in (data["M init (g)"]-0.03) - (data["M final (g)"]+0.03) - data["sys mass max (g)"] ])

	data.insert(11, "dM/dt min", [ round(x,3) for x in data["dM min"] / data["Time (s)"] ])
	data.insert(12, "dM/dt max", [ round(x,3) for x in data["dM max"] / data["Time (s)"] ])
	data.insert(13, "dM/dt nom", [ round(x,2) for x in (data["M init (g)"] - data["M final (g)"]) / data["Time (s)"] ])
	data.insert(14, "dM/dt err", [ x for x in data["dM/dt max"]-data["dM/dt min"] ])

	data.insert(15, "Pressure (psia)", [ round(x,2) for x in (data["P init (psig)"] + data["P fin (psig)"])/2 + 14.7 ])

	return data



def massflow_data_singlescale(filename):
	data = pd.read_csv(str(filename))
	data.insert(1, "Pressure (psia)", [x+14.7 for x in data["Pressure (psig)"]], True)

	data.insert(4, "dM max", [ round(x,2) for x in (data["M init (g)"]+0.3) - (data["M final (g)"]-0.3) ])
	data.insert(5, "dM min", [ round(x,2) for x in (data["M init (g)"]-0.3) - (data["M final (g)"]+0.3) ])

	data.insert(7, "dM/dt min", [ round(x,3) for x in data["dM min"] / data["Time (s)"] ])
	data.insert(8, "dM/dt max", [ round(x,3) for x in data["dM max"] / data["Time (s)"] ])
	data.insert(9, "dM/dt nom", [ round(x,2) for x in (data["M init (g)"] - data["M final (g)"]) / data["Time (s)"] ])
	data.insert(10, "dM/dt err", [ x for x in data["dM/dt max"]-data["dM/dt min"] ])

	return data



def massflow_data(filename):
	data = pd.read_csv(str(filename))
	data.insert(1, "Pressure (psia)", [ x+14.7 for x in data["P init (psig)"] ], True)
	data.insert(5, "Beta min init", [ round(x,4) for x in (data["Ma init (g)"]-0.3) / (data["Mb init (g)"]+0.03) ])
	data.insert(6, "Beta max init", [ round(x,4) for x in (data["Ma init (g)"]+0.3) / (data["Mb init (g)"]-0.03) ])

	data.insert(9, "Beta min final", [ round(x,4) for x in (data["Ma final (g)"]-0.3) / (data["Mb final (g)"]+0.03) ])
	data.insert(10, "Beta max final", [ round(x,4) for x in (data["Ma final (g)"]+0.3) / (data["Mb final (g)"]-0.03) ])

	data.insert(11, "dM min", [ round(x,2) for x in (data["Beta min init"]+1)*data["Mb init (g)"] + 
													(data["Beta min init"]-1)*0.03 -
													(data["Beta max final"]+1)*data["Mb final (g)"] +
													(data["Beta max final"]-1)*0.03 ])

	data.insert(12, "dM max", [ round(x,2) for x in (data["Beta max init"]+1)*data["Mb init (g)"] - 
													(data["Beta max init"]-1)*0.03 -
													(data["Beta min final"]+1)*data["Mb final (g)"] -
													(data["Beta min final"]-1)*0.03 ])

	data.insert(13, "dM/dt min", [ round(x,3) for x in data["dM min"]/data["Time (s)"] ])
	data.insert(14, "dM/dt max", [ round(x,3) for x in data["dM max"]/data["Time (s)"] ])
	data.insert(15, "dM/dt nom", [ round(x,2) for x in (data["Ma init (g)"] + data["Mb init (g)"] - data["Ma final (g)"] - data["Mb final (g)"])/data["Time (s)"] ])
	data.insert(16, "dM/dt err", [ x for x in data["dM/dt max"]-data["dM/dt min"] ])	

	return data



def thrust_data(filename, mass_offset):
	data = pd.read_csv(str(filename))
	data.insert(1, "Pressure (psia)", [x+14.7 for x in data["Pressure (psig)"]], True)

	time_offset = data["Timestamp (ms)"][0]

	data.insert(3, "Time corrected (ms)", [x-time_offset for x in data["Timestamp (ms)"]])
	data.insert(5, "Mass corrected (g)", [x-mass_offset for x in data["Mass (g)"]])
	data.insert(6, "Thrust (mN)", [round(9.80665*x,3) for x in data["Mass corrected (g)"]])

	return data