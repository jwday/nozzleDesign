import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

def all_data(prefix):
	data_float_psi = pd.read_csv(str(prefix + '_float_data.csv'))
	data_float_psi.insert(2, "Float Pressure (psia)", [x+14.7 for x in data_float_psi["Float Pressure (psig)"]])

	data_prop_psi = pd.read_csv(str(prefix + '_prop_data.csv'))
	data_prop_psi.insert(2, "Prop Pressure (psia)", [x+14.7 for x in data_prop_psi["Prop Pressure (psig)"]])

	data_weight = pd.read_csv(str(prefix + '_loadcell_data.csv'))
	data_weight.insert(2, "Thrust (mN)", [x*9.81 for x in data_weight["Weight (?)"]])

	data_temp = pd.read_csv(str(prefix + '_temp_data.csv'))
	data_temp.insert(2, "Temperature (K)", [x+273.15 for x in data_temp["Exit Temperature (Celsius)"]])

	return [data_float_psi, data_prop_psi, data_weight, data_temp]



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