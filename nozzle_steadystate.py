# Steady-state Isentropic Nozzle Flow: Pressure sweep for constant temperature

from nozzle_code import nozzle
import matplotlib.pyplot as plt
from scipy import stats
import math
import pandas as pd
import numpy as np
import sys

# Gas conditions
P_amb = 14.7  # Ambient Pressure (psia)
T_t = -20 + 273.15  # Total Temperature (K)

# Nozzle geometry
d_star = 0.6  # Throat diameter (mm)
half_angle = 10  # Conical Nozzle expansion angle (degrees)
expansion_ratio = 1.3225
# # PSI  ----- Exp. Ratio
# # 114.7 --- 1.8048
# # 65 --- 1.3225
# # 

# Gas Choices: R236fa, R134a, N2, CO2, air
gas_type = 'CO2'


list_of_P_ts = [x+14.7 for x in list(range(101))]
list_of_T_ts = [T_t]*101
list_of_mdots = []
list_of_thrusts = []


for P_t in list_of_P_ts:
    m_dot, M_crit_sub, M_crit_sup, PR_crit_sub, PR_crit_sup, PR_exit_shock, M_exit_behindshock, M_exit, P_exit, v_exit, F, P_star, T_star, rho_star, T_exit, rho_exit = nozzle(P_t, T_t, P_amb, d_star, expansion_ratio, half_angle, gas_type)
    
    list_of_mdots.append(m_dot*1000)
    list_of_thrusts.append(F)



# Import experimental data
data1 = pd.read_csv('10312019_test5.csv')
data1.insert(1, "Pressure (psia)", [x+14.7 for x in data1["Pressure (psig)"]], True)
data1.insert(5, "Beta min init", [round(x,4) for x in (data1["Ma init (g)"]-0.3)/(data1["Mb init (g)"]+0.03)])
data1.insert(6, "Beta max init", [round(x,4) for x in (data1["Ma init (g)"]+0.3)/(data1["Mb init (g)"]-0.03)])

data1.insert(9, "Beta min final", [round(x,4) for x in (data1["Ma final (g)"]-0.3)/(data1["Mb final (g)"]+0.03)])
data1.insert(10, "Beta max final", [round(x,4) for x in (data1["Ma final (g)"]+0.3)/(data1["Mb final (g)"]-0.03)])

data1.insert(11, "dM min", [round(x,2) for x in (data1["Beta min init"]+1)*data1["Mb init (g)"] + 
												(data1["Beta min init"]-1)*0.03 -
												(data1["Beta max final"]+1)*data1["Mb final (g)"] +
												(data1["Beta max final"]-1)*0.03])

data1.insert(12, "dM max", [round(x,2) for x in (data1["Beta max init"]+1)*data1["Mb init (g)"] - 
												(data1["Beta max init"]-1)*0.03 -
												(data1["Beta min final"]+1)*data1["Mb final (g)"] -
												(data1["Beta min final"]-1)*0.03])

data1.insert(13, "dM/dt min", [ round(x,3) for x in data1["dM min"]/data1["Time (s)"] ])
data1.insert(14, "dM/dt max", [ round(x,3) for x in data1["dM max"]/data1["Time (s)"] ])
data1.insert(15, "dM/dt nom", [ round(x,2) for x in (data1["Ma init (g)"] + data1["Mb init (g)"] - data1["Ma final (g)"] - data1["Mb final (g)"])/data1["Time (s)"] ])
data1.insert(16, "dM/dt err", [ x for x in data1["dM/dt max"]-data1["dM/dt min"] ])



# Setup theoretical vs. experimental plot
    # Blue: #1f77b4
    # Orange: #ff7f0e
    # Green: #2ca02c

fig1, ax1 = plt.subplots(figsize=(8.5, 5), dpi=90)
ax1.set_xlabel('Pressure (psia)', color='#413839')
ax1.set_ylabel('Mass Flow Rate (g/s)', color='#413839')

ax1.plot(list_of_P_ts, list_of_mdots, color='#1f77b4', label='isentropic')
ax1.errorbar(data1["Pressure (psia)"], data1["dM/dt nom"], xerr=2, yerr=data1["dM/dt err"]/2, color='#2ca02c', label='experimental', linestyle='none', marker='x')

slope, intercept, r_value, p_value, std_err = stats.linregress(data1["Pressure (psia)"], data1["dM/dt nom"])
ax1.plot(data1["Pressure (psia)"], slope*data1["Pressure (psia)"]+intercept, color='#2ca02c', linestyle='--', label='_nolegend_')

ax1.tick_params(colors='#413839')
ax1.grid(which='major', axis='both', linestyle='--')


ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

ax2.set_ylabel('Temperature (K)', color='#413839')
ax2.plot(list_of_P_ts, list_of_T_ts, color='#ff7f0e', label='Critical Temperature (K)')

box = ax1.get_position()
ax1.set_position([box.x0, box.y0 + box.height*0.1, box.width, box.height*0.9])

fig1.legend(['Inviscid', 'Experimental', 'Temperature'], loc='center', bbox_to_anchor=(0.5, 0.03), ncol=3, frameon=False )
plt.title('Mass Flow Rate Comparison ({} mm)'.format(d_star), y=1.03, color='#413839')

# ---- Plot Everything
plt.show()