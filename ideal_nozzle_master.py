# Nozzle master control
# Given a nozzle geometry and initial conditions, this code will sweep through a range of stagnation pressures and output the exit conditions
# It's interesting to note that the critical conditions are dependent ONLY on geometry and not stagnation conditions
from nozzle_code import nozzle
import matplotlib.pyplot as plt
from scipy import stats
import math
import pandas as pd
import numpy as np
import sys
import seaborn as sns
from data_handling_funcs import *


## ---- OPTIONS --------------------------------------------------------------------

# Gas initial conditions
P_amb = 14.7  # Ambient Pressure (psia)
list_of_P_ts = list(range(100))
list_of_Tts = [218, 293]



# Nozzle geometry
d_star = 0.6  # Throat diameter (mm)
expansion_ratio = 1.3225  # 1.8048 for ideal expansion at 114.7 psi supply, 2.2447 for 164.7, 1.3225 for 0.2mm and 64.7 psi, 1.1235 for 0.3 and 44.7 psi
# # PSI  ----- Exp. Ratio
# # 114.7 --- 1.8048
# # 65 --- 1.3225
# # 
half_angle = 10  # Conical Nozzle expansion angle (degrees)


# Gas Choices: R236fa, R134a, N2, CO2, air
gas_type = 'CO2'
k = 1.289
R = 8.314/0.04401  # Specific gas constant (J/kg-K)

data_dict = {}

## ---- DO THE THING ----------------------------------------------------------------
for T in list_of_Tts:
    list_of_mdots = []
    list_of_thrusts = []
    list_of_T_exits = []

    for P in list_of_P_ts:
        m_dot, M_crit_sub, M_crit_sup, PR_crit_sub, PR_crit_sup, PR_exit_shock, M_exit_behindshock, M_exit, P_exit, v_exit, F, P_star, T_star, rho_star, T_exit, rho_exit = nozzle(P, T, P_amb, d_star, expansion_ratio, half_angle, gas_type)

        list_of_mdots.append(m_dot*1000)
        list_of_thrusts.append(F)
        list_of_T_exits.append(T_exit)

    single_iter = pd.DataFrame(data={'mdot': list_of_mdots, 'thrust': list_of_thrusts})
    label = str(T)
    data_dict[label] = single_iter




# Plot Mass Flow Rate (using 16g cartridges)
massflow_data11 = massflow_data_test11('Test Data/11132019_test11.csv', 0.0008809746, 0.00079844206)    

fig1, ax1 = plt.subplots(figsize=(6.5, 4), dpi=90)
ax1.set_xlabel('Pressure (psia)', color='#413839')
ax1.set_ylabel('Mass Flow Rate (g/s)', color='#413839')

colors = ['#1f77b4', '#ff7f0e']
    # Blue: #1f77b4
    # Green: #2ca02c
    # Orange: #ff7f0e
	# Gray: #808080
linestyles = ['-', '--']

for i, T in enumerate(list_of_Tts):
    label = 'Inviscid (' + str(T-273) + ' C)'
    ax1.plot(list_of_P_ts, data_dict[str(T)]['mdot'], color=colors[i], label=label, linestyle=linestyles[i])

massflow_data = [massflow_data11]

for name in massflow_data:
	ax1.errorbar(name["P init (psia)"], name["dM/dt nom"], xerr=2, yerr=name["dM/dt err"]/2, color='#2ca02c', label='Experimental', linestyle='none', marker='x')

ax1.tick_params(colors='#413839')
ax1.grid(which='major', axis='both', linestyle='--')
box = ax1.get_position()
ax1.set_position([box.x0, box.y0 + box.height*0.1, box.width, box.height*0.9])

# legend((line1, line2, line3), ('label1', 'label2', 'label3'))
fig1.legend(loc='center', bbox_to_anchor=(0.5, 0.03), ncol=3, frameon=False)
plt.title('Mass Flow Rate ({} mm Nozzle)'.format(d_star), y=1.03, color='#413839')

ax1.set_xlim([10, 110])
ax1.set_ylim([0, 0.8])



# Thrust vs. Pressure (Inviscid and Experimental)
# test8_trial1 = thrust_data('Test Data/11062019_thrust_test1.csv', 97.56)
# thrust_data2 = thrust_data('11062019_thrust_test2.csv', 97.61)
# thrust_data3 = thrust_data('11062019_thrust_test3.csv', 97.51)
test8_trial7 = thrust_data('Test Data/11062019_test8_thrust_trial7.csv', 98.18)
test8_trial10 = thrust_data('Test Data/11062019_test8_thrust_trial10.csv', 97.71)
test8_trial11 = thrust_data('Test Data/11062019_test8_thrust_trial11.csv', 97.62)
# test8_trial12 = thrust_data('Test Data/11062019_test8_thrust_trial12.csv', 97.66)

# test9_trial1 = thrust_data('Test Data/11072019_test9_thrust_trial1.csv', 144.33)
# test9_trial5 = thrust_data('Test Data/11072019_test9_thrust_trial5.csv', 144.65)
test9_trial9 = thrust_data('Test Data/11072019_test9_thrust_trial9.csv', 145.07)

data_points = [test8_trial7, test8_trial10, test8_trial11, test9_trial9]

fig8, ax1 = plt.subplots(figsize=(6.5, 4), dpi=90)
	# Blue: #1f77b4 (Inviscid)
	# Green: #2ca02c
	# Orange: #ff7f0e
	# Gray: #808080
ax1.set_xlabel('Pressure (psia)', color='#413839')
ax1.set_ylabel('Thrust (mN)', color='#413839')
for i, T in enumerate(list_of_Tts):
    label = 'Inviscid (' + str(T-273) + ' C)'
    ax1.plot(list_of_P_ts, data_dict[str(T)]['thrust']*1000, color=colors[i], label=label, linestyle=linestyles[i])

# for name in data_points:
# 	ax1.plot(name["Pressure (psia)"], name["Thrust (mN)"], color='#2ca02c', label='Experimental', linestyle='none', marker='x')

all_thrust_data = pd.concat(data_points)
ax1.plot(all_thrust_data["Pressure (psia)"], all_thrust_data["Thrust (mN)"], color='#2ca02c', label='Experimental', linestyle='none', marker='x')

slope, intercept, r_value, p_value, std_err = stats.linregress(all_thrust_data["Pressure (psia)"], all_thrust_data["Thrust (mN)"])
ax1.plot(all_thrust_data["Pressure (psia)"], slope*all_thrust_data["Pressure (psia)"]+intercept, color='#2ca02c', linestyle='--', label='_nolegend_')

ax1.tick_params(colors='#413839')
ax1.grid(which='major', axis='both', linestyle='--')

box = ax1.get_position()
ax1.set_position([box.x0, box.y0 + box.height*0.1, box.width, box.height*0.9])

fig8.legend(loc='center', bbox_to_anchor=(0.5, 0.03), ncol=3, frameon=False )
# plt.title('Thrust Comparison ({} mm)'.format(d_star), y=1.03, color='#413839')

ax1.set_xlim([10, 110])



# ---- Plot Everything
plt.show()
