# plot_ss_pres_and_mdot.py
# Use simplied nozzle code (nozzle_code_2)

from nozzle_code_2 import nozzle
import matplotlib.pyplot as plt
from scipy import stats
import math
import pandas as pd
import numpy as np
import sys
from data_handling_funcs import *

## ==================================================================================
## ---- USER OPTIONS ----------------------------------------------------------------
## ==================================================================================

P_amb = 14.7 * 6894.76  			# Ambient Pressure, units of Pa (psia * 6894.76)
d_star = 0.6 / 1000  			# Nozzle throat diameter, units of m (mm / 1000)
expansion_ratio = 1.17			# Nozzle expansion ratio (Exit Area / Throat Area)
half_angle = 10  				# (Conical) Nozzle expansion angle (degrees)
gas_type = 'CO2'				# Gas Choices: R236fa, R134a, N2, CO2, H2, air



k = 1.289						# Specific heat ratio (NOT thermal conductivity)
R = 8.314/0.04401  				# Specific gas constant (J/kg-K)


## ==================================================================================
## ---- MACH-AREA RELATION ----------------------------------------------------------
## ==================================================================================

# Nozzle Geometry
A_star = np.pi*(d_star**2)/4  					# Throat area
A_exit = A_star*expansion_ratio  				# Exit area
d_exit = np.sqrt(4*A_exit/np.pi)  				# Exit diameter (m)

# Common isentropic nozzle relationships rewritten in compact terms
P = (k+1)/2
L = (k-1)/2
W = k/(k-1)
Q = P/L  # aka (k+1)/(k-1)
S = (A_exit/A_star)**2
Y = np.sqrt(k/R)

def machAreaRelation(X, S):
	f = (1 + L*X)**Q - S*X*(P**Q)  # Where 'X' is M^2 and S is AREA RATIO SQUARED (so that this same function can be used for different area ratios i.e. along the length of the nozzle)
	return f

x0 = np.array([0.00001, 20])	# List of M^2 to solve for roots over
sol0 = opti.fsolve(machAreaRelation, x0, args=S, maxfev=100000, full_output=False, xtol=0.000000001)	# Solve for M^2 given S as an extra argument

M_crit_sub = np.sqrt(sol0[0])  # Subsonic critical mach no.
M_crit_sup = np.sqrt(sol0[1])  # Supersonic critical mach no.




## ==================================================================================
## ---- INIT DATA LISTS -------------------------------------------------------------
## ==================================================================================

list_of_P_ts = [(x+14.7)*6894.76 for x in list(range(101))]
list_of_T_ts = [220, 260, 293]
list_of_mdots = []
list_of_thrusts = []



## ==================================================================================
## ---- EXECUTE LOOP ----------------------------------------------------------------
## ==================================================================================

for P_t in list_of_P_ts:
    # m_dot, M_crit_sub, M_crit_sup, PR_crit_sub, PR_crit_sup, PR_exit_shock, M_exit_behindshock, M_exit, P_exit, v_exit, F, P_star, T_star, rho_star, Re_star, T_exit, rho_exit = nozzle(P_t, list_of_T_ts[-1], P_amb, d_star, expansion_ratio, half_angle, gas_type)
	P_star, T_star, rho_star, Re_star, M_star, v_star, P_exit, T_exit, rho_exit, M_exit, v_exit, c_exit, m_dot, F, CF, flow_regime, area_ratio_at_shock, F_mdotv, F_pdiff = nozzle(k, R, M_crit_sub, M_crit_sup, P_t, T_t_inlet[-1], rho_t_inlet[-1], P_amb, d_star, expansion_ratio, half_angle, gas_type, visc_func, r_from_PT_gas_func)
    
    list_of_mdots.append(m_dot*1000)
    list_of_thrusts.append(F)

sim_results = pd.DataFrame(data={'Pressure (psig)': [(x/6894.76)-14.7 for x in list_of_P_ts], 'Thrust (mN)': [x*1000 for x in list_of_thrusts], 'mdot (g/s)': list_of_mdots})


# fig1, axs = plt.subplots(1, sharex=True, dpi=150, figsize=(6, 4))
# fig1.suptitle('Steady-State Thrust vs. Pressure', y=0.98, fontsize=12)
# axs.set_title(r'CO$_2$, Nozzle $\varnothing$0.6 mm, $\lambda$=1.17', fontsize=8, color='dimgray', y=1)

# sns.lineplot(ax=axs, x="Pressure (psig)", y="Thrust (mN)", data=sim_results)


# Import experimental mass flow rate data
data1 = pd.read_csv('/home/josh/nozzleDesign/testdata/11132019_test11.csv')
data1.insert(1, "P init (psia)", [x+14.7 for x in data1["P init (psig)"]], True)
data1.insert(3, "P fin (psia)", [x+14.7 for x in data1["P fin (psig)"]], True)
data1.insert(4, "P avg (psia)", [x for x in (data1["P init (psia)"] + data1["P fin (psia)"])/2], True)

# data1.insert(8, "Beta min init", [round(x,4) for x in (data1["Ma init (g)"]-0.3)/(data1["Mb init (g)"]+0.03)])
# data1.insert(9, "Beta max init", [round(x,4) for x in (data1["Ma init (g)"]+0.3)/(data1["Mb init (g)"]-0.03)])

# data1.insert(12, "Beta min final", [round(x,4) for x in (data1["Ma final (g)"]-0.3)/(data1["Mb final (g)"]+0.03)])
# data1.insert(13, "Beta max final", [round(x,4) for x in (data1["Ma final (g)"]+0.3)/(data1["Mb final (g)"]-0.03)])

data1.insert(10, "dM", [round(x,2) for x in data['M final (g)']-data['M init (g)']])
# Adding or subtracting uncertain data means the uncertanties add *in quadrature* (sqaure-root of the sum of the squares)
# This assumes the uncertainy provided by the manufacturer is given in terms of standard deviation
# Thus, the error is: +/- sqrt(0.02) grams for all data points

# data1.insert(10, "dM min", [round(x,2) for x in (data1["Beta min init"]+1)*data1["Mb init (g)"] + 
# 												(data1["Beta min init"]-1)*0.03 -
# 												(data1["Beta max final"]+1)*data1["Mb final (g)"] +
# 												(data1["Beta max final"]-1)*0.03])

# data1.insert(11, "dM max", [round(x,2) for x in (data1["Beta max init"]+1)*data1["Mb init (g)"] - 
# 												(data1["Beta max init"]-1)*0.03 -
# 												(data1["Beta min final"]+1)*data1["Mb final (g)"] -
# 												(data1["Beta min final"]-1)*0.03])

# data1.insert(16, "dM/dt min", [ round(x,3) for x in data1["dM min"]/data1["Time (s)"] ])
# data1.insert(17, "dM/dt max", [ round(x,3) for x in data1["dM max"]/data1["Time (s)"] ])
# data1.insert(18, "dM/dt nom", [ round(x,2) for x in (data1["Ma init (g)"] + data1["Mb init (g)"] - data1["Ma final (g)"] - data1["Mb final (g)"])/data1["Time (s)"] ])
# data1.insert(19, "dM/dt err", [ x for x in data1["dM/dt max"]-data1["dM/dt min"] ])

data1.insert(11, "dM/dt", [round(x,3) for x in data1["dM"]/data1["Time (s)"]])
# Multiplying or dividing uncertain data is a little more complicated
# All measurements must be put in the form of *fractional* uncertainties, then added in quadrature
# Thus the error is: +/- sqrt((dm_err/dm)**2 + (dt_err/dt)**2) grams/sec for each data point




# Setup theoretical vs. experimental plot
    # Blue: #1f77b4
    # Orange: #ff7f0e
    # Green: #2ca02c

fig1, ax1 = plt.subplots(figsize=(8.5, 5), dpi=90)
ax1.set_xlabel('Pressure (psia)', color='#413839')
ax1.set_ylabel('Mass Flow Rate (g/s)', color='#413839')

ax1.plot("Pressure (psig)", "mdot (g/s)", data=sim_results, color='#1f77b4', label='isentropic')
ax1.errorbar(data1["P avg (psia)"], data1["dM/dt nom"], xerr=3, yerr=data1["dM/dt err"]/2, color='#2ca02c', label='experimental', linestyle='none', marker='x')

slope, intercept, r_value, p_value, std_err = stats.linregress(data1["P avg (psia)"], data1["dM/dt nom"])
ax1.plot(data1["P avg (psia)"], slope*data1["P avg (psia)"]+intercept, color='#2ca02c', linestyle='--', label='_nolegend_')

ax1.tick_params(colors='#413839')
ax1.grid(which='major', axis='both', linestyle='--')


ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

# ax2.set_ylabel('Temperature (K)', color='#413839')
# ax2.plot(list_of_P_ts, list_of_T_ts, color='#ff7f0e', label='Critical Temperature (K)')

box = ax1.get_position()
ax1.set_position([box.x0, box.y0 + box.height*0.1, box.width, box.height*0.9])

fig1.legend(['Inviscid', 'Experimental', 'Temperature'], loc='center', bbox_to_anchor=(0.5, 0.03), ncol=3, frameon=False )
plt.title('Mass Flow Rate Comparison ({} mm)'.format(d_star), y=1.03, color='#413839')







# # Thrust vs. Pressure (Inviscid and Experimental)
# # test8_trial1 = thrust_data('testdata/11062019_thrust_test1.csv', 97.56)
# # thrust_data2 = thrust_data('11062019_thrust_test2.csv', 97.61)
# # thrust_data3 = thrust_data('11062019_thrust_test3.csv', 97.51)
# test8_trial7 = thrust_data('testdata/11062019_test8_thrust_trial7.csv', 98.18)
# test8_trial10 = thrust_data('testdata/11062019_test8_thrust_trial10.csv', 97.71)
# test8_trial11 = thrust_data('testdata/11062019_test8_thrust_trial11.csv', 97.62)
# # test8_trial12 = thrust_data('testdata/11062019_test8_thrust_trial12.csv', 97.66)

# # test9_trial1 = thrust_data('testdata/11072019_test9_thrust_trial1.csv', 144.33)
# # test9_trial5 = thrust_data('testdata/11072019_test9_thrust_trial5.csv', 144.65)
# test9_trial9 = thrust_data('testdata/11072019_test9_thrust_trial9.csv', 145.07)

# data_points = [test8_trial7, test8_trial10, test8_trial11, test9_trial9]

# fig8, ax1 = plt.subplots(figsize=(6.5, 4), dpi=90)
# 	# Blue: #1f77b4 (Inviscid)
# 	# Green: #2ca02c
# 	# Orange: #ff7f0e
# 	# Gray: #808080
# ax1.set_xlabel('Pressure (psia)', color='#413839')
# ax1.set_ylabel('Thrust (mN)', color='#413839')
# ax1.plot(list_of_P_ts, [x*1000 for x in list_of_thrusts], color='#1f77b4', label='isentropic')

# for name in data_points:
# 	ax1.plot(name["Pressure (psia)"], name["Thrust (mN)"], color='#2ca02c', label='trial1', linestyle='none', marker='x')

# all_thrust_data = pd.concat(data_points)

# slope, intercept, r_value, p_value, std_err = stats.linregress(all_thrust_data["Pressure (psia)"], all_thrust_data["Thrust (mN)"])
# ax1.plot(all_thrust_data["Pressure (psia)"], slope*all_thrust_data["Pressure (psia)"]+intercept, color='#2ca02c', linestyle='--', label='_nolegend_')

# ax1.tick_params(colors='#413839')
# ax1.grid(which='major', axis='both', linestyle='--')

# box = ax1.get_position()
# ax1.set_position([box.x0, box.y0 + box.height*0.1, box.width, box.height*0.9])

# fig8.legend(['Inviscid', 'Experimental'], loc='center', bbox_to_anchor=(0.5, 0.03), ncol=3, frameon=False )
# # plt.title('Thrust Comparison ({} mm)'.format(d_star), y=1.03, color='#413839')

# ---- Plot Everything
plt.show()