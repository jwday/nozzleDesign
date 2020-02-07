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


## ==================================================================================
## ---- USER OPTIONS -----------------------------------------------------------------
## ==================================================================================

gas_type = 'CO2'				# Gas Choices: R236fa, R134a, N2, CO2, air
P_t_init = 114.7 * 6894.76  	# Init Total Pressure, units of Pa (psia * 6894.76)
P_amb = 14.7 * 6894.76  		# Ambient Pressure, units of Pa (psia * 6894.76)
T_t_init = -30 + 273.15  			# Init Total Temperature, units of K (C + 273.15)
vol = 30 / 10**6  				# Plenum volume, units of m^3 (cm^3 / 10^6)
time_step = 0.001				# Simulation time step
d_star = 0.4 / 1000  			# Nozzle throat diameter, units of m (mm / 1000)
half_angle = 10  				# (Conical) Nozzle expansion angle (degrees)
expansion_ratio = 1.3225		# Nozzle expansion ratio (Exit Area / Throat Area)
								# 	Inlet PSI ------- Ideal Expansion Ratio
								# 		114.7 ------- 1.8048
								# 		80 ---------- 1.6173
								# 		65 ---------- 1.3225
								#	Impulse is maximized when Exp Ratio = ~1.1850
figsize = (6, 4.5)				# Figure size (in)
dpi = 150						# Figure dpi


## ==================================================================================
## ---- PRE-CALCULATIONS ------------------------------------------------------------
## ==================================================================================

if gas_type == 'R236fa':
	k = 1.083  # Heat capacity ratio (Cp/Cv)
	R = 8.314/0.152039  # Specific gas constant (J/kg-K)

if gas_type == 'R134a':
	k = 1.127  # Override value to compare results to MST paper
	R = 8.314/0.10203  # Specific gas constant (J/kg-K)

if gas_type == 'N2':
	k = 1.039
	R = 8.314/0.028014  # Specific gas constant (J/kg-K)
	
if gas_type == 'CO2':
	k = 1.289
	R = 8.314/0.04401  # Specific gas constant (J/kg-K)

if gas_type == 'air':
	k = 1.401
	R = 8.314/0.0289645  # Specific gas constant (J/kg-K)

density_init = P_t_init/(R*(T_t_init))  # Units of kg/m^3 (or g/l)
m_init = density_init*vol  # Units of kg
dia = 2*(vol*(3/4)/math.pi)**(1/3)  # Plenum diatmer , units of m


## ==================================================================================
## ---- SET UP DATA LISTS -----------------------------------------------------------
## ==================================================================================

# list_of_P_ts = list(np.linspace (P_t_max, P_amb, no_of_points))
list_of_P_ts = [P_t_init]
list_of_T_ts = [T_t_init]
list_of_chamber_densities = [density_init]
list_of_mdots = []

list_of_P_exits = []
list_of_v_exits = []
list_of_M_exits = []
list_of_thrusts = []
list_of_pressure_ratios = []
# list_of_exit_shock_PRs = []
list_of_P_stars = []
list_of_T_stars = []
list_of_rho_stars = []
list_of_T_exits = []
list_of_rho_exits = []
list_of_Re_stars = []

time = [0]
m_gas = [m_init]



## ==================================================================================
## ---- EXECUTE LOOP ----------------------------------------------------------------
## ==================================================================================
# 0. Start with init P, T, m_gas
# 1. Run nozzle given P, T
# 2. Return m_dot, use to update m_gas assuming a specific time step
# 3. Use m_gas to determine new density
# 4. Use new density to update P, T assuming polytropic process + isentropic + ideal gas # NO NO NO. 'isentropic' assumes NO EXCHANGE OF MATTER. THIS IS INVALID.
# 5. Repeat 1-4 until P < 35

# for P_t in list_of_P_ts:
while list_of_P_ts[-1] > P_amb:
	m_dot, M_crit_sub, M_crit_sup, PR_crit_sub, PR_crit_sup, PR_exit_shock, M_exit_behindshock, M_exit, P_exit, v_exit, F, P_star, T_star, rho_star, Re_star, T_exit, rho_exit = nozzle(list_of_P_ts[-1], list_of_T_ts[-1], P_amb, d_star, expansion_ratio, half_angle, gas_type)

	list_of_mdots.append(m_dot*1000)  # Units of g/s
	list_of_P_exits.append(P_exit)
	list_of_v_exits.append(v_exit)
	list_of_M_exits.append(M_exit)
	list_of_thrusts.append(F)
	list_of_pressure_ratios.append(P_exit/list_of_P_ts[-1])
	list_of_P_stars.append(P_star)
	list_of_T_stars.append(T_star)
	list_of_rho_stars.append(rho_star)
	list_of_Re_stars.append(Re_star)
	list_of_T_exits.append(T_exit)
	list_of_rho_exits.append(rho_exit)

	time.append(time[-1] + time_step)
	m_gas.append(m_gas[-1] - m_dot*time_step)
	list_of_chamber_densities.append(m_gas[-1]/vol)

	# THIS IS INVALD. You have to...
	# 1. Determine the NEW enthalpy after some of the mass has left.
	# list_of_T_ts.append( list_of_T_ts[-1]*(list_of_chamber_densities[-1]/list_of_chamber_densities[-2])**(k-1) )
	list_of_T_ts.append( list_of_T_ts[-1] )  # Ideal gas model states that, for FREE EXPANSION OF A GAS, temperature is constant
	# list_of_P_ts.append( list_of_P_ts[-1]*(list_of_chamber_densities[-1]/list_of_chamber_densities[-2])**k )
	list_of_P_ts.append( list_of_chamber_densities[-1]*R*list_of_T_ts[-1] )


 # By the nature of this loop, anything that has an init value will end up with one extra element in its list
 # So we must manually remove the last element once all is said and done in order to make all the array lengths the same
 # This is due to the fact that the last element will be calculated for Pt < P_amb, which is not realistic.
del list_of_P_ts[-1], list_of_T_ts[-1], list_of_chamber_densities[-1], time[-1], m_gas[-1]



## ==================================================================================
## ---- POST-CALCULATIONS -----------------------------------------------------------
## ==================================================================================

cumulative_impulse = []
cumulative_mass = []
average_thrust = []
time_offset = []
ISP = []

for i in range(0, len(time)-1):
	time_offset.append( np.average([time[i], time[i+1]]) )
	if i == 0:
		cumulative_impulse.append( time_step*np.average([list_of_thrusts[i], list_of_thrusts[i+1]]) )
		cumulative_mass.append( time_step*np.average([list_of_mdots[i], list_of_mdots[i+1]]) )
	else:
		cumulative_impulse.append( time_step*np.average([list_of_thrusts[i], list_of_thrusts[i+1]]) + cumulative_impulse[i-1] )
		cumulative_mass.append( time_step*np.average([list_of_mdots[i], list_of_mdots[i+1]]) + cumulative_mass[i-1] )

	average_thrust.append( np.average(list_of_thrusts[0:i+1]) )

for i in range(0, len(time)):
	ISP.append( 1000*list_of_thrusts[i]/(9.81*list_of_mdots[i]) )




## ==================================================================================
## ---- PLOT ------------------------------------------------------------------------
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

test_nos = [ 
			'20191204_143449',  # 3x trials @ 5x pressures, 0.6 mm nozzle, raw data in g (multiply by 9.81)
			'20191204_143538',
			'20191204_143636',

			'20191204_143735',
			'20191204_143817',
			'20191204_143931',

			'20191204_144142',
			'20191204_144207',
			'20191204_144256',

			'20191204_144715',
			'20191204_144754',
			'20191204_144840',

			'20191204_145400',
			'20191204_145419',
			'20191204_145442'
			]  
steady_state = True
mult_by_g = True

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




## ---- Single plenum discharge tests -----------------------------------------------

# test_nos = [
# 			'20191130_131419', # 114.7 psia, 0.6 mm nozzle, raw data in g (multiply by 9.81)
# 			'20191130_131515',
# 			'20191130_131607',
# 			'20191130_131624',
# 			'20191130_131644'
# 			]  
# steady_state = False
# mult_by_g = True


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

linewidth = 2
fontsize = 12
if steady_state:
	data_marker = 'None'
else:
	data_marker = 'o'

fig1, axs = plt.subplots(2, 1, figsize=figsize, dpi=dpi, sharex='col')
fig1.suptitle('Pressure & Thrust vs. Time ({} Trials)'.format(len(test_nos)), fontsize=fontsize)

td1 = []
for trial, test_no in enumerate(test_nos):
	test_data = all_data('Test Data/' + test_no, mult_by_g)[0]  # All the Float Pressure data
	test_data['trial'] = trial  # Used to label the data for showing individually
	td1.append(test_data)
td1 = pd.concat(td1)
td1['Time (s)'] = td1['Time (s)'].round(1)
sns.lineplot(ax=axs[0],
			 x='Time (s)',
			 y='Float Pressure (psia)',
			 data=td1,
			 hue='trial',  estimator=None, # Show each trial individually instead of an aggregate
			 marker=data_marker)
axs[0].set_ylabel('Pressure, psia', color='#413839', fontsize=fontsize)
axs[0].tick_params(colors='#413839')
axs[0].grid(which='major', axis='both', linestyle='--')

if not steady_state:
	axs[0].plot([x+0.57 for x in time], [x / 6894.76 for x in list_of_P_ts], color='#ff7f0e', label='pres_inlet', linestyle='-', linewidth=linewidth)
	axs[0].set_xlim(left=0, right=5)

box0 = axs[0].get_position()
# axs[0].set_position([box0.x0 + box0.width * 0.05, box0.y0 + box0.height * 0.05, box0.width, box0.height])


td2 = []
for trial, test_no in enumerate(test_nos):
	test_data = all_data('Test Data/' + test_no, mult_by_g)[2]  # All the Thrust data
	test_data['trial'] = trial  # Used to label the data for showing individually
	td2.append(test_data)
td2 = pd.concat(td2)
td2['Time (s)'] = td2['Time (s)'].round(1)
sns.lineplot(ax=axs[1],
			 x='Time (s)',
			 y='Thrust Corrected (mN)',
			 data=td2,
			 hue='trial', estimator=None,  # Show each trial individually instead of an aggregate
			 marker=data_marker)
axs[1].set_xlabel('Time, s', color='#413839', fontsize=fontsize)
axs[1].set_ylabel('Thrust, mN', color='#413839', fontsize=fontsize)
axs[1].tick_params(colors='#413839')
axs[1].grid(which='major', axis='both', linestyle='--')

if not steady_state:
	axs[1].plot([x+0.57 for x in time], [x * 1000 for x in list_of_thrusts], color='#ff7f0e', label='thrust', linestyle='-', linewidth=linewidth)
	axs[1].set_xlim(left=0, right=12)

box1 = axs[1].get_position()
# axs[1].set_position([box1.x0 + box1.width * 0.05, box1.y0 + box1.height * 0.05, box1.width, box1.height])

fig1.align_ylabels()

# ---- Plot Everything
plt.show()
