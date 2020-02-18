# Nozzle Impulse Optimization
# This code will run the standard nozzle_master sequence and sweep through a range of expansion ratios to observe the effect on Net Impulse
from nozzle_code import nozzle
import matplotlib.pyplot as plt
from scipy import stats
import math
import pandas as pd
import numpy as np
import sys
import seaborn as sns
from data_handling_funcs import *
from matplotlib.lines import Line2D


## ==================================================================================
## ---- USER OPTIONS -----------------------------------------------------------------
## ==================================================================================

gas_type = 'CO2'				# Gas Choices: R236fa, R134a, N2, CO2, air
P_t_init = 114.7 * 6894.76  	# Init Total Pressure, units of Pa (psia * 6894.76)
P_amb = 14.7 * 6894.76  		# Ambient Pressure, units of Pa (psia * 6894.76)
T_t_init = 0 + 273.15  			# Init Total Temperature, units of K (C + 273.15)
vol = 26 / 10**6  				# Plenum volume, units of m^3 (cm^3 / 10^6)
time_step = 0.01				# Simulation time step
d_star = 0.18 / 1000  			# Nozzle throat diameter, units of m (mm / 1000)
half_angle = 10  				# (Conical) Nozzle expansion angle (degrees)
bit_tip_dia = 0.1 / 1000		# (Conical) Engraving bit tip diameter, used to determine drill depth for optimized nozzle expansion ratio

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


list_of_expansion_ratios = [x/100 for x in np.arange(100, 200, 2).tolist()]  # Gotta do it like this to circumvent floating point precision errors w/ the .tolist method

list_of_exit_areas = [x*(np.pi*d_star**2)/4 for x in list_of_expansion_ratios]
list_of_exit_diameters = [d_star*np.sqrt(x) for x in list_of_expansion_ratios]
list_of_nozzle_lengths = [((x-d_star)/2)/np.tan(np.radians(half_angle)) for x in list_of_exit_diameters]
list_of_drill_depths = [((x-bit_tip_dia)/2)/np.tan(np.radians(half_angle)) for x in list_of_exit_diameters]

list_of_cumulative_impulses = []

## ==================================================================================
## ---- GO --------------------------------------------------------------------------
## ==================================================================================
for expansion_ratio in list_of_expansion_ratios:
	print('Running at \u03B5 = {}'.format(expansion_ratio))

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
	# 4. Use new density to update P, T assuming polytropic process + isentropic + ideal gas # NO NO NO. "isentropic" assumes NO EXCHANGE OF MATTER. THIS IS INVALID.
	# 5. Repeat 1-4 until P < 35

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

	list_of_cumulative_impulses.append(cumulative_impulse[-1])


## ==================================================================================
## ---- PLOT ------------------------------------------------------------------------
## ==================================================================================

# ---- Plot 2x3 [Thrust, Impulse, ISP, Re, Ma, Density] all vs. Inlet + Throat + Exit Pressure

linewidth = 2
fontsize = 12

fig1, axs = plt.subplots(2, 1, figsize=figsize, dpi=dpi, sharex='col')
fig1.suptitle('Net Impulse & Drill Depth vs. Expansion Ratio')

axs[0].plot(list_of_expansion_ratios, [x*1000 for x in list_of_cumulative_impulses], color='#ff7f0e', linestyle='-', linewidth=linewidth)
# axs[0].set_title('Required Drill Depth for Varying Expansion Ratio')
axs[0].set_ylabel('Net Impulse, mN-s', color='#413839', fontsize=fontsize)
axs[0].tick_params(colors='#413839')
axs[0].grid(which='major', axis='both', linestyle='--')
box0 = axs[0].get_position()
axs[0].set_position([box0.x0 + box0.width * 0.05, box0.y0 + box0.height * 0.05, box0.width, box0.height])

axs[1].plot(list_of_expansion_ratios, [x*1000 for x in list_of_drill_depths], color='#ff7f0e', linestyle='-', linewidth=linewidth)
axs[1].set_xlabel('Expansion Ratio, \u03B5', color='#413839', fontsize=fontsize)
axs[1].set_ylabel('Drill Depth, mm', color='#413839', fontsize=fontsize)
axs[1].tick_params(colors='#413839')
axs[1].grid(which='major', axis='both', linestyle='--')
box1 = axs[1].get_position()
axs[1].set_position([box1.x0 + box1.width * 0.05, box1.y0 + box1.height * 0.05, box1.width, box1.height])

fig1.align_ylabels()
# plt.tight_layout()
plt.show()