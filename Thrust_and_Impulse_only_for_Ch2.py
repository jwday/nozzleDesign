# Nozzle master control
# This code is only for plotting nozzle and impulse for Chapter 2 of thesis. Nothing else.
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
## ---- USER OPTIONS ----------------------------------------------------------------
## ==================================================================================

gas_type = 'CO2'				# Gas Choices: R236fa, R134a, N2, CO2, H2, air
P_t_init = 114.7 * 6894.76  	# Init Total Pressure, units of Pa (psia * 6894.76)
P_amb = 14.7 * 6894.76  		# Ambient Pressure, units of Pa (psia * 6894.76)
T_t_init = 0 + 273.15  			# Init Total Temperature, units of K (C + 273.15)
vol = 30 / 10**6  				# Plenum volume, units of m^3 (cm^3 / 10^6)
time_step = 0.01				# Simulation time step
d_star = 0.6 / 1000  			# Nozzle throat diameter, units of m (mm / 1000)
half_angle = 10  				# (Conical) Nozzle expansion angle (degrees)
expansion_ratio = 1.1850		# Nozzle expansion ratio (Exit Area / Throat Area)
								# 	Inlet PSI ------- Ideal Expansion Ratio
								# 		114.7 ------- 1.8048
								# 		80 ---------- 1.6173
								# 		65 ---------- 1.3225
								#	Impulse is maximized when Exp Ratio = ~1.1850
figsize = (6, 6)				# Figure size (in)
dpi = 300						# Figure dpi




## ==================================================================================
## ---- SETUP PROPELLANT ------------------------------------------------------------
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

if gas_type == 'H2':
	k = 1.410
	R = 8.314/0.002016  # Specific gas constant (J/kg-K)

if gas_type == 'air':
	k = 1.401
	R = 8.314/0.0289645  # Specific gas constant (J/kg-K)

density_init = P_t_init/(R*(T_t_init))  # Units of kg/m^3 (or g/l)
m_init = density_init*vol  # Units of kg
dia = 2*(vol*(3/4)/math.pi)**(1/3)  # Plenum diatmer , units of m




## ==================================================================================
## ---- INIT DATA LISTS -------------------------------------------------------------
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




## ==================================================================================
## ---- PLOT ------------------------------------------------------------------------
## ==================================================================================

# ---- Plot 2x3 [Thrust, Impulse, ISP, Re, Ma, Density] all vs. Inlet + Throat + Exit Pressure

linewidth = 2
fontsize = 12

data = 		{ 
			#   'pressure': [x/1000 for x in list_of_P_ts],
			  'thrust': [x*1000 for x in list_of_thrusts], 
			  'impulse': [x*1000 for x in cumulative_impulse], 	
			#   'isp': ISP, 			
			#   'mach_exit': list_of_M_exits, 		
			#   'rho_star': list_of_rho_stars, 		
			#   'reynolds': [x/1000 for x in list_of_Re_stars],
			#   't_star': list_of_T_stars
			  }

figname = 	{ 'pressure': 'Pressure',
			  'thrust': 'Thrust', 							
			  'impulse': 'Net Impulse', 							
			  'isp': '$I_{SP}$', 			
			  'mach_exit': 'Exit Mach Number', 	
			  'rho_star': 'Throat Density', 		
			  'reynolds': 'Throat Reynold\'s Number',
			  't_star': 'Throat Temperature' }

times = 	{ 'pressure': time,
			  'thrust': time, 								
			  'impulse': time_offset, 							
			  'isp': time, 			
			  'mach_exit': time, 					
			  'rho_star': time, 					
			  'reynolds': time,
			  't_star': time }

ylabels = 	{ 'pressure': 'Pressure, $kPa$',
			  'thrust': 'Thrust, $mN$', 						
			  'impulse': 'Impulse, $mN-s$', 						
			  'isp': '$I_{SP}$, $s$', 		
			  'mach_exit': 'Mach', 				
			  'rho_star': 'Density, $kg/m^3$', 	
			  'reynolds': 'Re x $10^3$',
			  't_star': 'Temperature, $K$' }

# legend = 	{ 'thrust': 'Thrust', 							
# 			  'impulse': 'Net Impulse', 							
# 			  'isp': '$I_{SP}$', 			
# 			  'mach_exit': 'Exit Mach Number', 	
# 			  'rho_star': 'Throat Density', 		
# 			  'reynolds': 'Throat Reynold\'s Number' }

# colors = 	{ 'thrust': '#ff7f0e', 							
# 			  'impulse': '#ff7f0e', 								
# 			  'isp': '#ff7f0e', 		
# 			  'mach_exit': '#ff7f0e', 			
# 			  'rho_star': '#ff7f0e', 				
# 			  'reynolds': '#ff7f0e' }


fig, axs = plt.subplots(3, 1, figsize=figsize, dpi=dpi, sharex=True)

axs[0].plot(time, [x/1000 for x in list_of_P_ts], color='#1f77b4', label='Reservoir Pressure', linestyle='-', linewidth=linewidth)
axs[0].set_ylabel('Pressure, $kPa$', color='#413839', fontsize=fontsize)
axs[0].tick_params(colors='#413839')
axs[0].grid(which='major', axis='both', linestyle='--')

axs[1].plot(time, [x*1000 for x in list_of_thrusts], color='#1f77b4', label='Instantaneous Thrust', linestyle='-', linewidth=linewidth)
axs[1].plot(time_offset, [x*1000 for x in average_thrust], label='Cumulative Average Thrust', linestyle='--')
axs[1].set_ylabel('Thrust, $mN$', color='#413839', fontsize=fontsize)
axs[1].tick_params(colors='#413839')
axs[1].grid(which='major', axis='both', linestyle='--')
handles, labels = axs[1].get_legend_handles_labels()
axs[1].legend(handles, labels, loc='upper right', fontsize=10)

axs[2].plot(time_offset, [x*1000 for x in cumulative_impulse], color='#1f77b4', label='Net Impulse', linestyle='-', linewidth=linewidth)
axs[2].set_xlabel('Time, $s$', color='#413839', fontsize=fontsize)
axs[2].set_ylabel('Net Impulse, $mN-s$', color='#413839', fontsize=fontsize)
axs[2].tick_params(colors='#413839')
axs[2].grid(which='major', axis='both', linestyle='--')

# box1 = axs[1].get_position()
# axs.set_position([box1.x0, box1.y0 + box1.height * 0.2, box1.width * 0.95, box1.height * 0.78])

# Plot indicator lines
# target_thrust_line = [11, 11]
# target_thrust_time = [-1, 28]
# target_avg_thrust = axs.plot(target_thrust_time, target_thrust_line, label='Target Avg Thrust', color='black', linestyle=':', linewidth=linewidth/2)


axs[0].set_xlim(left=0, right=1.5)
axs[0].set_ylim(bottom=0)
axs[1].set_ylim(bottom=0)
axs[2].set_ylim(bottom=0)

# handles2, labels2 = axs[1].get_legend_handles_labels()
# handles.insert(2, handles2[0])
# labels.insert(2, labels2[0])

fig.suptitle('Target Thrust, Impulse for Single Plenum Discharge\n'r'({} cm$^3$, $\varnothing${} mm, $\lambda$={})'.format(vol*10**6, d_star*1000, expansion_ratio))
fig.canvas.set_window_title('Nozzle Performance Metrics')

plt.show()