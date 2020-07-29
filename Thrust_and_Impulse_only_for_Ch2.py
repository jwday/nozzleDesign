# Nozzle master control
# This code is only for plotting nozzle and impulse for Chapter 2 of thesis. Nothing else.
from nozzle_code_2 import nozzle
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
vol = 60 / 10**6  				# Plenum volume, units of m^3 (cm^3 / 10^6)
d_star = 0.6 / 1000  			# Nozzle throat diameter, units of m (mm / 1000)
cutoff_cond = 0.0001			# Cutoff condition, defined by the fractional change in pressure (relative to P_t_init) per second, units of 1/sec
half_angle = 10  				# (Conical) Nozzle expansion angle (degrees)
expansion_ratio = 1.17			# Nozzle expansion ratio (Exit Area / Throat Area)
								# 	Inlet PSI ------- Ideal Expansion Ratio
								# 		114.7 ------- 1.8048
								# 		80 ---------- 1.6173
								# 		65 ---------- 1.3225
								#	Impulse is maximized when Exp Ratio = ~1.1850 (1.17?)
figsize = (7.5, 4)				# Figure size (in)
dpi = 150						# Figure dpi




## ==================================================================================
## ---- SETUP PROPELLANT ------------------------------------------------------------
## ==================================================================================

if gas_type == 'R236fa':
	k = 1.083  # Heat capacity ratio (Cp/Cv)
	R = 8.314/0.152039  # Specific gas constant (J/kg-K)

if gas_type == 'R134a':
	gas_label = gas_type
	right_limit = 28
	k = 1.127  # Override value to compare results to MST paper
	R = 8.314/0.10203  # Specific gas constant (J/kg-K)
	T_trip = 169.85  # Triple point temperature (K)
	P_trip = 389.56  # Triple point pressure (Pa)

if gas_type == 'N2':
	k = 1.039
	R = 8.314/0.028014  # Specific gas constant (J/kg-K)
	
if gas_type == 'CO2':
	gas_label = 'CO$_{2}$'
	right_limit = 1.4*2
	k = 1.289
	R = 8.314/0.04401  # Specific gas constant (J/kg-K)
	T_trip = 216.58  # Triple point temperature (K)
	P_trip = 518500  # Triple point pressure (Pa)

if gas_type == 'H2':
	k = 1.410
	R = 8.314/0.002016  # Specific gas constant (J/kg-K)

if gas_type == 'air':
	k = 1.401
	R = 8.314/0.0289645  # Specific gas constant (J/kg-K)

density_init = P_t_init/(R*(T_t_init))  # Units of kg/m^3 (or g/l)
m_init = density_init*vol  # Units of kg
dia = 2*(vol*(3/4)/math.pi)**(1/3)  # Plenum diatmer , units of m

time_step_init = 15.28*vol*(P_t_init-P_amb)/((P_t_init*d_star)**2)


## ==================================================================================
## ---- INIT DATA LISTS -------------------------------------------------------------
## ==================================================================================

# list_of_P_ts = list(np.linspace (P_t_max, P_amb, no_of_points))
list_of_P_ts = [P_t_init]  # These need to be precribed at the first time step so the nozzle function can calculate the rest of the parameters at the same time step
list_of_T_ts = [T_t_init]
list_of_chamber_densities = [density_init]
m_gas = [m_init]

time = [0]

list_of_mdots = []
list_of_P_exits = []
list_of_v_exits = []
list_of_M_exits = []
list_of_thrusts = []
list_of_dthrust = []
list_of_pressure_ratios = []
list_of_P_stars = []
list_of_T_stars = []
list_of_rho_stars = []
list_of_Re_stars = []
list_of_T_exits = []
list_of_rho_exits = []


average_thrust = []
cumulative_impulse = []
ISP = []
list_of_dthrust = []
cumulative_mass = []
delta_pres = 1

## ==================================================================================
## ---- EXECUTE LOOP ----------------------------------------------------------------
## ==================================================================================
# 0. Start with init P, T, m_gas
# 1. Run nozzle given P, T
# 2. Return m_dot, use to update m_gas assuming a specific time step
# 3. Use m_gas to determine new density
# 4. Use new density to update P, T assuming polytropic process + isentropic + ideal gas
# 5. Repeat 1-4 until P < 35
i = 0
n_max = 0
end_loop_flag = False
time_step = time_step_init

while delta_pres > cutoff_cond and list_of_P_ts[-1] > P_amb:
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

	average_thrust.append( np.average(list_of_thrusts) )
	ISP.append( 1000*list_of_thrusts[-1]/(9.81*list_of_mdots[i]) )

	if i == 0:  # If we're on the first time step...
		list_of_dthrust.append(list_of_thrusts[-1]/time_step)  # Thrust starts here, so the abs. value will be positive and large. Doesn't matter what this actually is, though, because the end loop conditional won't ever see it
		cumulative_impulse.append( time_step*list_of_thrusts[-1])
	else:  # If we're on any other time step...
		list_of_dthrust.append((list_of_thrusts[-1] - list_of_thrusts[-2])/time_step)
		cumulative_impulse.append(time_step*list_of_thrusts[-1] + cumulative_impulse[-1])
	
	# print('P_t: ' + str(round(list_of_P_ts[-1]/6894.76, 1)) + ' psia at ' + str(round(time[-1], 3)) + ' sec', end='\r', flush=True)

	# Calculate these properties in preparation for the next loop. Loop will end if the newly calculated pressure drops below the ambient pressure.
	m_gas.append(m_gas[-1] - m_dot*time_step) 
	list_of_chamber_densities.append(m_gas[-1]/vol)
	list_of_T_ts.append( list_of_T_ts[-1]*(list_of_chamber_densities[-1]/list_of_chamber_densities[-2])**(k-1) )
	list_of_P_ts.append( list_of_P_ts[-1]*(list_of_chamber_densities[-1]/list_of_chamber_densities[-2])**k )
	time.append(i*time_step)  # The first iteration is at t=0, so the first time[] entry will be 0.

	# print('Avg thrust: ' + str(round(average_thrust[-1]*1000, 3)) + ' mN at ' + str(round(time[-1], 3)) + ' sec', end='\r', flush=True)

	##### END LOOP CONDITIONALS TO CHECK FOR EXIT SHOCK#####
	# if i>=2:  # Don't even bother checking end loop conditionals until after the first two iterations
	# 	# Except for the first time step, dthrust should always be NEGATIVE. Also, dthrust should continue to get LESS NEGATIVE as time goes on (because chamber pressure is decreasing therefore thrust is decreasing too)
	# 	if list_of_dthrust[-1] < list_of_dthrust[-2]:  # This should only occur ONCE, when a shock appears at the exit. It's at this point we expect the thrust to suddenly drop a LOT
	# 		end_loop_flag = True  # Set this to true
	# 		n_max = int(i + 5/time_step)
	# 		print('Shock found!\n')


	if i>=2:
		delta_pres = np.absolute((list_of_P_ts[-2] - list_of_P_ts[-1])/P_t_init)/time_step

	# 	if 1000*cumulative_impulse[-1] > 110 and 1000*cumulative_impulse[-2] < 110:
	# 		n_target = i
	# 		end_loop_flag = True
	# 		n_max = int(i + 5/time_step)
	# 		print('Target impulse reached!\n')

	# if not end_loop_flag:  # If the end loop flag has not been set, then continue on as normal and increment the max steps
	# 	n_max+=1
	# else:
	# 	pass

	i+=1

	print('Time step: ' + str(round(time_step, 6)) + ' sec, P_t: ' + str(round(list_of_P_ts[-1]/ 6894.76, 1)) + ' psia, Change in pres: ' + str(round(delta_pres*100, 3)) + '%/sec ' + str(round(time[-1], 4)) + ' sec', end='\r', flush=True)

print('\n')

# By the nature of this loop, anything that has an init value will end up with one extra element in its list
# So we must manually remove the last element once all is said and done in order to make all the array lengths the same
del m_gas[-1], list_of_chamber_densities[-1], list_of_T_ts[-1], list_of_P_ts[-1], time[-1]



## ==================================================================================
## ---- POST-CALCULATIONS -----------------------------------------------------------
## ==================================================================================



# for i in range(0, len(time)-1):
# 	time_offset.append( np.average([time[i], time[i+1]]) )
# 	if i == 0:
# 		cumulative_impulse.append( time_step*np.average([list_of_thrusts[i], list_of_thrusts[i+1]]) )
# 		# cumulative_mass.append( time_step*np.average([list_of_mdots[i], list_of_mdots[i+1]]) )
# 	else:
# 		# cumulative_mass.append( time_step*np.average([list_of_mdots[i], list_of_mdots[i+1]]) + cumulative_mass[i-1] )

# 	average_thrust.append( np.average(list_of_thrusts[0:i+1]) )

# 	print('Avg thrust: ' + str(round(average_thrust[-1]*1000, 3)) + ' mN at ' + str(round(time[i], 3)) + ' sec', end='\r', flush=True)

# print('\n')
# for i in range(0, len(time)):

# 	print('ISP: ' + str(round(ISP[-1], 1)) + ' mN at ' + str(round(time[i], 3)) + ' sec', end='\r', flush=True)



## ==================================================================================
## ---- PLOT ------------------------------------------------------------------------
## ==================================================================================

# ---- Plot 2x3 [Thrust, Impulse, ISP, Re, Ma, Density] all vs. Inlet + Throat + Exit Pressure

linewidth = 2
fontsize = 12

data = 		{ 
			#   'pressure': [x/1000 for x in list_of_P_ts],
			  'thrust': [x*1000 for x in list_of_thrusts], 
			#   'impulse': [x*1000 for x in cumulative_impulse], 	
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
			  'impulse': time, 							
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


fig, axs = plt.subplots(2, 1, figsize=figsize, dpi=dpi, sharex=True)

axs[0].plot(time, [x/1000 for x in list_of_P_ts], color='#1f77b4', label='Reservoir Pressure', linestyle='-', linewidth=linewidth)
axs[0].set_ylabel('Pressure, $kPa$', color='#413839', fontsize=fontsize)
axs[0].tick_params(colors='#413839')
axs[0].grid(which='major', axis='both', linestyle='--')

# axs[1].plot(time, [x for x in ISP], color='#1f77b4', label='Specific Impulse', linestyle='-', linewidth=linewidth)
# axs[1].set_xlabel('Time, $s$', color='#413839', fontsize=fontsize)
# axs[1].set_ylabel('Specific Impulse, $s$', color='#413839', fontsize=fontsize)
# axs[1].tick_params(colors='#413839')
# axs[1].grid(which='major', axis='both', linestyle='--')

axs[1].plot(time, [x*1000 for x in list_of_thrusts], color='#1f77b4', label='Instantaneous Thrust', linestyle='-', linewidth=linewidth)
axs[1].plot(time, [x*1000 for x in average_thrust], label='Cumulative Average Thrust', linestyle='--')
axs[1].set_ylabel('Thrust, $mN$', color='#413839', fontsize=fontsize)
axs[1].tick_params(colors='#413839')
axs[1].grid(which='major', axis='both', linestyle='--')
handles, labels = axs[1].get_legend_handles_labels()
axs[1].legend(handles, labels, loc='best', fontsize=10)

# axs[2].plot(time, [x*1000 for x in cumulative_impulse], color='#1f77b4', label='Net Impulse', linestyle='-', linewidth=linewidth)
# axs[2].set_xlabel('Time, $s$', color='#413839', fontsize=fontsize)
# axs[2].set_ylabel('Net Impulse, $mN-s$', color='#413839', fontsize=fontsize)
# axs[2].tick_params(colors='#413839')
# axs[2].grid(which='major', axis='both', linestyle='--')

# box1 = axs[1].get_position()
# axs.set_position([box1.x0, box1.y0 + box1.height * 0.2, box1.width * 0.95, box1.height * 0.78])

# Plot indicator lines
# target_thrust_line = [11, 11]
# target_thrust_time = [-1, 28]
# target_avg_thrust = axs.plot(target_thrust_time, target_thrust_line, label='Target Avg Thrust', color='black', linestyle=':', linewidth=linewidth/2)


axs[0].set_xlim(left=0, right=right_limit)
axs[0].set_ylim(bottom=0)
axs[1].set_ylim(bottom=0)
# axs[2].set_ylim(bottom=0)

# handles2, labels2 = axs[1].get_legend_handles_labels()
# handles.insert(2, handles2[0])
# labels.insert(2, labels2[0])

fig.suptitle('In-Space Specific Impulse Change', y=0.98)
axs[0].set_title(r'({} at $T_0$={} K, $V_{{p}}=${} cm$^3$, Nozzle $\varnothing${} mm, $\lambda$={})'.format(gas_label, T_t_init, vol*10**6, d_star*1000, expansion_ratio), fontsize=9)
fig.canvas.set_window_title('Nozzle Performance Metrics')

print('')
print('t_step: ' + str(time_step) + ' sec')
# print('Exit shock P_t: ' + str(round(list_of_P_ts[int(n_max - 5/time_step)], 3)))
print('Net impulse: ' + str(round(cumulative_impulse[-1]*1000, 3)) + ' mN-s')
print('Net avg thrust: ' + str(round(average_thrust[-1]*1000, 3)) + ' mN')
print('Total time: ' + str(round(time[-1], 3)) + ' s')
if end_loop_flag:
	print('Time at 110 mN-s: ' + str(time[n_target]) + ' sec')
	print('Thrust at 110 mN-s: ' + str(1000*average_thrust[n_target]) + ' mN')
print('')
plt.show()
