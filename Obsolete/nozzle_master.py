# Nozzle master control
# This code will estimate the performance of a nozzle + plenum assembly for given starting conditions and nozzle geometry
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
P_t_init = 114.7 * 6894.76  		# Init Total Pressure, units of Pa (psia * 6894.76)
P_amb = 14.7 * 6894.76  			# Ambient Pressure, units of Pa (psia * 6894.76)
T_t_init = 0 + 273.15  		# Init Total Temperature, units of K (C + 273.15)
vol = 30 / 10**6  			# Plenum volume, units of m^3 (cm^3 / 10^6)
d_star = 0.6 / 1000  			# Nozzle throat diameter, units of m (mm / 1000)
cutoff_cond = 0.0001			# Cutoff condition, defined by the fractional change in pressure (relative to P_t_init) per second, units of 1/sec
half_angle = 10  				# (Conical) Nozzle expansion angle (degrees)
expansion_ratio = 1.17			# Nozzle expansion ratio (Exit Area / Throat Area)
								# 	Inlet PSI ------- Ideal Expansion Ratio
								# 		114.7 ------- 1.8048
								# 		80 ---------- 1.6173
								# 		65 ---------- 1.3225
								#	Impulse is maximized (90.6 mN-s) when Exp Ratio = ~1.17 for CO2 @ 114.7psi in, 14.7 out, 0C, 25cc vol, 0.4mm dia, 10 deg half-angle
								# 		1.17	90.6	100.0%
								#		1.30	90.3	 99.7%
								#		1.40	89.5	 98.8%
								#		1.50	88.6	 97.8%
								#		1.60	87.5	 96.6%
								#		1.70	86.3	 95.3%
								#		1.80	85.0	 93.8%
								#		1.90	83.7	 
								#		2.0		82.3
figsize = (6, 5)				# Figure size (in)
dpi = 150						# Figure dpi




## ==================================================================================
## ---- SETUP PROPELLANT ------------------------------------------------------------
## ==================================================================================

if gas_type == 'R236fa':
	k = 1.083  # Heat capacity ratio (Cp/Cv)
	R = 8.314/0.152039  # Specific gas constant (J/kg-K)

if gas_type == 'R134a':
	gas_label = gas_type
	right_limit = 40
	k = 1.127  # Override value to compare results to MST paper
	R = 8.314/0.10203  # Specific gas constant (J/kg-K)
	T_trip = 169.85  # Triple point temperature (K)
	P_trip = 389.56  # Triple point pressure (Pa)

if gas_type == 'N2':
	k = 1.039
	R = 8.314/0.028014  # Specific gas constant (J/kg-K)
	
if gas_type == 'CO2':
	gas_label = 'CO$_{2}$'
	right_limit = 1.3
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
list_of_P_ts = [P_t_init]  # These need to be prescribed at the first time step so the nozzle function can calculate the rest of the parameters at the same time step
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
# 4. Use new density to update P, T assuming polytropic process + isentropic + ideal gas # NO NO NO. "isentropic" assumes NO EXCHANGE OF MATTER. THIS IS INVALID.
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
	ISP.append( 1000*list_of_thrusts[-1]/(9.81*list_of_mdots[-1]) )

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
# This is due to the fact that the last element will be calculated for Pt < P_amb, which is not realistic.
del m_gas[-1], list_of_chamber_densities[-1], list_of_T_ts[-1], list_of_P_ts[-1], time[-1]



## ==================================================================================
## ---- POST-CALCULATIONS -----------------------------------------------------------
## ==================================================================================

# cumulative_impulse = []
# cumulative_mass = []
# average_thrust = []
# time_offset = []
# ISP = []

# for i in range(0, len(time)-1):
# 	time_offset.append( np.average([time[i], time[i+1]]) )
# 	if i == 0:
# 		cumulative_impulse.append( time_step*np.average([list_of_thrusts[i], list_of_thrusts[i+1]]) )
# 		cumulative_mass.append( time_step*np.average([list_of_mdots[i], list_of_mdots[i+1]]) )
# 	else:
# 		cumulative_impulse.append( time_step*np.average([list_of_thrusts[i], list_of_thrusts[i+1]]) + cumulative_impulse[i-1] )
# 		cumulative_mass.append( time_step*np.average([list_of_mdots[i], list_of_mdots[i+1]]) + cumulative_mass[i-1] )

# 	average_thrust.append( np.average(list_of_thrusts[0:i+1]) )

# for i in range(0, len(time)):
# 	ISP.append( 1000*list_of_thrusts[i]/(9.81*list_of_mdots[i]) )




## ==================================================================================
## ---- PLOT ------------------------------------------------------------------------
## ==================================================================================

# ---- Plot 2x3 [Thrust, Impulse, ISP, Re, Ma, Density] all vs. Inlet + Throat + Exit Pressure

linewidth = 2
fontsize = 12

data = 		{ 
			  'pressure': [x/1000 for x in list_of_P_ts],
			  'thrust': [x*1000 for x in list_of_thrusts], 
			  'impulse': [x*1000 for x in cumulative_impulse], 	
			#   'isp': ISP, 			
			#   'mach_exit': list_of_M_exits, 		
			#   'P_exit': [x/1000 for x in list_of_P_exits],	
			#   'reynolds': [x/1000 for x in list_of_Re_stars],
			#   't_exit': list_of_T_exits,
			#   'm_dot': list_of_mdots
			  }

figname = 	{
			  'pressure': 'Pressure',
			  'thrust': 'Thrust', 							
			  'impulse': 'Net Impulse', 							
			  'isp': '$I_{SP}$', 			
			  'mach_exit': 'Exit Mach Number', 	
			  'P_exit': 'Exit Pressure', 		
			  'reynolds': 'Throat Reynold\'s Number',
			  't_exit': 'Exit Temperature',
			  'm_dot': 'Mass Flow Rate'
			}

times = 	{ 
			  'pressure': time,
			  'thrust': time, 								
			  'impulse': time, 							
			  'isp': time, 			
			  'mach_exit': time, 					
			  'P_exit': time, 					
			  'reynolds': time,
			  't_exit': time,
			  'm_dot': time
			}

ylabels = 	{ 
			  'pressure': 'Pressure, $kPa$',
			  'thrust': 'Thrust, $mN$', 						
			  'impulse': 'Impulse, $mN-s$', 						
			  'isp': '$I_{SP}$, $s$', 		
			  'mach_exit': 'Mach', 				
			  'P_exit': 'Pressure, $kPa$', 	
			  'reynolds': 'Re x $10^3$',
			  't_exit': 'Temperature, $K$',
			  'm_dot': 'Mass Flow Rate, $g/s$'
			}

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

num_rows = 3
num_cols = 1
# figsize = ((8.5/2)*num_cols, (11/4)*num_rows+0.5)				# Figure size (in)

fig, axs = plt.subplots(num_rows, num_cols, figsize=figsize, dpi=dpi)
for i, j in enumerate(data.items()):
	if num_rows > 1 and num_cols > 1:
		row = math.floor(i/2 % (len(data)/2))
		col = i % num_cols
		row_col = (row, col)
	else:
		row_col = (i)

	# --- For Pressure on Secondary Axis ----
	axs[row_col].plot(times[j[0]], data[j[0]], color='#1f77b4', label=ylabels[j[0]], linestyle='-', linewidth=linewidth)
	axs[row_col].set_title(figname[j[0]])
	axs[row_col].set_xlabel('Time, $s$', color='#413839', fontsize=fontsize)
	axs[row_col].set_ylabel(ylabels[j[0]], color='#413839', fontsize=fontsize)
	axs[row_col].tick_params(colors='#413839')
	axs[row_col].set_ylim(bottom=0)
	
	if j[0] == 'P_exit':
		axs[row_col].plot([times[j[0]][0], times[j[0]][-1]], [P_trip/1000, P_trip/1000], color='#FF0000', label='Critical Pressure', linestyle='-', linewidth=linewidth)
	if j[0] == 't_exit':
		axs[row_col].plot([times[j[0]][0], times[j[0]][-1]], [T_trip, T_trip], color='#FF0000', label='Critical Temperature', linestyle='-', linewidth=linewidth)
	if j[0] == 'thrust':
		axs[row_col].plot(times[j[0]], [x*1000 for x in average_thrust], color='#FF0000', label='Avg Thrust', linestyle='--', linewidth=linewidth)

	# second_ax = axs[row_col].twinx()
	# line1, = second_ax.plot(time, [x/1000 for x in list_of_P_ts], color='#1f77b4', label='pres_inlet', linestyle='-', linewidth=linewidth)
	# line2, = second_ax.plot(time, [x/1000 for x in list_of_P_stars], color='#1f77b4', label='pres_throat', linestyle=':', linewidth=linewidth)
	# line3, = second_ax.plot(time, [x/1000 for x in list_of_P_exits], color='#1f77b4', label='pres_exit', linestyle='--', linewidth=linewidth)
	# second_ax.set_ylabel('Pressure, kPa', color='#413839', fontsize=fontsize)
	# second_ax.tick_params(colors='#413839')


	# --- For Pressure on Primary Axis ----
	# axs[row_col].plot(time, [x/1000 for x in list_of_P_ts], color='#1f77b4', label='Inlet Pressure, kPa', linestyle='-', linewidth=linewidth)
	# axs[row_col].plot(time, [x/1000 for x in list_of_P_stars], color='#1f77b4', label='Throat Pressure, kPa', linestyle=':', linewidth=linewidth)
	# axs[row_col].plot(time, [x/1000 for x in list_of_P_exits], color='#1f77b4', label='Exit Pressure, kPa', linestyle='--', linewidth=linewidth)
	# axs[row_col].set_title(figname[j[0]])
	# axs[row_col].set_xlabel('Time, $s$', color='#413839', fontsize=fontsize)
	# axs[row_col].set_ylabel('Pressure, kPa', color='#413839', fontsize=fontsize)
	# axs[row_col].tick_params(colors='#413839')

	# second_ax = axs[row_col].twinx()
	# second_ax.plot(times[j[0]], data[j[0]], color='#ff7f0e', label=ylabels[j[0]], linestyle='-.', linewidth=linewidth)
	# second_ax.set_ylabel(ylabels[j[0]], color='#413839', fontsize=fontsize)
	# second_ax.tick_params(colors='#413839')
	# second_ax.set_ylim(bottom=0)

	# box = axs[row_col].get_position()
	# axs[row_col].set_position([box.x0, box.y0, box.width * 1.05, box.height * 1.1])
	# axs[row_col].grid(which='major', axis='both', linestyle='--')

# if num_rows > 1 and num_cols > 1:
# 	axs[0, 0].plot(time, [x*1000 for x in average_thrust], linestyle='--')
# else:
# 	axs[0].plot(time, [x*1000 for x in average_thrust], linestyle='--')

fig.canvas.set_window_title('Nozzle Performance Metrics')
fig.suptitle('Performance Metrics for Single Plenum Discharge\n'r'({} cm$^3$, $\varnothing${} mm, $\lambda$={})'.format(vol*10**6, d_star*1000, expansion_ratio))
# fig.legend((line1, line2, line3), ('Inlet Pressure', 'Throat Pressure', 'Exit Pressure'), loc='center', bbox_to_anchor=(0.5, 0.05), ncol=3, frameon=True, fontsize=fontsize, edgecolor='255', borderpad=1)
plt.tight_layout()
plt.subplots_adjust(bottom=0.15)
plt.subplots_adjust(top=0.8)
plt.show()

# print('Nozzle Exit Pressure: ', list_of_P_exits[0] / 6894.76)


# ---- 1) Plot Inlet Pressures and Exit Velocity

# fig1, ax1 = plt.subplots(figsize=(8, 5), dpi=90)

# ax1.set_xlabel('Time (s)', color='#413839')
# ax1.set_ylabel('Pressure (psia)', color='#413839')
# ax1.plot(time, list_of_P_ts, color='#1f77b4', label='Inlet Pressure (psia)')
# ax1.plot(time, list_of_P_exits, color='#ff7f0e', label='Exit Pressure (psia)')
# ax1.tick_params(colors='#413839')
# # ax1.set_ylim(0, 90)

# ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

# ax2.set_ylabel('Velocity (m/s)', color='#413839')  # we already handled the x-label with ax1
# ax2.plot(time, list_of_v_exits, color='#2ca02c', label='Exit Velocity (m/s)')
# # ax2.set_ylim(0, 350)

# box = ax1.get_position()
# ax1.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
# fig1.legend(['Inlet Pressure (psia)', 'Exit Pressure (psia)', 'Exit Velocity (m/s)'], loc='center', bbox_to_anchor=(0.5, 0.03), ncol=3, frameon=False )

# ax1.grid(which='major', axis='both', linestyle='--')

# plt.title('Chamber + Exit Pres. & Exit Velocity for Single Plenum Discharge', y=1.03, color='#413839')



# ---- 2) Plot Instantaneous Thrust, Average Thrust, Total Impulse, and ISP

# fig2, ax1 = plt.subplots(figsize=(8, 5), dpi=90)
# ax1.set_xlabel('Time (s)', color='#413839')
# ax1.set_ylabel('Thrust (mN)', color='#413839')
# ax1.plot(time, [x*1000 for x in list_of_thrusts], color='#1f77b4', label='Instantaneous Thrust (mN)')
# ax1.plot(time_offset, [x*1000 for x in average_thrust], color='#1f77b4', label='Cumulative Avg Thrust (mN)', linestyle='--')
# ax1.tick_params(colors='#413839')
# # ax1.set_ylim(0, 0.5)

# ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

# ax2.set_ylabel('Total Impulse (mN-s)', color='#413839')  # we already handled the x-label with ax1
# ax2.plot(time_offset, [x*1000 for x in cumulative_impulse], color='#2ca02c', label='Total Impulse (mN-s)')
# # ax2.set_ylim(0, 0.1)

# ax3 = ax1.twinx()
# ax3.set_ylabel('ISP', color='#413839')
# ax3.plot(time, ISP, color='#ff7f0e', label='ISP')
# ax3.spines['right'].set_position(('outward', 60))      

# box = ax1.get_position()
# ax1.set_position([box.x0, box.y0 + box.height * 0.1, box.width*0.88, box.height * 0.9])
# fig2.legend(['Instantaneous Thrust (mN)', 'Cumulative Avg Thrust (mN)', 'Total Impulse (mN-s)', 'ISP'], loc='center', bbox_to_anchor=(0.5, 0.03), ncol=4, frameon=False )

# ax1.grid(which='major', axis='both', linestyle='--')

# plt.title('Thrust & Total Impulse for Single Plenum Discharge', y=1.03, color='#413839')



# ---- Plot Mass Flow Rate and Total Mass

# fig3, ax1 = plt.subplots(figsize=(8, 5), dpi=90)
# ax1.set_xlabel('Time (s)', color='#413839')
# ax1.set_ylabel('Mass Flow Rate (g/s)', color='#413839')
# ax1.plot(time, list_of_mdots, color='#1f77b4', label='Mass Flow Rate (g/s)')
# ax1.tick_params(colors='#413839')
# # ax1.set_ylim(0, 1.3)

# ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

# ax2.set_ylabel('Propellant Consumed (g)', color='#413839')  # we already handled the x-label with ax1
# ax2.plot(time_offset, cumulative_mass, color='#2ca02c', label='Propellant Consumed (g)')
# # ax2.set_ylim(0, 0.3)

# box = ax1.get_position()
# ax1.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
# fig3.legend(['Mass Flow Rate (g/s)', 'Propellant Consumed (g)'], loc='center', bbox_to_anchor=(0.5, 0.03), ncol=2, frameon=False )

# ax1.grid(which='major', axis='both', linestyle='--')

# plt.title('Mass Flow Rate & Total Propellant Mass Consumed for Single Plenum Discharge', y=1.03, color='#413839')



# ---- Plot Critical Pressure, Critical Temperature, and Critical Density

# fig4, ax1 = plt.subplots(figsize=(8.5, 5), dpi=90)
# ax1.set_xlabel('Time (s)', color='#413839')
# ax1.set_ylabel('Temperature (K)', color='#413839')
# ax1.tick_params(colors='#413839')
# ax1.plot(time, list_of_T_stars, color='#1f77b4', label='Critical Temperature (K)')
# # ax1.plot(time, list_of_T_exits, color='#1f77b4', label='Exit Temperature (K)')

# ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
# ax2.set_ylabel('Density (g/mL)', color='#413839')  # we already handled the x-label with ax1
# ax2.plot(time, [x/1000 for x in list_of_rho_stars], color='#2ca02c', label='Critical Density (g/mL)')
# # ax2.plot(time, list_of_rho_exits, color='#2ca02c', label='Exit Density (g/mL)')

# ax3 = ax1.twinx()
# ax3.set_ylabel('Pressure (psi)', color='#413839')
# ax3.plot(time, list_of_P_stars, color='#ed9926', label='Critical Pressure (psi)')
# ax3.spines['right'].set_position(('outward', 60))      

# box = ax1.get_position()
# ax1.set_position([box.x0, box.y0 + box.height * 0.1, box.width*0.9, box.height * 0.9])
# fig4.legend(['Critical Temperature (K)', 'Critical Density (g/mL)', 'Critical Pressure (psi)'], loc='center', bbox_to_anchor=(0.5, 0.03), ncol=3, frameon=False )

# ax1.grid(which='major', axis='both', linestyle='--')

# plt.title('Temperature and Density for Single Plenum Discharge', y=1.03, color='#413839')





# figs = {'fig1': fig1, 'fig2': fig2, 'fig3': fig3, 'fig4': fig4}

# ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
# ax2.set_ylabel('Thrust, mN', color='#413839', fontsize=fontsize)
# ax2.plot(time, [x*1000 for x in list_of_thrusts], color='#ff7f0e', label='Instantaneous Thrust (mN)', linestyle='--', linewidth=linewidth)
# # ax2.plot(time_offset, [x*1000 for x in average_thrust], color='#ff7f0e', label='Cumulative Avg Thrust (mN)', linestyle='-.')
# ax2.tick_params(colors='#413839')
# ax2.set_ylim(0, 300)

# ax3 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
# ax3.set_ylabel('Total Impulse, mN-s', color='#413839', fontsize=fontsize)  # we already handled the x-label with ax1
# ax3.plot(time_offset, [x*1000 for x in cumulative_impulse], color='#2ca02c', label='Total Impulse (N-s)', linestyle=(0, (3, 1, 1, 1, 1, 1)), linewidth=linewidth)
# ax3.spines['right'].set_position(('outward', 60))

# ax4 = ax1.twinx()
# ax4.set_ylabel('ISP', color='#413839', fontsize=fontsize)
# ax4.plot(time, ISP, color='#cc0000', label='ISP', linestyle=(0, (3, 1, 1, 1)), linewidth=linewidth)

# box = ax1.get_position()
# ax1.set_position([box.x0, box.y0 + box.height * 0.15, box.width * 0.95, box.height * 0.95])
# fig5.legend(['Inlet Pressure', 'Exit Pressure', 'Instantaneous Thrust', 'Cumulative Avg Thrust', 'Total Impulse'], loc='center', bbox_to_anchor=(0.5, 0.08), ncol=3, frameon=False )
# fig5.legend(['Inlet Pressure', 'Exit Pressure', 'ISP'], loc='center', bbox_to_anchor=(0.5, 0.03), ncol=3, frameon=False, fontsize=fontsize)

# ax1.grid(which='major', axis='both', linestyle='--')
# plt.title('Single Plenum Discharge\n{0} mm Nozzle, {1} cm^3 Plenum'.format(d_star, vol*(100**3)), y=1.0, color='#413839')
# fig5.suptitle('Single Plenum Discharge, Pressure and Impulse', fontsize=14, y=0.985)
# fig5.suptitle('\t{0} mm Nozzle, {1} $cm^3$ Plenum'.format(d_star, vol*(100**3)), fontsize=10)
# plt.savefig('/mnt/d/OneDrive - UC Davis/HRVIP/Writing/AIAA SciTech 2019 Paper/Images/Sim Results/image.png')



# Inlet Plot Thrust vs. Pressure
# fig6, ax1 = plt.subplots(figsize=(8.5, 5), dpi=90)
# ax1.set_xlabel('Pressure (psia)', color='#413839')
# ax1.set_ylabel('Thrust (N)', color='#413839')
# ax1.plot(list_of_P_ts, list_of_thrusts, color='#1f77b4', label='Thrust (N)')
# ax1.tick_params(colors='#413839')

# ax1.grid(which='major', axis='both', linestyle='--')
# plt.title('Thrust vs. Inlet Pressure ({} mm)'.format(d_star), y=1.03, color='#413839')




# # Mass Flow Rate vs. Pressure (Isentropic vs. Experimental)

# massflow_data3 = massflow_data_singlescale('Test Data/10222019_test3_singlescale.csv')
# massflow_data4 = massflow_data('Test Data/10302019_test4.csv')
# massflow_data5 = massflow_data('Test Data/10312019_test5.csv')
# massflow_data6 = massflow_data('Test Data/10312019_test6.csv')
# massflow_data7 = massflow_data('Test Data/11012019_test7.csv')

# massflow_data = [massflow_data3, massflow_data4, massflow_data5, massflow_data6, massflow_data7]

# fig7, ax1 = plt.subplots(figsize=(8.5, 5), dpi=90)
# 	# Blue: #1f77b4
# 	# Orange: #ff7f0e
# 	# Green: #2ca02c
# ax1.set_xlabel('Pressure (psia)', color='#413839')
# ax1.set_ylabel('Mass Flow Rate (g/s)', color='#413839')
# ax1.plot(list_of_P_ts, list_of_mdots, color='#1f77b4', label='isentropi	c')

# for name in massflow_data:
# 	ax1.errorbar(name["Pressure (psia)"], name["dM/dt nom"], xerr=2, yerr=name["dM/dt err"]/2, color='#2ca02c', label='experimental', linestyle='none', marker='x')

# all_massflow_data = pd.concat(massflow_data, sort=False)

# slope, intercept, r_value, p_value, std_err = stats.linregress(all_massflow_data["Pressure (psia)"], all_massflow_data["dM/dt nom"])
# ax1.plot(all_massflow_data["Pressure (psia)"], slope*all_massflow_data["Pressure (psia)"]+intercept, color='#2ca02c', linestyle='--', label='_nolegend_')

# ax1.tick_params(colors='#413839')
# ax1.grid(which='major', axis='both', linestyle='--')


# # ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

# # ax2.set_ylabel('Temperature (K)', color='#413839')
# # ax2.plot(list_of_P_ts, list_of_T_ts, color='#ff7f0e', label='Critical Temperature (K)')

# box = ax1.get_position()
# ax1.set_position([box.x0, box.y0 + box.height*0.1, box.width, box.height*0.9])

# fig7.legend(['Inviscid', 'Experimental'], loc='center', bbox_to_anchor=(0.5, 0.03), ncol=3, frameon=False )
# plt.title('Mass Flow Rate Comparison ({} mm)'.format(d_star), y=1.03, color='#413839')

# ax1.set_ylim([0, 0.8])
# # ax2.set_ylim([170, 280])



# Mass Flow Rate vs. Pressure (Isentropic vs. Experimental)

# massflow_data11 = [massflow_data_test11('Test Data/11132019_test11.csv', 0.12, 104.7)]

# fig7, ax1 = plt.subplots(figsize=(8.5, 5), dpi=90)
# 	# Blue: #1f77b4
# 	# Orange: #ff7f0e
# 	# Green: #2ca02c
# ax1.set_xlabel('Pressure (psia)', color='#413839')
# ax1.set_ylabel('Mass Flow Rate (g/s)', color='#413839')
# ax1.plot(list_of_P_ts, list_of_mdots, color='#1f77b4', label='isentropic')

# for name in massflow_data11:
# 	ax1.errorbar(name["Pressure (psia)"], name["dM/dt nom"], xerr=2, yerr=name["dM/dt err"]/2, color='#2ca02c', label='experimental', linestyle='none', marker='x')

# all_massflow_data = pd.concat(massflow_data11, sort=False)

# slope, intercept, r_value, p_value, std_err = stats.linregress(all_massflow_data["Pressure (psia)"], all_massflow_data["dM/dt nom"])
# ax1.plot(all_massflow_data["Pressure (psia)"], slope*all_massflow_data["Pressure (psia)"]+intercept, color='#2ca02c', linestyle='--', label='_nolegend_')

# ax1.tick_params(colors='#413839')
# ax1.grid(which='major', axis='both', linestyle='--')


# # ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

# # ax2.set_ylabel('Temperature (K)', color='#413839')
# # ax2.plot(list_of_P_ts, list_of_T_ts, color='#ff7f0e', label='Critical Temperature (K)')

# box = ax1.get_position()
# ax1.set_position([box.x0, box.y0 + box.height*0.1, box.width, box.height*0.9])

# fig7.legend(['Inviscid', 'Experimental'], loc='center', bbox_to_anchor=(0.5, 0.03), ncol=3, frameon=False )
# plt.title('Mass Flow Rate Comparison ({} mm)'.format(d_star), y=1.03, color='#413839')

# ax1.set_ylim([0, 0.8])
# ax2.set_ylim([170, 280])



# Thrust vs. Pressure (Isentropic vs. Experimental)
# test8_trial1 = thrust_data('Test Data/11062019_thrust_test1.csv', 97.56)
# thrust_data2 = thrust_data('11062019_thrust_test2.csv', 97.61)
# thrust_data3 = thrust_data('11062019_thrust_test3.csv', 97.51)
# test8_trial7 = thrust_data('Test Data/11062019_test8_thrust_trial7.csv', 98.18)
# test8_trial10 = thrust_data('Test Data/11062019_test8_thrust_trial10.csv', 97.71)
# test8_trial11 = thrust_data('Test Data/11062019_test8_thrust_trial11.csv', 97.62)
# test8_trial12 = thrust_data('Test Data/11062019_test8_thrust_trial12.csv', 97.66)

# test9_trial1 = thrust_data('Test Data/11072019_test9_thrust_trial1.csv', 144.33)
# test9_trial5 = thrust_data('Test Data/11072019_test9_thrust_trial5.csv', 144.65)
# test9_trial9 = thrust_data('Test Data/11072019_test9_thrust_trial9.csv', 145.07)

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
# plt.title('Thrust Comparison ({} mm)'.format(d_star), y=1.03, color='#413839')




# Single Plenum Discharge
# plenum_trial13 = thrust_data('Test Data/11072019_test9_thrust_trial13.csv', 145.19)
# plenum_trial14 = thrust_data('Test Data/11072019_test9_thrust_trial14.csv', 145.07)

# data_points = [plenum_trial13, plenum_trial14]
# colors = ['#2ca02c', '#ff7f0e']
# labels = ['Trial 1', 'Trial 2']

# fig9, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(6.5, 6), dpi=90)
# ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

# data_prop_psi = pd.read_csv('Test Data/20191130_131515_prop_data.csv', header=1)
# data_prop_psi.insert(2, "Prop Pressure (psia)", [x+14.7 for x in data_prop_psi["Prop Pressure (psig)"]])

# data_float_psi = pd.read_csv('Test Data/20191130_131515_float_data.csv', header=1)
# data_float_psi.insert(2, "Float Pressure (psia)", [x+14.7 for x in data_float_psi["Float Pressure (psig)"]])

# data_weight = pd.read_csv('Test Data/20191130_131515_loadcell_data.csv', header=1)
# data_weight.insert(2, "Thrust (mN)", [x*9.81 for x in data_weight["Weight (?)"]])

# ax1.plot(data_prop_psi['Time (s)'], data_prop_psi['Prop Pressure (psia)'], 
# 	color='#1f77b4', 
# 	label='Upstream Pres.', 
# 	marker='o', 
# 	markersize='5', 
# 	fillstyle='none', 
# 	linestyle='none', 
# 	linewidth='1')
# ax1.plot(data_float_psi['Time (s)'], data_float_psi['Float Pressure (psia)'], 
# 	color='#1f77b4', 
# 	label='Downstream Pres.', 
# 	marker='o', 
# 	markersize='5', 
# 	fillstyle='none', 
# 	linestyle='none', 
# 	linewidth='1')


# ax2.plot(data_weight['Time (s)'], data_weight['Thrust (mN)'], 
# 	color='#ff7f0e', 
# 	label='Measured Thrust', 
# 	marker='o',
# 	markersize='5',
# 	fillstyle='none', 
# 	linestyle='none', 
# 	linewidth='1')


# ax1.plot([x+0.663 for x in time], list_of_P_ts,
# 	color='#1f77b4',
# 	linestyle='-', 
# 	label='Inviscid Pressure')

# ax2.plot([x+0.663 for x in time], [i*1000 for i in list_of_thrusts], 
# 	color='#ff7f0e', 
# 	linestyle='-', 
# 	label='Inviscid Thrust')



# for i, data in enumerate(data_points):
# 	ax1.errorbar(data["Timestamp (ms)"], data["Thrust (mN)"], xerr=0.0625, yerr=18.2466/2, color=colors[i], label=labels[i], linestyle='none', marker='x')


# all_plenum_data = pd.concat(data_points)

# slope, intercept, r_value, p_value, std_err = stats.linregress(all_data["Pressure (psia)"], all_data["Thrust (mN)"])
# ax1.plot(all_data["Pressure (psia)"], slope*all_data["Pressure (psia)"]+intercept, color='#2ca02c', linestyle='--', label='_nolegend_')

# Blue: #1f77b4 (Inviscid)
# Green: #2ca02c
# Orange: #ff7f0e
# Gray: #808080
# ax2.set_xlabel('Time (s)', color='#413839')
# ax1.set_ylabel('Pressure (psia)', color='#413839')
# ax2.set_ylabel('Thrust (mN)', color='#413839')

# ax1.set_xlim([0.5, 2.5])


# ax1.tick_params(colors='#413839')
# ax2.tick_params(colors='#413839')
# ax1.grid(which='major', axis='both', linestyle='--')
# ax2.grid(which='major', axis='both', linestyle='--')

# box1 = ax1.get_position()
# box2 = ax2.get_position()
# ax1.set_position([box1.x0, box1.y0 + box1.height*0.1, box1.width, box1.height*0.9])
# ax2.set_position([box2.x0, box2.y0 + box2.height*0.2, box2.width, box2.height*0.9])

# legend_elements = [Line2D([0], [0], marker='o', color='k', label='Experimental', markersize=5, fillstyle='none', linestyle='none'),
# 				   Line2D([0], [0], linestyle='-', color='k', label='Inviscid')]

# fig9.legend(handles=legend_elements, loc='center', bbox_to_anchor=(0.5, 0.06), ncol=3, frameon=False)
# fig9.legend(['Experimental', 'Inviscid'], loc='center', bbox_to_anchor=(0.5, 0.08), ncol=2, frameon=False)

# ax1.set_title('Experimental vs. Inviscid Pressure & Thrust \n({0} $mm$ Nozzle, {1} $cm^3$ Plenum)'.format(d_star, vol*10**6), y=1.06, color='#413839')
#ax1.set_title('Pressure')
#ax2.set_title('Thrust')



# ---- Plot Everything
# plt.show()
