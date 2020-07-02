# Nozzle master control
# Given a nozzle geometry and initial conditions, this code will sweep through a range of stagnation pressures and output the exit conditions
# It's interesting to note that the critical conditions are dependent ONLY on geometry and not stagnation conditions
from nozzle_code_2 import nozzle
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
P_t_init = 114.7 * 6894.76  		# Init Total Pressure, units of Pa (psia * 6894.76)
P_amb = 14.7 * 6894.76  			# Ambient Pressure, units of Pa (psia * 6894.76)
T_t_init = 0 + 273.15  		# Init Total Temperature, units of K (C + 273.15)
# list_of_Tts = [223, 293]
vol = 25 / 10**6  			# Plenum volume, units of m^3 (cm^3 / 10^6)
cutoff_cond = 0.00001
list_of_time_steps = [0.01]
d_star = 0.5 / 1000  			# Nozzle throat diameter, units of m (mm / 1000)
half_angle = 10  				# (Conical) Nozzle expansion angle (degrees)
expansion_ratio = 1.1850			# Nozzle expansion ratio (Exit Area / Throat Area)
								# 	Inlet PSI ------- Ideal Expansion Ratio
								# 		114.7 ------- 1.8048
								# 		80 ---------- 1.6173
								# 		65 ---------- 1.3225
								#	Impulse is maximized when Exp Ratio = ~1.1850
figsize = (6, 6)				# Figure size (in)
dpi = 150						# Figure dpi


## ==================================================================================
## ---- PRE-CALCULATIONS ------------------------------------------------------------
## ==================================================================================

if gas_type == 'R236fa':
	k = 1.083  # Heat capacity ratio (Cp/Cv)
	R = 8.314/0.152039  # Specific gas constant (J/kg-K)

if gas_type == 'R134a':
	gas_label = gas_type
	right_limit = 30
	k = 1.127  # Override value to compare results to MST paper
	R = 8.314/0.10203  # Specific gas constant (J/kg-K)
	T_trip = 169.85  # Triple point temperature (K)
	P_trip = 389.56  # Triple point pressure (Pa)

if gas_type == 'N2':
	k = 1.039
	R = 8.314/0.028014  # Specific gas constant (J/kg-K)
	
if gas_type == 'CO2':
	gas_label = 'CO$_{2}$'
	right_limit = 3.
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


data_dict = pd.DataFrame()

for time_step_init in list_of_time_steps:

	## ==================================================================================
	## ---- SET UP DATA LISTS -----------------------------------------------------------
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
	# 4. Use new density to update P, T assuming polytropic process + isentropic + ideal gas # NO NO NO. 'isentropic' assumes NO EXCHANGE OF MATTER. THIS IS INVALID.
	# 5. Repeat 1-4 until P < 35

	# for P_t in list_of_P_ts:
	i = 0
	n_max = 0
	end_loop_flag = False
	time_step = time_step_init

	while delta_pres > cutoff_cond and list_of_P_ts[-1] > P_amb and time[-1] < right_limit:
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
			delta_pres = np.absolute((list_of_P_ts[-1] - list_of_P_ts[-2])/list_of_P_ts[-2])/time_step

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

		print('Time step: ' + str(round(time_step, 6)) + ' sec, P_t: ' + str(round(list_of_P_ts[-1]/ 6894.76, 1)) + ' psia, Change in pres: ' + str(round(delta_pres, 5)) + ' ' + str(round(time[-1], 4)) + ' sec', end='\r', flush=True)

	# By the nature of this loop, anything that has an init value will end up with one extra element in its list
	# So we must manually remove the last element once all is said and done in order to make all the array lengths the same
	# This is due to the fact that the last element will be calculated for Pt < P_amb, which is not realistic.
	del m_gas[-1], list_of_chamber_densities[-1], list_of_T_ts[-1], list_of_P_ts[-1], time[-1]

	single_iter = pd.DataFrame(data={'time': [x for x in time], 'P_t': [(x / 1000) for x in list_of_P_ts], 'thrust': [x * 1000 for x in list_of_thrusts], 'net impulse': [x * 1000 for x in cumulative_impulse]})
	single_iter['Time Step'] = '{} s'.format(time_step)
	data_dict = data_dict.append(single_iter)
	# label = str(time_step)
	# data_dict[label] = single_iter


# ## ==================================================================================
# ## ---- POST-CALCULATIONS -----------------------------------------------------------
# ## ==================================================================================

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
## ---- PRESSURE & THRUST VS TIME ---------------------------------------------------
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

# test_nos = [ 
# 			'20191204_143449',  # 3x trials @ 5x pressures, 0.6 mm nozzle, raw data in g (multiply by 9.81)
# 			'20191204_143538',
# 			'20191204_143636',

# 			'20191204_143735',
# 			'20191204_143817',
# 			'20191204_143931',

# 			'20191204_144142',
# 			'20191204_144207',
# 			'20191204_144256',

# 			'20191204_144715',
# 			'20191204_144754',
# 			'20191204_144840',

# 			'20191204_145400',
# 			'20191204_145419',
# 			'20191204_145442'
# 			]  
# steady_state = True
# mult_by_g = True

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




# ## ---- Single plenum discharge tests -----------------------------------------------

test_nos = [
			'20191130_131419', # 114.7 psia, 0.6 mm nozzle, raw data in g (multiply by 9.81)
			'20191130_131515',
			'20191130_131607',
			'20191130_131624',
			'20191130_131644'
			]  
steady_state = False
mult_by_g = True


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
fontsize = 10
if steady_state:
	data_marker = 'None'
else:
	data_marker = 'o'

fig1, axs = plt.subplots(3, 1, figsize=figsize, dpi=dpi, sharex='col')
if steady_state:
	fig1.suptitle('Steady-State Pressure & Thrust Measurements\n({0} Trials) ($\\varnothing${1} mm, $\\lambda$={2})'.format(len(test_nos), d_star*1000, expansion_ratio), fontsize=fontsize)
else:
	# fig1.suptitle('Single Plenum Discharge Pressure & Thrust Measurements\n({0} Trials) ($\\varnothing${1} mm, $\\lambda$={2})'.format(len(test_nos), d_star*1000, expansion_ratio), fontsize=fontsize)
	fig1.suptitle('Simulation Results Convergence for Varying Time Step', y=0.95)
	axs[0].set_title(r'({}, $V_{{p}}=${} cm$^3$, Nozzle $\varnothing${} mm, $\lambda$={})'.format(gas_label, vol*10**6, d_star*1000, expansion_ratio), fontsize=9)

td1 = []  # Pressure
# for trial, test_no in reversed(list(enumerate(test_nos))):
# 	test_data = all_data('Test Data/' + test_no, mult_by_g)[0]  # All the Float Pressure data
# 	test_data['Trial'] = trial  # Used to label the data for showing individually
# 	test_data['Setpoint'] = '{} psig'.format(int(test_data['Float Pressure (psig)'][0].round(-1)))  # Used to label the data for showing individually
# 	td1.append(test_data)
# td1 = pd.concat(td1)
# td1['Time (s)'] = td1['Time (s)'].round(1)
# sns.lineplot(ax=axs[0],
# 			 x='Time (s)',
# 			 y='Float Pressure (psig)',
# 			 data=td1,
# 			#  hue='Setpoint', 
# 			#  style='Setpoint',  # Show each trial individually instead of an aggregate
# 			#  estimator=np.mean,  # Show each trial individually instead of an aggregate
# 			 marker=data_marker
# 			 )
			 
if not steady_state:
	sns.lineplot(ax=axs[0],
				 data=data_dict,
				 x='time',
				 y='P_t',
				 hue='Time Step'
	)
	# axs[0].plot([x+0.57 for x in data_dict[str(time_step)]['time']], [(x / 6894.76) - 14.7 for x in data_dict[str(time_step)]['P_t']], color='#ff7f0e', label='pres_inlet', linestyle='-', linewidth=linewidth)
	axs[0].set_xlim(left=0, right=right_limit)
	axs[0].set_ylim(bottom=0)

axs[0].set_ylabel('Pressure, $kPa$', color='#413839', fontsize=fontsize)
axs[0].tick_params(colors='#413839')
axs[0].grid(which='major', axis='both', linestyle='--')
# axs[0].set_ylim(bottom=0)
axs[0].legend(loc=1, fontsize=9)


box0 = axs[0].get_position()
# axs[0].set_position([box0.x0 + box0.width * 0.05, box0.y0 + box0.height * 0.05, box0.width, box0.height])


td2 = []  # Thrust
# for trial, test_no in reversed(list(enumerate(test_nos))):
# 	test_data = all_data('Test Data/' + test_no, mult_by_g)[2]  # All the Thrust data
# 	test_data['Trial'] = trial  # Used to label the data for showing individually
# 	test_data['Setpoint'] = '{} psig'.format(int(all_data('Test Data/' + test_no, mult_by_g)[0]['Float Pressure (psig)'][0].round(-1)))  # Used to label the data for showing individually
# 	td2.append(test_data)
# td2 = pd.concat(td2)
# td2['Time (s)'] = td2['Time (s)'].round(1)
# sns.lineplot(ax=axs[1],
# 			 x='Time (s)',
# 			 y='Thrust Corrected (mN)',
# 			 data=td2,
# 			#  hue='Setpoint',  # Show each trial individually instead of an aggregate
# 			#  style='Setpoint',  # Show each trial individually instead of an aggregate
# 			#  estimator=np.mean,  
# 			 marker=data_marker
# 			 )

if not steady_state:
	sns.lineplot(ax=axs[1],
				x='time',
				y='thrust',
				data=data_dict,
				hue='Time Step'
	)
	# axs[1].plot([x+0.57 for x in data_dict[str(time_step)]['time']], [x *1000 for x in data_dict[str(time_step)]['thrust']], color='#ff7f0e', label='thrust', linestyle='-', linewidth=linewidth)
	axs[1].set_xlim(left=0, right=right_limit)
	axs[1].set_ylim(bottom=0)

axs[1].set_xlabel('Time, s', color='#413839', fontsize=fontsize)
axs[1].set_ylabel('Thrust, $mN$', color='#413839', fontsize=fontsize)
axs[1].tick_params(colors='#413839')
axs[1].grid(which='major', axis='both', linestyle='--')
axs[1].legend(loc=1, fontsize=9)


# Net impulse
if not steady_state:
	sns.lineplot(ax=axs[2],
				x='time',
				y='net impulse',
				data=data_dict,
				hue='Time Step'
	)
	# axs[1].plot([x+0.57 for x in data_dict[str(time_step)]['time']], [x *1000 for x in data_dict[str(time_step)]['thrust']], color='#ff7f0e', label='thrust', linestyle='-', linewidth=linewidth)
	axs[2].set_xlim(left=0, right=right_limit)
	axs[2].set_ylim(bottom=0)

axs[2].set_xlabel('Time, s', color='#413839', fontsize=fontsize)
axs[2].set_ylabel('Impulse, $mN-s$', color='#413839', fontsize=fontsize)
axs[2].tick_params(colors='#413839')
axs[2].grid(which='major', axis='both', linestyle='--')
axs[2].legend(loc=1, fontsize=9)


box1 = axs[2].get_position()
# axs[1].set_position([box1.x0 + box1.width * 0.05, box1.y0 + box1.height * 0.05, box1.width, box1.height])

fig1.align_ylabels()




## ==================================================================================
## ---- THRUST VS PRESSURE ----------------------------------------------------------
## ==================================================================================
# test_nos = [ 
# 			'20191204_143449',  # 3x trials @ 5x pressures, 0.6 mm nozzle, raw data in g (multiply by 9.81)
# 			'20191204_143538',
# 			'20191204_143636',

# 			'20191204_143735',
# 			'20191204_143817',
# 			'20191204_143931',

# 			'20191204_144142',
# 			'20191204_144207',
# 			'20191204_144256',

# 			'20191204_144715',
# 			'20191204_144754',
# 			'20191204_144840',

# 			'20191204_145400',
# 			'20191204_145419',
# 			'20191204_145442'
# 			]  
# steady_state = True
# mult_by_g = True

# linewidth = 2
# fontsize = 12
# if steady_state:
# 	data_marker = 'None'
# else:
# 	data_marker = 'o'

# td1 = []  # Pressure
# for trial, test_no in reversed(list(enumerate(test_nos))):
# 	test_data = all_data('Test Data/' + test_no, mult_by_g)[0]  # All the Float Pressure data
# 	test_data['trial'] = trial  # Used to label the data for showing individually
# 	test_data['Setpoint'] = '{} psig'.format(int(test_data['Float Pressure (psig)'][0].round(-1)))  # Used to label the data for showing individually
# 	td1.append(test_data)
# td1 = pd.concat(td1)
# pressure_data = pd.concat([td1['Float Pressure Resampled (psig)'][12], td1['Float Pressure Resampled (psig)'][22]]) # psig, at 1.2 and 2.2 seconds

# td2 = []  # Thrust
# for trial, test_no in reversed(list(enumerate(test_nos))):
# 	test_data = all_data('Test Data/' + test_no, mult_by_g)[2]  # All the Thrust data
# 	test_data['trial'] = trial  # Used to label the data for showing individually
# 	test_data['Setpoint'] = '{} psig'.format(int(all_data('Test Data/' + test_no, mult_by_g)[0]['Float Pressure (psig)'][0].round(-1)))  # Used to label the data for showing individually
# 	td2.append(test_data)
# td2 = pd.concat(td2)
# thrust_data = pd.concat([td2['Thrust Resampled (mN)'][12], td2['Thrust Resampled (mN)'][22]])  # mN, at 1.2 and 2.2 seconds


# fig8, ax1 = plt.subplots(figsize=figsize, dpi=dpi)
# 	# Blue: #1f77b4 (Inviscid)
# 	# Green: #2ca02c
# 	# Orange: #ff7f0e
# 	# Gray: #808080
# ax1.set_xlabel('Pressure (psig)', color='#413839')
# ax1.set_ylabel('Thrust (mN)', color='#413839')
# colors = ['#1f77b4', '#ff7f0e']
# linestyles = ['-', ':']

# for i, T in enumerate(list_of_Tts):
# 	label = 'Inviscid (' + str(T-273) + ' C)'
# 	pres = [(x/6894.76)-14.7 for x in list_of_P_ts]
# 	thrust = list(data_dict[str(T)]['thrust']*1000)
# 	ax1.plot(pres, thrust, color=colors[i], label=label, linestyle=linestyles[i])

# ax1.plot([x for x in pressure_data], thrust_data, color='#2ca02c', linestyle='none', label='Experimental', marker='x')

# ax1.tick_params(colors='#413839')
# ax1.grid(which='major', axis='both', linestyle='--')
# ax1.set_xlim(left=0)
# ax1.set_ylim(bottom=0)

# box = ax1.get_position()
# ax1.set_position([box.x0, box.y0 + box.height*0.1, box.width, box.height*0.9])
# fig8.legend(loc='center', bbox_to_anchor=(0.5, 0.03), ncol=3, frameon=False )

# # fig8.legend(['Inviscid', 'Experimental'], loc='center', bbox_to_anchor=(0.5, 0.03), ncol=3, frameon=False )
# plt.title('Ideal & Measured Thrust vs. Pressure\n($\\varnothing${0} mm, $\\lambda$={1})'.format(d_star*1000, expansion_ratio), y=1.0, color='#413839')



# ---- Plot Everything
plt.show()
