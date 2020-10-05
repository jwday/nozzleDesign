# plot_optimum_impulse.py
# This code will run the nozzle_code sequence and sweep through a range of expansion ratios to observe the effect on net impulse

from nozzle_code_3 import nozzle
import matplotlib.pyplot as plt
from scipy import stats
import math
import pandas as pd
import numpy as np
import sys
import seaborn as sns
from data_handling_funcs import *
from matplotlib.lines import Line2D
import matplotlib as mpl
from nozzle_helperFuncs import *


## ==================================================================================
## ---- USER OPTIONS -----------------------------------------------------------------
## ==================================================================================

gas_type = 'CO2'				# Gas Choices: R236fa, R134a, N2, CO2, H2, air


P_t_init = 82.9 * 6894.76  	# Init Total Pressure, units of Pa (psia * 6894.76)
P_amb = 0 * 6894.76  		# Ambient Pressure, units of Pa (psia * 6894.76)
T_t_init = 20 + 273.15  			# Init Total Temperature, units of K (C + 273.15)
vol = 10 / 10**6  				# Plenum volume, units of m^3 (cm^3 / 10^6)
cutoff_cond = 0.0001				# Cutoff condition, defined by the fractional change in pressure (relative to P_t_init) per second, units of 1/sec
d_star = 0.2 / 1000  			# Nozzle throat diameter, units of m (mm / 1000)
list_of_expansion_ratios = [x/100 for x in np.arange(100, 1000, 10).tolist()]  # Gotta do it like this to circumvent floating point precision errors w/ the .tolist method
list_of_expansion_ratios = [x for x in np.arange(1.0, 10.0, 0.1)]	# Use for CO2
# list_of_expansion_ratios = [x for x in np.logspace(0, 1.6, 20)]	# Use for R134a
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
half_angle = 10  				# (Conical) Nozzle expansion half-angle (degrees)
bit_tip_dia = 0.1 / 1000		# (Conical) Engraving bit tip diameter, used to determine drill depth for optimized nozzle expansion ratio

# figsize = (6,3.5)			# Figure size (in)
dpi = 300						# Figure dpi





## ==================================================================================
## ---- SETUP -----------------------------------------------------------------------
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
	fluid_props = pd.read_excel('R134a_props_NIST.xlsx', sheet_name=None)
	fluid_props_vol = pd.read_excel('R134a_props_const_vol_NIST.xlsx', sheet_name=None)

if gas_type == 'N2':
	k = 1.039
	R = 8.314/0.028014  # Specific gas constant (J/kg-K)
	
if gas_type == 'CO2':
	gas_label = 'CO$_{2}$'
	right_limit = 1.4
	k = 1.289
	R = 8.314/0.04401  # Specific gas constant (J/kg-K)
	T_trip = 216.58  # Triple point temperature (K)
	P_trip = 518500  # Triple point pressure (Pa)
	fluid_props = pd.read_excel('CO2_props_NIST.xlsx', sheet_name=None)
	fluid_props_vol = pd.read_excel('CO2_props_costVol_NIST.xlsx', sheet_name=None)

if gas_type == 'H2':
	k = 1.410
	R = 8.314/0.002016  # Specific gas constant (J/kg-K)

if gas_type == 'air':
	k = 1.401
	R = 8.314/0.0289645  # Specific gas constant (J/kg-K)



h_from_PT_gas_func = create_h_from_PT_gas_func(fluid_props)		# Returns kJ/kg for (Pa, K) input
u_from_PT_gas_func = create_u_from_PT_gas_func(fluid_props)		# Returns kJ/kg for (Pa, K) input
r_from_PT_gas_func = create_r_from_PT_gas_func(fluid_props)		# Returns g/ml for (Pa, K) input

P_from_ru_func, T_from_ru_func = create_PT_from_ru_gas_func(fluid_props_vol)
P_from_rh_func, T_from_rh_func = create_PT_from_rh_gas_func(fluid_props_vol)
P_from_rT_func = create_P_from_rT_gas_func(fluid_props_vol)

# --------------------------------------------------------------------------------
# Calculate initial properties
rho_t_init = r_from_PT_gas_func(P_t_init, T_t_init)[0]*1000				# Initial density, kg/m^3 (or g/l), using REAL DATA
mu_t_init = 1/rho_t_init  												# Initial specific volume, m^3/kg (or l/g)
m_init = rho_t_init*vol  												# Initial propellant mass, kg
u_sp_init = u_from_PT_gas_func(P_t_init, T_t_init)[0]*1000				# Initial specific internal energy, J/kg
h_sp_init = h_from_PT_gas_func(P_t_init, T_t_init)[0]*1000				# Initial specific enthalpy, J/kg
dia = 2*(vol*(3/4)/np.pi)**(1/3)  										# Plenum diameter, m
time_step_init = 15.28*vol*(P_t_init-P_amb)/((P_t_init*d_star)**2)		# Time step, s

# Recalculate P_init based on rho_t_init and h_sp_init. I know it shouldn't be different, but it is, based on whichever function is being used (h_from_PT vs P_from_rh)
P_t_init = P_from_ru_func(rho_t_init, u_sp_init/1000)[0]*1000000
T_t_init = T_from_ru_func(rho_t_init, u_sp_init/1000)[0]



# --------------------------------------------------------------------------------
list_of_exit_areas = [x*(np.pi*d_star**2)/4 for x in list_of_expansion_ratios]
list_of_exit_diameters = [d_star*np.sqrt(x) for x in list_of_expansion_ratios]
list_of_nozzle_lengths = [((x-d_star)/2)/np.tan(np.radians(half_angle)) for x in list_of_exit_diameters]
list_of_nondim_nozzle_lengths = [x/(np.pi*(d_star**2)/4) for x in list_of_nozzle_lengths]
list_of_drill_depths = [((x-bit_tip_dia)/2)/np.tan(np.radians(half_angle)) for x in list_of_exit_diameters]

list_of_cumulative_impulses = []

## ==================================================================================
## ---- GO --------------------------------------------------------------------------
## ==================================================================================
for expansion_ratio in list_of_expansion_ratios:
	print('Running at \u03B5 = {}'.format(expansion_ratio), flush=True)

	## ==================================================================================
	## ---- INIT DATA LISTS -------------------------------------------------------------
	## ==================================================================================

	# list_of_P_ts = list(np.linspace (P_t_max, P_amb, no_of_points))
	list_of_P_ts = [P_t_init]  # These need to be precribed at the first time step so the nozzle function can calculate the rest of the parameters at the same time step
	list_of_T_ts = [T_t_init]
	list_of_rho_ts = [rho_t_init]
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
	time_step = time_step_init  # Just in case you want to try variable time stepping in the future. It didn't work in the past, but also that might be because you were just messing around and didn't know what you were doing.

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
		list_of_rho_ts.append(m_gas[-1]/vol)
		list_of_T_ts.append( list_of_T_ts[-1]*(list_of_rho_ts[-1]/list_of_rho_ts[-2])**(k-1) )
		list_of_P_ts.append( list_of_P_ts[-1]*(list_of_rho_ts[-1]/list_of_rho_ts[-2])**k )
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
	del m_gas[-1], list_of_rho_ts[-1], list_of_T_ts[-1], list_of_P_ts[-1], time[-1]



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

	list_of_cumulative_impulses.append(cumulative_impulse[-1])


## ==================================================================================
## ---- PLOT ------------------------------------------------------------------------
## ==================================================================================
class ScalarFormatterForceFormat(mpl.ticker.ScalarFormatter):
		def _set_format(self):  # Override function that finds format to use.
			self.format = "%1.1f"  # Give format here
sns.axes_style("white")
sns.set_style("whitegrid", {"xtick.major.size": 0, "ytick.major.size": 0, 'grid.linestyle': '--'})
sns.set_context("paper", font_scale = 1, rc={"grid.linewidth": .5})

# ---- Plot 2x3 [Thrust, Impulse, ISP, Re, Ma, Density] all vs. Inlet + Throat + Exit Pressure

linewidth = 2
fontsize = 8

# fig, axs = plt.subplots(2, 1, figsize=figsize, dpi=dpi, sharex='col')
# # fig.suptitle( '          Net Impulse & Drill Depth vs. Expansion Ratio')
# fig.suptitle( 'Net Impulse & Nozzle Length vs. Expansion Ratio', y=0.98)
# # axs[0].set_title(r'({}, $V_{{p}}=${} cm$^3$, Nozzle $\varnothing${} mm)'.format(gas_label, vol*10**6, d_star*1000), fontsize=9)
# axs[0].set_title(r'($P_0$={} kPa, $P_{{amb}}$={} kPa, $V_{{p}}=${} cm$^3$, Nozzle $\varnothing${} mm, Prop: {})'.format(round(P_t_init/1000, 1), round(P_amb/1000, 1), vol*10**6, d_star*1000, gas_label), fontsize=7)
# # axs[0].set_title(r'$P_0 = $' + str(int(P_t_init/1000)) + r' kPa, $P_{amb} = $' + str(int(P_amb/1000)) + ' kPa, Propellent: ' + gas_type, fontsize=10)

# axs[0].plot(list_of_expansion_ratios, [x*1000 for x in list_of_cumulative_impulses], color='#ff7f0e', linestyle='-', linewidth=linewidth)
# # axs[0].set_title('Required Drill Depth for Varying Expansion Ratio')
# axs[0].set_ylabel('Net Impulse, mN-s', color='#413839', fontsize=fontsize)
# axs[0].tick_params(colors='#413839')
# axs[0].grid(which='major', axis='both', linestyle='--')
# box0 = axs[0].get_position()
# axs[0].set_position([box0.x0 + box0.width * 0.05, box0.y0 + box0.height * 0.05, box0.width, box0.height])

# axs[1].plot(list_of_expansion_ratios, [x*1000 for x in list_of_nozzle_lengths], color='#ff7f0e', linestyle='-', linewidth=linewidth)
# axs[1].set_xlabel('Expansion Ratio, \u03B5', color='#413839', fontsize=fontsize)
# # axs[1].set_ylabel('Drill Depth, mm', color='#413839', fontsize=fontsize)
# axs[1].set_ylabel('Nozzle Length, mm', color='#413839', fontsize=fontsize)
# axs[1].tick_params(colors='#413839')
# axs[1].grid(which='major', axis='both', linestyle='--')
# box1 = axs[1].get_position()
# axs[1].set_position([box1.x0 + box1.width * 0.05, box1.y0 + box1.height * 0.05, box1.width, box1.height])

yfmt = ScalarFormatterForceFormat()
# # yfmt.set_powerlimits((0,0))
# axs[0].yaxis.set_major_formatter(yfmt)
# axs[1].yaxis.set_major_formatter(yfmt)
# # axs[i].ticklabel_format(axis='y', style='sci', scilimits=(0,0), useMathText=True)
# axs[0].tick_params(axis='y', labelsize=6, pad=0)
# axs[1].tick_params(axis='y', labelsize=6, pad=0)
# axs[0].yaxis.offsetText.set_fontsize(6)
# axs[1].yaxis.offsetText.set_fontsize(6)

# fig.align_ylabels()
# plt.tight_layout()
# plt.subplots_adjust(top=0.88,
# 					bottom=0.153,
# 					left=0.101,
# 					right=0.976,
# 					hspace=0.116,
# 					wspace=0.2)



fig, axs = plt.subplots(2,1, figsize=(6, 3.5), dpi=dpi, sharex=True)
fig.suptitle( 'On-Ground Net Impulse & Nozzle Length vs. Expansion Ratio', y=0.98)
axs[0].set_title(r'({} at $P_0$={} kPa, $P_{{amb}}$={} kPa, $V_{{p}}=${} cm$^3$, Nozzle $\varnothing =${} mm)'.format(gas_label, round(P_t_init/1000, 1), round(P_amb/1000, 1), vol*10**6, d_star*1000), fontsize=8)
fig.canvas.set_window_title('Sim Imp+NozzleLength_vs_ExpRatio {}'.format(gas_type))

axs[0].plot(list_of_expansion_ratios, [x*1000 for x in list_of_cumulative_impulses], color='#1f77b4', linestyle='-')
axs[1].plot(list_of_expansion_ratios, [x*1000 for x in list_of_nozzle_lengths], color='#1f77b4', linestyle='-')

for ax in axs:
	ax.yaxis.set_major_formatter(yfmt)
	ax.yaxis.offsetText.set_fontsize(7)
	ax.tick_params(axis='x', colors='#413839', labelsize=7, pad=0)
	ax.tick_params(axis='y', labelsize=7, pad=0)
	ax.grid(which='major', axis='both', linestyle='--')
	ax.xaxis.label.set_size(8)
	ax.yaxis.label.set_size(8)

axs[0].set(ylabel=r'Net Impulse, $mN-s$')
axs[1].set(ylabel=r'Nozzle Length, $mm$')
axs[-1].set(xlabel=r'Expansion Ratio, $\epsilon$')

plt.tight_layout()
plt.subplots_adjust(top=0.85)
fig.align_ylabels()

plt.show()