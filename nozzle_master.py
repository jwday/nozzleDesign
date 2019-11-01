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


## ---- OPTIONS --------------------------------------------------------------------

# Gas initial conditions
P_t_init = 114.7  # Max Total Pressure (psia)
P_amb = 14.7  # Ambient Pressure (psia)
T_t_init = 15 + 273.15  # Total Temperature (K)
vol = 35/10**6  # Plenum volume, units of m^3 (cm^3 / 10^6)


# no_of_points = 30
time_step = 0.001

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

density_init = (P_t_init*(137.9/20)*1000)/(R*(T_t_init))  # Units of kg/m^3 (or g/l)
m_init = density_init*(vol)  # Units of kg


## ---- DO THE THING ----------------------------------------------------------------

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

time = [0]
m_gas = [m_init]

# 0. Start with init P, T, m_gas
# 1. Run nozzle given P, T
# 2. Return m_dot, use to update m_gas assuming a specific time step
# 3. Use m_gas to determine new density
# 4. Use new density to update P, T assuming polytropic process + isentropic + ideal gas
# 5. Repeat 1-4 until P < 35

# for P_t in list_of_P_ts:
while list_of_P_ts[-1] > P_amb:
	m_dot, M_crit_sub, M_crit_sup, PR_crit_sub, PR_crit_sup, PR_exit_shock, M_exit_behindshock, M_exit, P_exit, v_exit, F, P_star, T_star, rho_star, T_exit, rho_exit = nozzle(list_of_P_ts[-1], list_of_T_ts[-1], P_amb, d_star, expansion_ratio, half_angle, gas_type)

	list_of_mdots.append(m_dot*1000)  # Units of g/s
	list_of_P_exits.append(P_exit)
	list_of_v_exits.append(v_exit)
	list_of_M_exits.append(M_exit)
	list_of_thrusts.append(F)
	list_of_pressure_ratios.append(P_exit/list_of_P_ts[-1])
	list_of_P_stars.append(P_star)
	list_of_T_stars.append(T_star)
	list_of_rho_stars.append(rho_star)
	list_of_T_exits.append(T_exit)
	list_of_rho_exits.append(rho_exit)

	time.append(time[-1] + time_step)
	m_gas.append(m_gas[-1] - m_dot*time_step)
	list_of_chamber_densities.append(m_gas[-1]/vol)
	list_of_T_ts.append( list_of_T_ts[-1]*(list_of_chamber_densities[-1]/list_of_chamber_densities[-2])**(k-1) )
	list_of_P_ts.append( list_of_P_ts[-1]*(list_of_chamber_densities[-1]/list_of_chamber_densities[-2])**k )


 # By the nature of this loop, anything that has an init value will end up with one extra element in its list
 # So we must manually remove the last element once all is said and done in order to make all the array lengths the same
 # This is due to the fact that the last element will be calculated for Pt < P_amb, which is not realistic.
del list_of_P_ts[-1], list_of_T_ts[-1], list_of_chamber_densities[-1], time[-1], m_gas[-1]

## ---- UPDATE THE THING ------------------------------------------------------------
# Search through the list of Pressure Ratios to find where the two critical conditions and the nozzle shock condition fit in
# Then at each of those Pressure Ratios, determine the Total Pressure and run it through nozzle() to get the exit conditions
# Then update the lists by inserting the exit conditions at the appropriate location

for i in range(len(list_of_pressure_ratios)):
	if ((1/PR_crit_sub) > list_of_pressure_ratios[i] and (1/PR_crit_sub) < list_of_pressure_ratios[i+1]) or ((1/PR_crit_sub) < list_of_pressure_ratios[i] and (1/PR_crit_sub) > list_of_pressure_ratios[i+1]):
		P_t = PR_crit_sub*P_amb

		# Use linear interpolation to estimate the temperature at which this event occurs
		m = (list_of_T_ts[i+1]-list_of_T_ts[i])/(list_of_P_ts[i+1]-list_of_P_ts[i])
		T_t = m*(P_t - list_of_P_ts[i]) + list_of_T_ts[i]

		# Use linear interpolation to estimate the time at which this event occurs
		m = (time[i+1]-time[i])/(list_of_P_ts[i+1]-list_of_P_ts[i])
		time_now = m*(P_t - list_of_P_ts[i]) + time[i]

		m_dot, M_crit_sub, M_crit_sup, PR_crit_sub, PR_crit_sup, PR_exit_shock, M_exit_behindshock, M_exit, P_exit, v_exit, F, P_star, T_star, rho_star, T_exit, rho_exit = nozzle(P_t, T_t, P_amb, d_star, expansion_ratio, half_angle, gas_type)

		list_of_P_ts.insert(i+1, P_t)
		list_of_T_ts.insert(i+1, T_t)
		time.insert(i+1, time_now)
		list_of_mdots.insert(i+1, m_dot*1000)
		list_of_P_exits.insert(i+1, P_exit)
		list_of_pressure_ratios.insert(i+1, P_exit/P_t)
		list_of_P_stars.append(P_star)
		list_of_M_exits.insert(i+1, M_exit)
		list_of_v_exits.insert(i+1, v_exit)
		list_of_thrusts.insert(i+1, F)
		list_of_T_stars.append(T_star)
		list_of_rho_stars.append(rho_star)
		list_of_T_exits.append(T_exit)
		list_of_rho_exits.append(rho_exit)
		print("Flow conditions calculated at subsonic critical pressure ratio!")
		break
	else:
		pass

for i in range(len(list_of_pressure_ratios)):
	if (1/PR_exit_shock > list_of_pressure_ratios[i] and 1/PR_exit_shock < list_of_pressure_ratios[i+1]) or (1/PR_exit_shock < list_of_pressure_ratios[i] and 1/PR_exit_shock > list_of_pressure_ratios[i+1]):
		P_t = PR_exit_shock*P_amb
		# Use linear interpolation to estimate the temperature at which this event occurs
		m = (list_of_T_ts[i+1]-list_of_T_ts[i])/(list_of_P_ts[i+1]-list_of_P_ts[i])
		T_t = m*(P_t - list_of_P_ts[i]) + list_of_T_ts[i]

		# Use linear interpolation to estimate the time at which this event occurs
		m = (time[i+1]-time[i])/(list_of_P_ts[i+1]-list_of_P_ts[i])
		time_now = m*(P_t - list_of_P_ts[i]) + time[i]

		m_dot, M_crit_sub, M_crit_sup, PR_crit_sub, PR_crit_sup, PR_exit_shock, M_exit_behindshock, M_exit, P_exit, v_exit, F, P_star, T_star, rho_star, T_exit, rho_exit = nozzle(P_t, T_t, P_amb, d_star, expansion_ratio, half_angle, gas_type)

		list_of_P_ts.insert(i+1, P_t)
		list_of_T_ts.insert(i+1, T_t)
		time.insert(i+1, time_now)
		list_of_mdots.insert(i+1, m_dot*1000)
		list_of_P_exits.insert(i+1, P_exit)
		list_of_pressure_ratios.insert(i+1, P_exit/P_t)
		list_of_P_stars.append(P_star)
		list_of_M_exits.insert(i+1, M_exit)
		list_of_v_exits.insert(i+1, v_exit)
		list_of_thrusts.insert(i+1, F)
		list_of_T_stars.append(T_star)
		list_of_rho_stars.append(rho_star)
		list_of_T_exits.append(T_exit)
		list_of_rho_exits.append(rho_exit)
		print("Flow conditions calculated at shock-at-exit pressure ratio!")
		break
	else:
		pass

for i in range(len(list_of_P_exits)):
	if (list_of_P_exits[i] > P_amb and list_of_P_exits[i+1] < P_amb) or (list_of_P_exits[i] < P_amb and list_of_P_exits[i+1] > P_amb):
		P_t = PR_crit_sup*P_amb
		# Use linear interpolation to estimate the temperature at which this event occurs
		m = (list_of_T_ts[i+1]-list_of_T_ts[i])/(list_of_P_ts[i+1]-list_of_P_ts[i])
		T_t = m*(P_t - list_of_P_ts[i]) + list_of_T_ts[i]

		# Use linear interpolation to estimate the time at which this event occurs
		m = (time[i+1]-time[i])/(list_of_P_ts[i+1]-list_of_P_ts[i])
		time_now = m*(P_t - list_of_P_ts[i]) + time[i]

		m_dot, M_crit_sub, M_crit_sup, PR_crit_sub, PR_crit_sup, PR_exit_shock, M_exit_behindshock, M_exit, P_exit, v_exit, F, P_star, T_star, rho_star, T_exit, rho_exit = nozzle(P_t, T_t, P_amb, d_star, expansion_ratio, half_angle, gas_type)

		list_of_P_ts.insert(i+1, P_t)
		list_of_T_ts.insert(i+1, T_t)
		time.insert(i+1, time_now)
		list_of_mdots.insert(i+1, m_dot*1000)
		list_of_P_exits.insert(i+1, P_exit)
		list_of_pressure_ratios.insert(i+1, P_exit/P_t)
		list_of_P_stars.append(P_star)
		list_of_M_exits.insert(i+1, M_exit)
		list_of_v_exits.insert(i+1, v_exit)
		list_of_thrusts.insert(i+1, F)
		list_of_T_stars.append(T_star)
		list_of_rho_stars.append(rho_star)
		list_of_T_exits.append(T_exit)
		list_of_rho_exits.append(rho_exit)
		print("Flow conditions calculated at supersonic critical pressure ratio!")
		break
	else:
		pass



## ---- CALC SOME OTHER STUFF ------------------------------------------------------
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


dia = 2*(vol*(3/4)/math.pi)**(1/3)  # Plenum diatmer , units of m



## ---- OUTPUT THE STUFF -----------------------------------------------------------

list_of_NaNs = [np.NaN for i in range(len(list_of_P_ts)-1)]

list_of_M_crit_subs = [M_crit_sub]
list_of_M_crit_sups = [M_crit_sup]
list_of_PR_crit_subs = [1/PR_crit_sub]
list_of_PR_crit_sups = [1/PR_crit_sup]
list_of_PR_exit_shocks = [1/PR_exit_shock]
list_of_M_exit_behindshock = [M_exit_behindshock]

for list in [list_of_M_crit_subs, list_of_M_crit_sups, list_of_PR_crit_subs, list_of_PR_crit_sups, list_of_PR_exit_shocks, list_of_M_exit_behindshock]:
        list.extend(list_of_NaNs)

out =  {'Time (s)': time,
		'Chamber Pressure (psi)': list_of_P_ts,
		'Chamber Temp (K)': list_of_T_ts,
		'Mass Flow Rate (g/s)': list_of_mdots,
		'Exit Pressure (psi)': list_of_P_exits,
		'Pressure Ratio': list_of_pressure_ratios,
		'Exit Mach Number': list_of_M_exits,
		'Exit Velocity (m/s)': list_of_v_exits,
		'Instantaneous Thrust (N)': list_of_thrusts,
		'': [np.NaN for i in range(len(list_of_P_ts))],
		'Critical Mach, Sub': list_of_M_crit_subs,
		'Critical Mach, Sup': list_of_M_crit_sups,
		'Critical PR, Sub': list_of_PR_crit_subs,
		'Critical PR, Sup': list_of_PR_crit_subs,
		'': [np.NaN for i in range(len(list_of_P_ts))],
		'PR for Shock at Exit': list_of_PR_exit_shocks,
		'Mach Behind Exit Shock': list_of_M_exit_behindshock}

df = pd.DataFrame.from_dict(out)
df.to_csv("nozzle_pressure.csv", index=False)


## ---- PLOT THE STUFF -----------------------------------------------------------

# ---- 1) Plot Inlet Pressures and Exit Velocity
fig1, ax1 = plt.subplots(figsize=(8, 5), dpi=90)

ax1.set_xlabel('Time (s)', color='#413839')
ax1.set_ylabel('Pressure (psia)', color='#413839')
ax1.plot(time, list_of_P_ts, color='#1f77b4', label='Inlet Pressure (psia)')
ax1.plot(time, list_of_P_exits, color='#ff7f0e', label='Exit Pressure (psia)')
ax1.tick_params(colors='#413839')
# ax1.set_ylim(0, 90)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

ax2.set_ylabel('Velocity (m/s)', color='#413839')  # we already handled the x-label with ax1
ax2.plot(time, list_of_v_exits, color='#2ca02c', label='Exit Velocity (m/s)')
# ax2.set_ylim(0, 350)

box = ax1.get_position()
ax1.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
fig1.legend(['Inlet Pressure (psia)', 'Exit Pressure (psia)', 'Exit Velocity (m/s)'], loc='center', bbox_to_anchor=(0.5, 0.03), ncol=3, frameon=False )

ax1.grid(which='major', axis='both', linestyle='--')

plt.title('Chamber + Exit Pres. & Exit Velocity for Single Plenum Discharge', y=1.03, color='#413839')


# ---- 2) Plot Instantaneous Thrust, Average Thrust, Total Impulse, and ISP
fig2, ax1 = plt.subplots(figsize=(8, 5), dpi=90)
ax1.set_xlabel('Time (s)', color='#413839')
ax1.set_ylabel('Thrust (N)', color='#413839')
ax1.plot(time, list_of_thrusts, color='#1f77b4', label='Instantaneous Thrust (N)')
ax1.plot(time_offset, average_thrust, color='#1f77b4', label='Cumulative Avg Thrust (N)', linestyle='--')
ax1.tick_params(colors='#413839')
# ax1.set_ylim(0, 0.5)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

ax2.set_ylabel('Total Impulse (N-s)', color='#413839')  # we already handled the x-label with ax1
ax2.plot(time_offset, cumulative_impulse, color='#2ca02c', label='Total Impulse (N-s)')
# ax2.set_ylim(0, 0.1)

ax3 = ax1.twinx()
ax3.set_ylabel('ISP', color='#413839')
ax3.plot(time, ISP, color='#ff7f0e', label='ISP')
ax3.spines['right'].set_position(('outward', 60))      

box = ax1.get_position()
ax1.set_position([box.x0, box.y0 + box.height * 0.1, box.width*0.88, box.height * 0.9])
fig2.legend(['Instantaneous Thrust (N)', 'Cumulative Avg Thrust (N)', 'Total Impulse (N-s)', 'ISP'], loc='center', bbox_to_anchor=(0.5, 0.03), ncol=4, frameon=False )

ax1.grid(which='major', axis='both', linestyle='--')

plt.title('Thrust & Total Impulse for Single Plenum Discharge', y=1.03, color='#413839')


# ---- Plot Mass Flow Rate and Total Mass
fig3, ax1 = plt.subplots(figsize=(8, 5), dpi=90)
ax1.set_xlabel('Time (s)', color='#413839')
ax1.set_ylabel('Mass Flow Rate (g/s)', color='#413839')
ax1.plot(time, list_of_mdots, color='#1f77b4', label='Mass Flow Rate (g/s)')
ax1.tick_params(colors='#413839')
# ax1.set_ylim(0, 1.3)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

ax2.set_ylabel('Propellant Consumed (g)', color='#413839')  # we already handled the x-label with ax1
ax2.plot(time_offset, cumulative_mass, color='#2ca02c', label='Propellant Consumed (g)')
# ax2.set_ylim(0, 0.3)

box = ax1.get_position()
ax1.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
fig3.legend(['Mass Flow Rate (g/s)', 'Propellant Consumed (g)'], loc='center', bbox_to_anchor=(0.5, 0.03), ncol=2, frameon=False )

ax1.grid(which='major', axis='both', linestyle='--')

plt.title('Mass Flow Rate & Total Propellant Mass Consumed for Single Plenum Discharge', y=1.03, color='#413839')


# ---- Plot Critical Pressure, Critical Temperature, and Critical Density
fig4, ax1 = plt.subplots(figsize=(8.5, 5), dpi=90)
ax1.set_xlabel('Time (s)', color='#413839')
ax1.set_ylabel('Temperature (K)', color='#413839')
ax1.tick_params(colors='#413839')
ax1.plot(time, list_of_T_stars, color='#1f77b4', label='Critical Temperature (K)')
# ax1.plot(time, list_of_T_exits, color='#1f77b4', label='Exit Temperature (K)')

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
ax2.set_ylabel('Density (g/mL)', color='#413839')  # we already handled the x-label with ax1
ax2.plot(time, [x/1000 for x in list_of_rho_stars], color='#2ca02c', label='Critical Density (g/mL)')
# ax2.plot(time, list_of_rho_exits, color='#2ca02c', label='Exit Density (g/mL)')

ax3 = ax1.twinx()
ax3.set_ylabel('Pressure (psi)', color='#413839')
ax3.plot(time, list_of_P_stars, color='#ed9926', label='Critical Pressure (psi)')
ax3.spines['right'].set_position(('outward', 60))      

box = ax1.get_position()
ax1.set_position([box.x0, box.y0 + box.height * 0.1, box.width*0.9, box.height * 0.9])
fig4.legend(['Critical Temperature (K)', 'Critical Density (g/mL)', 'Critical Pressure (psi)'], loc='center', bbox_to_anchor=(0.5, 0.03), ncol=3, frameon=False )

ax1.grid(which='major', axis='both', linestyle='--')

plt.title('Temperature and Density for Single Plenum Discharge', y=1.03, color='#413839')



# ---- Plot 2x Thrust(s), 2x Pressure(s), Impulse
fig5, ax1 = plt.subplots(figsize=(8.5, 5), dpi=90)
ax1.set_xlabel('Time (s)', color='#413839')
ax1.set_ylabel('Pressure (psia)', color='#413839')
ax1.plot(time, list_of_P_ts, color='#1f77b4', label='Inlet Pressure (psia)')
ax1.plot(time, list_of_P_exits, color='#1f77b4', label='Exit Pressure (psia)', linestyle='--')
ax1.tick_params(colors='#413839')

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

ax2.set_ylabel('Thrust (N)', color='#413839')
ax2.plot(time, list_of_thrusts, color='#ff7f0e', label='Instantaneous Thrust (N)')
ax2.plot(time_offset, average_thrust, color='#ff7f0e', label='Cumulative Avg Thrust (N)', linestyle='--')
ax2.tick_params(colors='#413839')

ax3 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

ax3.set_ylabel('Total Impulse (N-s)', color='#413839')  # we already handled the x-label with ax1
ax3.plot(time_offset, cumulative_impulse, color='#2ca02c', label='Total Impulse (N-s)')

ax3.spines['right'].set_position(('outward', 60))

box = ax1.get_position()
ax1.set_position([box.x0, box.y0 + box.height*0.1, box.width*0.88, box.height*0.88])
fig5.legend(['Inlet Pressure', 'Exit Pressure', 'Instantaneous Thrust', 'Cumulative Avg Thrust', 'Total Impulse'], loc='center', bbox_to_anchor=(0.5, 0.05), ncol=3, frameon=False )

ax1.grid(which='major', axis='both', linestyle='--')
plt.title('General Thruster Characteristics ({} mm)'.format(d_star), y=1.03, color='#413839')


# Inlet Plot Pressure vs. Thrust
fig6, ax1 = plt.subplots(figsize=(8.5, 5), dpi=90)
ax1.set_xlabel('Pressure (psia)', color='#413839')
ax1.set_ylabel('Thrust (N)', color='#413839')
ax1.plot(list_of_P_ts, list_of_thrusts, color='#1f77b4', label='Thrust (N)')
ax1.tick_params(colors='#413839')

ax1.grid(which='major', axis='both', linestyle='--')
plt.title('Thrust vs. Inlet Pressure ({} mm)'.format(d_star), y=1.03, color='#413839')



# Mass Flow Rate vs. Pressure (Isentropic vs. Experimental)
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



# data1.insert(4, "Nominal M used (g)", [round(x,1) for x in data1["M init (g)"]-data1["M final (g)"]], True)
# data1.insert(6, "Nominal Mdot (g/s)", [round(x, 4) for x in (data1["Nominal M used (g)"])/(data1["Time (s)"])])
# data1.insert(7, "Max Mdot (g/s)", [round(x, 4) for x in (data1["Nominal M used (g)"]+0.13)/(data1["Time (s)"]-0.0001)])
# data1.insert(8, "Min Mdot (g/s)", [round(x, 4) for x in (data1["Nominal M used (g)"]-0.13)/(data1["Time (s)"]+0.0001)])
# data1.insert(9, "Mdot Error (g/s)", [round(x, 4) for x in data1["Max Mdot (g/s)"]-data1["Min Mdot (g/s)"]])

# Setup theoretical vs. experimental plot
fig7, ax1 = plt.subplots(figsize=(8.5, 5), dpi=90)
    # Blue: #1f77b4
    # Orange: #ff7f0e
    # Green: #2ca02c
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

fig7.legend(['Inviscid', 'Experimental', 'Temperature'], loc='center', bbox_to_anchor=(0.5, 0.03), ncol=3, frameon=False )
plt.title('Mass Flow Rate Comparison ({} mm)'.format(d_star), y=1.03, color='#413839')



# ---- Plot Everything
plt.show()
