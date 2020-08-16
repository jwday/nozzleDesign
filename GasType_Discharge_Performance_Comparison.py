# Nozzle master control
# This code is only for plotting nozzle and impulse for Chapter 2 of thesis. Nothing else.
from nozzle_code_3 import nozzle
import matplotlib.pyplot as plt
from scipy import stats
import pandas as pd
import numpy as np
import seaborn as sns
from data_handling_funcs import *
import seaborn as sns
import matplotlib as mpl
from scipy.interpolate import interp1d
from scipy.interpolate import interp2d
import scipy.optimize as opti

## ==================================================================================
## ---- USER OPTIONS ----------------------------------------------------------------
## ==================================================================================

gas_types = ['R134a', 'CO2']		# Gas choices: R236fa, R134a, N2, CO2, H2, air
								# Gas choice will determine geometry based on desired output that was determined in the course of this project
cutoff_cond = 0.0001			# Cutoff condition, defined by the fractional change in pressure (relative to P_t_init) per second, units of 1/sec
figsize = (7.5, 4)				# Figure size (in)
dpi = 150						# Figure dpi



## ==================================================================================
## ---- SETUP LOOP ------------------------------------------------------------------
## ==================================================================================

## --------------------------------------------------------------------------------
# Set up Pandas DataFrame to store aggregate data
all_data = pd.DataFrame(columns = [	'time',
									'P_t',
									'T_t',
									'rho_t',
									'm_gas',
									'mdot',
									'P_exit',
									'v_exit',
									'M_exit',
									'thrust',
									'P_star',
									'T_star',
									'rho_star',
									'Re_star',
									'T_exit',
									'rho_exit',
									'thrust_coeff',
									'visc_loss',
									'avg_thrust',
									'cum_impulse',
									'ISP',
									'gas_type'])


## --------------------------------------------------------------------------------
# Create a 2D function to return viscosity of a fluid at a given (P,T) based on NIST data and using linear interpolation
# Currently not working
def create_visc_func(data_table_loc):
	# fluid_props = pd.read_excel(data_table_loc, sheet_name=None)																			# .xlxs file containing the fluid properties at variable T for fixed P per sheet
	# visc_fluid = pd.DataFrame(columns=['Temperature (K)'])
	# visc_pres = []
	# for i in list(fluid_props.keys()):																										# Iterate over sheet names aka Fixed Pressures
	# 	data = fluid_props[i][['Temperature (K)', 'Viscosity (uPa*s)']]																		# Acquire Temp+Visc data from the current pressure -- BUT WHAT IF A PHASE CHANGE OCCURS IN THIS DATA!?
	# 	# Now you gotta handle what to do with phase changes, because if one exists, then you'll have two entries for the same P-T but different viscosity
	# 	# Also, those Temps are unique so there won't be viscosity data at other pressures for those temps, so you'll have to interpolate, but that should be easy (use df.interpolate())
	# 	visc_fluid = pd.merge(visc_fluid, data, how='outer', on='Temperature (K)')															# Merge the Temp+Visc data for the current pressure into visc_fluid. 'Outer' join will add a row for non-shared temperatures.
	# 	visc_fluid.rename(columns={'Viscosity (uPa*s)': i}, inplace=True)																	# Rename the viscosity column to correspond to the pressure, labeled as 'XX MPa'
	# 	visc_pres.append(fluid_props[i]['Pressure (MPa)'][0])																				# Append to a list of numerical values of pressure
	# visc_fluid = visc_fluid.sort_values('Temperature (K)').reset_index(drop=True)
	# visc_func = interp2d([x*1E6 for x in visc_pres], visc_fluid['Temperature (K)'].values, visc_fluid.iloc[:,1:].values)

	fluid_props = pd.read_excel(data_table_loc, sheet_name=None)	# ALL the fluid properties 
	visc_fluid = pd.DataFrame()
	temps = fluid_props[list(fluid_props.keys())[0]]['Temperature (K)']
	visc_fluid = pd.concat([visc_fluid, temps], axis=1)
	pressures = []
	for pres in list(fluid_props.keys()):	# Iterate over sheet names
		f = fluid_props[pres]['Viscosity (uPa*s)'].rename(pres, inplace=True)
		visc_fluid = pd.concat([visc_fluid, f], axis=1)
		pressures.append(fluid_props[pres]['Pressure (MPa)'][0]*1E6)
	visc_func = interp2d(pressures, visc_fluid['Temperature (K)'].values, visc_fluid.iloc[:,1:].values)
	return visc_func


## --------------------------------------------------------------------------------
# Create two 1D functions to return the saturated pressure given a temperature, and another to return the saturated temperature given pressure
def create_phase_funcs(data_table_loc):
	fluid_props = pd.read_excel(data_table_loc, sheet_name=None)
	phase_change = pd.DataFrame()
	saturated_temp = []
	saturated_pres = []
	for pres in list(fluid_props.keys()):
		fluid_props[pres]['Phase Change'] = fluid_props[pres]['Phase'].ne(fluid_props[pres]['Phase'].shift().bfill()).astype(bool)		# Add a column called 'Phase Change' to identify when a phase change takes place
		if fluid_props[pres]['Phase Change'].any():
			saturated_temp.append(fluid_props[pres][fluid_props[pres]['Phase Change'] == True]['Temperature (K)'].values[0])			# Returns float, is temperature of phase change at specified pressure
			saturated_pres.append(fluid_props[pres]['Pressure (MPa)'][0]*1E6)															# Returns float, is the specified pressure
		else:
			saturated_temp.append(np.nan)
			saturated_pres.append(np.nan)
	log_saturated_pres = [np.log(x) for x in saturated_pres]																			# Behavior can be modeled and interpolated logarithmically (T vs log(P) is very linear)
	f = interp1d(pd.DataFrame(saturated_temp).dropna()[0].values, pd.DataFrame(log_saturated_pres).dropna()[0].values, kind='quadratic', fill_value='extrapolate')		# Linearlly interpolate over T vs log(P), return a function f. Need to convert to df to drop na's to use quadratic fit
	g = interp1d(pd.DataFrame(log_saturated_pres).dropna()[0].values, pd.DataFrame(saturated_temp).dropna()[0].values, kind='quadratic', fill_value='extrapolate')		# Linearlly interpolate over T vs log(P), return a function g. Need to convert to df to drop na's to use quadratic fit

	def saturated_pres_from_temp(temp):
		h = np.e**f(temp)																												# When a temp is specified, use the function f to return a log(P), then exponentiate it to return P
		return h

	def saturated_temp_from_pres(pres):
		h = float(g(np.log(pres)))																										# When a temp is specified, supply the function g with a log(P) to return a T. g returns a 0D array -- convert to float
		return h
		
	return saturated_pres_from_temp, saturated_temp_from_pres


## ==================================================================================
## ---- BEGIN LOOP ------------------------------------------------------------------
## ==================================================================================

## --------------------------------------------------------------------------------
# Select gas properties
for gas_type in gas_types:
	if gas_type == 'CO2':
		P_t_init = 114.7 * 6894.76  	# Init Total Pressure, units of Pa (psia * 6894.76)
		P_amb = 14.7 * 6894.76  		# Ambient Pressure, units of Pa (psia * 6894.76)
		T_t_init = 20 + 273.15  			# Init Total Temperature, units of K (C + 273.15)
		vol = 30 / 10**6  				# Plenum volume, units of m^3 (cm^3 / 10^6)
		d_star = 0.6 / 1000  			# Nozzle throat diameter, units of m (mm / 1000)
		half_angle = 10  				# (Conical) Nozzle expansion angle (degrees)
		expansion_ratio = 1.17			# Nozzle expansion ratio (Exit Area / Throat Area)
		gas_label = 'CO$_{2}$'
		right_limit = 1.4*((vol * 10**6)/30)*(d_star * 1000)/(0.6) # For X-axis time scale
		visc_loss_param = 0.39			# As determined from NASA TN D-3056 (This stays pretty constant)
		k = 1.289
		R = 8.314/0.04401  				# Specific gas constant (J/kg-K)
		T_trip = 216.58  				# Triple point temperature (K)
		P_trip = 518500  				# Triple point pressure (Pa)
		fg_pres_from_temp, fg_temp_from_pres = create_phase_funcs('CO2_props_NIST.xlsx')

		# Create a 2D function to return viscosity of CO2 at a given Temp + Pres, based on NIST data and using linear interpolation
		# visc_CO2 = pd.read_csv('../DESKTOP/CO2_visc_PvT.csv')	# Viscosity in Pa*s
		# visc_CO2.iloc[:,1:] = visc_CO2.iloc[:,1:].mul(1E6) # Change all Viscosity in Pa*s to Viscosity in uPa*s
		# visc_func = interp2d([x*1E6 for x in [0.8,0.6,0.4,0.2,0.1]], visc_CO2['Temp (K)'].values, visc_CO2.iloc[:,1:].values)
		visc_func = create_visc_func('CO2_props_NIST.xlsx')

	elif gas_type == 'R134a':
		P_t_init = 82.9 * 6894.76  		# Init Total Pressure, units of Pa (psia * 6894.76)
		P_amb = 0 * 6894.76  			# Ambient Pressure, units of Pa (psia * 6894.76)
		T_t_init = 20 + 273.15  		# Init Total Temperature, units of K (C + 273.15)
		vol = 11.2 / 10**6  			# Plenum volume, units of m^3 (cm^3 / 10^6)
		d_star = 0.212 / 1000  			# Nozzle throat diameter, units of m (mm / 1000)
		half_angle = 10  				# (Conical) Nozzle expansion angle (degrees)
		expansion_ratio = 30			# Nozzle expansion ratio (Exit Area / Throat Area)
		gas_label = 'R134a'
		right_limit = 28*((vol * 10**6)/11.2)*(d_star * 1000)/(0.212) # For X-axis time scale
		visc_loss_param = 4.0			# As determined from NASA TN D-3056 (This varies a lot over Re numbers)
		k = 1.127  						# Override value to compare results to MST paper
		R = 8.314/0.10203  				# Specific gas constant (J/kg-K)
		T_trip = 169.85  				# Triple point temperature (K)
		P_trip = 389.56  				# Triple point pressure (Pa)
		fg_pres_from_temp, fg_temp_from_pres = create_phase_funcs('R134a_props_NIST.xlsx')

		# visc_func = create_visc_func('../DESKTOP/R134a_props_NIST.xlsx')

		# temps = fluid_props[list(fluid_props.keys())[0]]['Temperature (K)'].rename('Temp (K)', inplace=True)
		# pres = [fluid_props[x]['Pressure (MPa)'][0]*1E6 for x in list(fluid_props.keys())[0]]
		# visc_fluid = pd.concat([visc_fluid, temps], axis=1)

		# temps = pd.DataFrame(columns=['Temperature (K)'])
		# for i in list(fluid_props.keys()):	# Iterate over sheet names aka Fixed Pressures
		# 	foo = fluid_props[i]['Viscosity (uPa*s)'].rename(i, inplace=True)
		# 	visc_fluid = pd.concat([visc_fluid, foo], axis=1)
		# visc_func = interp2d([x*1E6 for x in [0.6,0.4,0.2,0.1]], visc_fluid['Temp (K)'].values, visc_fluid.iloc[:,1:].values)


		# Create a 2D function to return viscosity of R134a at a given Temp + Pres, based on NIST data and using linear interpolation
		# R134a_props = pd.read_excel('../DESKTOP/R134a_props_NIST.xlsx', sheet_name=None)	# ALL the R134a properties 
		# visc_R134a = pd.DataFrame()
		# temps = R134a_props[list(R134a_props.keys())[0]]['Temperature (K)'].rename('Temp (K)', inplace=True)
		# visc_R134a = pd.concat([visc_R134a, temps], axis=1)
		# for i in list(R134a_props.keys()):	# Iterate over sheet names
		# 	temp = R134a_props[i]['Viscosity (uPa*s)'].rename(i, inplace=True)
		# 	visc_R134a = pd.concat([visc_R134a, temp], axis=1)
		# visc_func = interp2d([x*1E6 for x in [0.8,0.6,0.4,0.2,0.1,0.01]], visc_R134a['Temp (K)'].values, visc_R134a.iloc[:,1:].values)

		visc_func = create_visc_func('R134a_props_NIST.xlsx')





	# ----- For future implementation-----
	# Combine the gas properties into a dictionary and select from there, rather than using the 'if' statements as above
	# gas_props = {
	# 	'P_t_init' : [x*6894.76 for x in [114.7, 82.9]]  	# Init Total Pressure, units of Pa (psia * 6894.76)
	# 	'P_amb' : [x*6894.76 for x in [14.7, 0]]  			# Ambient Pressure, units of Pa (psia * 6894.76)
	# 	'T_t_init' : [x+273.15 for x in [0, 20]]  			# Init Total Temperature, units of K (C + 273.15)
	# 	'vol' : [x/10**6 for x in [30, 11.2]]  				# Plenum volume, units of m^3 (cm^3 / 10^6)
	# 	'd_star' : [x/1000 for x in [0.6, 0.212]]  			# Nozzle throat diameter, units of m (mm / 1000)
	# 	'half_angle' : [10, 10]  							# (Conical) Nozzle expansion angle (degrees)
	# 	'expansion_ratio' : [1.17, 30]						# Nozzle expansion ratio (Exit Area / Throat Area)
	# 	'gas_label' : ['CO$_{2}$', 'R134a']
	# 	'visc_loss_param' : [0.39		# As determined from NASA TN D-3056
	# 	'k' : [1.289
	# 	'R' : [8.314/0.04401  			# Specific gas constant (J/kg-K)
	# 	'T_trip' : [216.58  				# Triple point temperature (K)
	# 	'P_trip' : [518500  				# Triple point pressure (Pa)
	# }
	# right_limit = [x * ((vol*10**6)/gas_props['vol'])*(d_star*1000)/(gas_props['d_star']) for x in [1.4, 28]] # For X-axis time scale

	rho_t_init = P_t_init/(R*(T_t_init))  # Units of kg/m^3 (or g/l)
	m_init = rho_t_init*vol  # Units of kg
	dia = 2*(vol*(3/4)/np.pi)**(1/3)  # Plenum diatmer , units of m
	time_step_init = 15.28*vol*(P_t_init-P_amb)/((P_t_init*d_star)**2)
	if gas_type == 'CO2':
		time_step_init = 15.28*vol*(P_t_init-P_amb)/((P_t_init*d_star)**2)



	## ==================================================================================
	## ---- INIT DATA LISTS -------------------------------------------------------------
	## ==================================================================================
	time = [0]
	list_of_P_ts = [P_t_init] 
	list_of_T_ts = [T_t_init]
	list_of_rho_ts = [rho_t_init]
	m_gas = [m_init]

	list_of_mdots = []
	list_of_P_exits = []
	list_of_v_exits = []
	list_of_c_exits = []
	list_of_M_exits = []
	list_of_thrusts = []
	list_of_pressure_ratios = []
	list_of_P_stars = []
	list_of_T_stars = []
	list_of_rho_stars = []
	list_of_Re_stars = []
	list_of_v_stars = []
	list_of_T_exits = []
	list_of_rho_exits = []
	list_of_thrust_coeffs = []
	list_of_visc_losses = []
	list_of_F_mdotv = []
	list_of_F_pdiff = []

	average_thrust = []
	cumulative_impulse = []
	ISP = []
	cumulative_mass = []

	list_of_flow_regimes = [] 
	list_of_area_ratios_at_shock = []

	list_of_P_fg_exit = []
	list_of_T_fg_exit = []

	## ==================================================================================
	## ---- EXECUTE LOOP ----------------------------------------------------------------
	## ==================================================================================
	#-1. Solve Mach Number-Area relation ONCE to determine the CRITICAL MACH NUMBERS (applies only to fully subsupersonic and fully subsonic flows)
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
	delta_pres = 1

	## --------------------------------------------------------------------------------
	# Nozzle Geometry
	A_star = np.pi*(d_star**2)/4  					# Throat area
	A_exit = A_star*expansion_ratio  				# Exit area
	d_exit = np.sqrt(4*A_exit/np.pi)  				# Exit diameter (m)

	# Establish thermodynamic relationships in more compact terms
	P = (k+1)/2
	L = (k-1)/2
	W = k/(k-1)
	Q = P/L  # aka (k+1)/(k-1)
	S = (A_exit/A_star)**2
	Y = np.sqrt(k/R)
	
	def pres_ratio(M):
		f = (1 + L*M**2)**(-W)  # P/P_t, your basic Isentropic Mach-Pressure relation.
		return f

	## --------------------------------------------------------------------------------
	# Solve the Isentropic Mach Number-Area relation for the (critical) Mach Numbers at the EXIT (defined by S above)
	# NOTE: This equation is only valid for fully subsonic OR fully supersonic flow throughout the entire nozzle
	# It is not valid for the regime which includes a normal shock within the nozzle
	def objective(X, ARsq):
		f = (1 + L*X)**Q - ARsq*X*(P**Q)  # Where 'X' is M^2 and ARsq is AREA RATIO SQUARED (so that this same function can be used for different area ratios i.e. along the length of the nozzle)
		return f

	x0 = np.array([0.00001, 20])	# List of M^2 to solve for roots over
	sol0 = opti.fsolve(objective, x0, args=S, maxfev=100000, full_output=False, xtol=0.000000001)	# Solve for M^2 given S as an extra argument

	M_crit_sub = np.sqrt(sol0[0])  # Subsonic critical mach no.
	M_crit_sup = np.sqrt(sol0[1])  # Supersonic critical mach no.

	## --------------------------------------------------------------------------------
	# This will be used to help us estimate the viscous losses
	# I'm going to define a series of Mach vs. Area Ratio, then use some linear interpolation to determine a Mach number at a given area ratio (for viscous loss integral). Should only have to do this one.
	range_of_subsonic_mach_nos = np.linspace(0.1, 1.1, 501)
	range_of_supersonic_mach_nos = np.linspace(1.1, 5.1, 501)

	range_of_subsonic_area_ratios = [np.sqrt( ((1 + L*x**2)**Q) / ((x**2)*(P**Q)) ) for x in range_of_subsonic_mach_nos]
	range_of_supersonic_area_ratios = [np.sqrt( ((1 + L*x**2)**Q) / ((x**2)*(P**Q)) ) for x in range_of_supersonic_mach_nos]

	subsonic_mach_anywhere = interp1d(range_of_subsonic_area_ratios, range_of_subsonic_mach_nos)
	supersonic_mach_anywhere = interp1d(range_of_supersonic_area_ratios, range_of_supersonic_mach_nos)

	def stream_props(area_ratio, P_t, T_t):
		mach_no = supersonic_mach_anywhere(area_ratio)
		pres = P_t * (1 + L*mach_no**2)**-W
		temp = T_t * (1 + L*mach_no**2)**-1
		return mach_no, pres, temp

	## --------------------------------------------------------------------------------
	# Begin loop, starting with a P and T, and running the nozzle code
	# Nozzle code will return throat and exit flow properties (along with thrust, ISP, impulse) based on 
	while delta_pres > cutoff_cond and list_of_P_ts[-1] > P_amb:
		P_star, T_star, rho_star, Re_star, v_star, P_exit, T_exit, rho_exit, M_exit, v_exit, c_exit, m_dot, F, CF, flow_regime, area_ratio_at_shock, F_mdotv, F_pdiff = nozzle(k, R, M_crit_sub, M_crit_sup, list_of_P_ts[-1], list_of_T_ts[-1], list_of_rho_ts[-1], P_amb, d_star, expansion_ratio, half_angle, gas_type, visc_func)

		# Append returned items to their respective lists
		list_of_P_stars.append(P_star)		# Units of Pa
		list_of_T_stars.append(T_star)		# Units of K
		list_of_rho_stars.append(rho_star)	# Units of kg/m^s
		list_of_Re_stars.append(Re_star)	# Nondimensional
		list_of_v_stars.append(v_star)

		list_of_P_exits.append(P_exit)		# Units of Pa
		list_of_T_exits.append(T_exit)		# Units of K
		list_of_rho_exits.append(rho_exit)	# Units of kg/m^3
		list_of_M_exits.append(M_exit)		# Nondimensional
		list_of_v_exits.append(v_exit)		# Units of m/s
		list_of_c_exits.append(c_exit)		# Units of m/s

		list_of_mdots.append(m_dot*1000)  	# Units of g/s
		list_of_thrusts.append(F)			# Units of N
		list_of_thrust_coeffs.append(CF)	# Nondimensional
		list_of_flow_regimes.append(flow_regime)
		list_of_area_ratios_at_shock.append(area_ratio_at_shock)
		list_of_F_mdotv.append(F_mdotv)
		list_of_F_pdiff.append(F_pdiff)

		average_thrust.append( np.average(list_of_thrusts) )	# Cumulative average
		ISP.append( 1000*list_of_thrusts[-1]/(9.81*list_of_mdots[i]) )

		if i == 0:  # If we're on the first time step...
			cumulative_impulse.append( time_step*list_of_thrusts[-1])
		else:  # If we're on any other time step...
			cumulative_impulse.append(time_step*list_of_thrusts[-1] + cumulative_impulse[-1])
		
		## --------------------------------------------------------------------------------
		# Maybe now is when we can check if a shock exists? Or rather...
		# Nozzle code already checks that. It's based on pressure ratios. Just have the nozzle code return a boolean for "shock in flow"
		# If the flag is 'False', return 'n/a'
		# If the flag is 'True', determine where the shock is located and report it in terms of the Area Ratio?
		# Done

		## --------------------------------------------------------------------------------
		# Now we should try to see if we can determine the actual PHASE mixture of the flow AT THE EXIT ONLY (for now)
		# First, let's try to verify that the entropy is actually constant maybe? Well no, because the equations you're using assume an isentropic process (adiabatic and internally reversable), and the entropy change for an adiabatic process is 0.
		# Maybe instead we should just try to plot the phase change lines along with the actual pressure.
		# That means writing a function that goes through every sheet in the R134a (and later CO2) data and returning the (P,T) at which the phase change takes place
		# Done
		list_of_P_fg_exit.append(fg_pres_from_temp(T_exit))
		list_of_T_fg_exit.append(fg_temp_from_pres(P_exit))


		## --------------------------------------------------------------------------------
		# Now let's calculate viscous losses
		# This is my temporary function which is only really valid for CO2 but I'll use it for both for the time being
		list_of_visc_losses.append(visc_loss_param/np.sqrt(Re_star*np.tan(np.deg2rad(half_angle))))

		# Integrate over expansion ratio of nozzle (from 1 to either 1.17 or 30, depending on if CO2 or R134a)
		# Need Area and Mach number at each point
		# Establish your domain
		# 
		# First you need a function that will give you the LOCAL stream properties based on the given area ratio (aka where along the nozzle you are)
		# The result will depend on if you're in part of a supersonic regime, part of a fully subsonic regime, or if you've passed a shock
		# So first thing to do is check if a shock exists inside the nozzle
		#
		# Oh shiz though if there's a phase change then that's going to wildly change what the viscosity is. Case in point: R134a appears to undergo exactly that process.
		# So then you'll need to figure out how to handle phase changes in the nozzle. THAT'S going to be fun.
		



		## --------------------------------------------------------------------------------
		# Calculate these properties in preparation for the next loop. Loop will end if the newly calculated pressure drops below the ambient pressure.
		m_gas.append(m_gas[-1] - m_dot*time_step) 
		list_of_rho_ts.append(m_gas[-1]/vol)
		list_of_T_ts.append( list_of_T_ts[-1]*(list_of_rho_ts[-1]/list_of_rho_ts[-2])**(k-1) )
		list_of_P_ts.append( list_of_P_ts[-1]*(list_of_rho_ts[-1]/list_of_rho_ts[-2])**k )
		time.append((i+1)*time_step)  # The first iteration is at t=0, so the first time[] entry will be 0.

		# Loop counter to check cutoff condition
		if i>=2:
			delta_pres = np.absolute((list_of_P_ts[-2] - list_of_P_ts[-1])/P_t_init)/time_step
		i+=1
		print('Time step: ' + str(round(time_step, 6)) + ' sec, P_t: ' + str(round(list_of_P_ts[-1]/ 6894.76, 1)) + ' psia, Change in pres: ' + str(round(delta_pres*100, 3)) + '%/sec ' + str(round(time[-1], 4)) + ' sec', end='\r', flush=True)

	print('\n')

	# By the nature of this loop, anything that has an init value will end up with one extra element in its list
	# So we must manually remove the last element once all is said and done in order to make all the array lengths the same
	del m_gas[-1], list_of_rho_ts[-1], list_of_T_ts[-1], list_of_P_ts[-1], time[-1]

	## ==================================================================================
	## ---- ALL DATA SAVED TO PD DATAFRAME ----------------------------------------------
	## ==================================================================================

	current_data = pd.DataFrame(zip(time,
									list_of_P_ts, 
									list_of_T_ts, 
									list_of_rho_ts, 

									list_of_P_stars, 
									list_of_T_stars, 
									list_of_rho_stars, 
									list_of_Re_stars,
									list_of_v_stars,

									list_of_P_exits, 
									list_of_T_exits,
									list_of_rho_exits,
									list_of_M_exits, 
									list_of_v_exits, 
									list_of_c_exits, 

									m_gas,
									list_of_mdots,
									list_of_F_mdotv,
									list_of_F_pdiff,
									list_of_thrusts,
									list_of_thrust_coeffs,
									list_of_visc_losses,
									average_thrust,
									cumulative_impulse,
									ISP,
									list_of_flow_regimes,
									list_of_area_ratios_at_shock,

									list_of_P_fg_exit,
									list_of_T_fg_exit),

						columns = [	'time',
									'P_t',
									'T_t',
									'rho_t',

									'P_star',
									'T_star',
									'rho_star',
									'Re_star',
									'v_star',

									'P_exit',
									'T_exit',
									'rho_exit',
									'M_exit',
									'v_exit',
									'c_exit',

									'm_gas',
									'mdot',
									'F_mdotv',
									'F_pdiff',
									'thrust',
									'thrust_coeff',
									'visc_loss',
									'avg_thrust',
									'cum_impulse',
									'ISP',
									'flow regimes',
									'area ratios at shock',
									
									'P_fg_exit',
									'T_fg_exit'])
	current_data['gas_type'] = gas_type
	current_data['P_trip'] = P_trip
	current_data['T_trip'] = T_trip

	all_data = all_data.append(current_data, ignore_index=True)



## ==================================================================================
## ---- PLOT ------------------------------------------------------------------------
## ==================================================================================

# ---- Plot 2x3 [Thrust, Impulse, ISP, Re, Ma, Density] all vs. Inlet + Throat + Exit Pressure

linewidth = 2
fontsize = 12

data = 	{ 
			# 'P_t': 			all_data['P_t'],
			# 'T_t': 			all_data['T_t'],
			# 'rho_t':			all_data['rho_t'],

			# 'P_star': 		all_data['T_t'],
			# 'T_star': 		all_data['T_t'],
			# 'rho_star': 		all_data['rho_star'],
			# 'Re_star': 		all_data['Re_star'],
			# 'v_star': 		all_data['v_star'],

			'P_exit': 		all_data['P_exit'],
			'T_exit': 		all_data['T_exit'],
			# 'rho_exit': 		all_data['rho_exit'],
			# 'M_exit': 		all_data['M_exit'], 
			# 'v_exit': 		all_data['v_exit'],
			# 'c_exit': 		all_data['c_exit'],

			# 'm_gas': 			all_data['m_gas'],
			# 'mdot': 			all_data['mdot'],
			# 'F_mdotv': 		all_data['F_mdotv'], 
			# 'F_pdiff': 		all_data['F_pdiff'], 
			# 'thrust': 		all_data['thrust'],
			# 'thrust_coeff': 	all_data['thrust_coeff'],
			# 'visc_losses': 	all_data['visc_loss'],
			# 'avg_thrust': 	all_data['avg_thrust'],
			# 'cum_impulse': 	all_data['cum_impulse'],
			# 'ISP': 			all_data['ISP'],
			# 'ARs at shock':	all_data['area ratios at shock']

			# 'P_fg_exit':		all_data['P_fg_exit'],
			# 'T_fg_exit':		all_data['T_fg_exit']
			  }

figname = { 
			'P_t': 				'Total Pressure',
			'T_t': 				'Total Temperature', 							
			'rho_t': 			'Total Density',

			'P_star': 			'Throat Pressure',
			'T_star': 			'Throat Temperature',
			'rho_star': 		'Throat Density', 		
			'Re_star': 			'Throat Reynold\'s No.',
			'v_star': 			'Throat Velocity',

			'P_exit': 			'Exit Pressure',
			'T_exit': 			'Exit Temperature',
			'rho_exit': 		'Exit Density', 		
			'M_exit': 			'Exit Mach No.', 	
			'v_exit': 			'Exit Velocity', 	
			'c_exit': 			'Exit Speed of Sound', 	

			'm_gas': 			'Propellant Mass',
			'mdot': 			'Mass Flow Rate',
			'F_mdotv': 			'Thrust From mdot',
			'F_pdiff': 			'Thrust From pdiff',
			'thrust': 			'Instantaneous Thrust', 							
			'thrust_coeff': 	'Thrust Coefficient',
			'visc_losses': 		'Viscous Loss Percentage',
			'avg_thrust': 		'Time Average Thrust',
			'cum_impulse': 		'Net Impulse',
			'ISP': 				'$I_{SP}$',
			'ARs at shock':		'Area Ratio at Shock',

			'P_fg_exit':		'Saturated Vapor Pressure at Exit',
			'T_fg_exit':		'Saturated Vapor Temperature at Exit'
		  }

times = {
			'P_t': 				all_data['time'],
			'T_t': 				all_data['time'],
			'rho_t':			all_data['time'],

			'P_star': 			all_data['time'],
			'T_star': 			all_data['time'],
			'rho_star': 		all_data['time'],
			'Re_star': 			all_data['time'],
			'v_star': 			all_data['time'],

			'P_exit': 			all_data['time'],
			'T_exit': 			all_data['time'],
			'rho_exit': 		all_data['time'],
			'M_exit': 			all_data['time'],
			'v_exit': 			all_data['time'],
			'c_exit': 			all_data['time'],

			'm_gas': 			all_data['time'],
			'mdot': 			all_data['time'],
			'F_mdotv': 			all_data['time'],
			'F_pdiff': 			all_data['time'],
			'thrust': 			all_data['time'],
			'thrust_coeff': 	all_data['time'],
			'visc_losses': 		all_data['time'],
			'avg_thrust': 		all_data['time'],
			'cum_impulse': 		all_data['time'],
			'ISP': 				all_data['time'],
			'ARs at shock': 	all_data['time'],

			'P_fg_exit': 	all_data['time'],
			'T_fg_exit': 	all_data['time']
		}

ylabels = {
			'P_t': 				'Total Pressure, $Pa$',
			'T_t': 				'Total Thrust, $K$', 						
			'rho_t': 			'Total Density, $kg/m^3$',

			'P_star': 			'Throat Pressure, $Pa$',
			'T_star': 			'Throat Temperature, $K$',
			'rho_star': 		'Throat Density, $kg/m^3$', 		
			'Re_star': 			'Throat Reynold\'s No.',
			'v_star': 			'Throat Velocity, $m/s$', 				

			'P_exit': 			'Exit Pressure, $Pa$',
			'T_exit': 			'Exit Temperature, $K$',
			'rho_exit': 		'Exit Density, $kg/m^3$',
			'M_exit': 			'Exit Mach No.', 				
			'v_exit': 			'Exit Velocity, $m/s$', 				
			'c_exit': 			'Exit Speed of Sound, $m/s$', 				

			'm_gas':			'Propellant Mass, $g$',
			'mdot': 		   r'Mass Flow Rate ($\dot{m}$), $g/s$',
			'F_mdotv': 			'Thrust from mdot, $mN$', 
			'F_pdiff': 			'Thrust from pdiff, $mN$', 
			'thrust': 			'Thrust, $mN$',
			'thrust_coeff':		'Thrust Coefficient, $C_f$',
			'visc_losses': 	   r'$\frac{C_f - C_{f_v}}{C_f}$ (%)',
			'avg_thrust':		'Time Average Thrust, $mN$',
			'cum_impulse': 		'Impulse, $mN-s$',
			'ISP': 				'$I_{SP}$, $s$', 		
			'ARs at shock':	   r'Area Ratio, $\lambda$',

			'P_fg_exit':		'Saturated Pressure, $Pa$',
			'T_fg_exit':		'Saturated Temperature, $K$',
		   }

titles = { 	
			'CO2': 'On-Ground Single Plenum Discharge Reynolds Number',
			'R134a': 'In-Space Single Plenum Discharge Reynolds Number'
		 }


# Just some formatting stuff 
class ScalarFormatterForceFormat(mpl.ticker.ScalarFormatter):
		def _set_format(self):  # Override function that finds format to use.
			self.format = "%1.1f"  # Give format here

# Let's see if we can plot exit pres & sat pres at exit on same plot, and also temp on another
fig_sat, axs = plt.subplots(2,1, figsize=figsize, dpi=dpi, sharex=True)
fig_sat.suptitle('Saturated Pressure & Temperature On-ground', y=0.98)
fig_sat.canvas.set_window_title('Saturated Pressure Stuff')

for i,key in enumerate(data):
	sns.lineplot(ax=axs[i], x=times[key], y=data[key], palette='colorblind', data=all_data[all_data['gas_type']=='CO2'], hue='flow regimes', legend='full')

	if key=='P_exit':
		# sns.lineplot(ax=axs[i], x=times['P_fg_exit'], y=all_data['P_fg_exit'], palette='colorblind', data=all_data[all_data['gas_type']=='CO2'], legend='full')
		axs[i].plot(all_data[all_data['gas_type']=='CO2']['time'], all_data[all_data['gas_type']=='CO2']['P_fg_exit'], color='red', label='Phase Change at Nozzle Temp')
		axs[i].plot(all_data[all_data['gas_type']=='CO2']['time'], all_data[all_data['gas_type']=='CO2']['P_trip'], color='green', linestyle='--', label='Triple Point')
	if key=='T_exit':
		# sns.lineplot(ax=axs[i], x=times['T_fg_exit'], y=all_data['T_fg_exit'], palette='colorblind', data=all_data[all_data['gas_type']=='CO2'], legend='full')
		axs[i].plot(all_data[all_data['gas_type']=='CO2']['time'], all_data[all_data['gas_type']=='CO2']['T_fg_exit'], color='red', label='Phase Change at Nozzle Pres')
		axs[i].plot(all_data[all_data['gas_type']=='CO2']['time'], all_data[all_data['gas_type']=='CO2']['T_trip'], color='green', linestyle='--', label='Triple Point')

	axs[i].set_ylabel(ylabels[key], color='#413839', fontsize=fontsize)
	axs[i].set_xlim(left=0)
	# axs[i].set(xscale="log")
	# axs[i].set_ylim(bottom=0)

	axs[i].legend(loc='upper right', fontsize=6, framealpha=0.9)

	yfmt = ScalarFormatterForceFormat()
	yfmt.set_powerlimits((0,0))
	axs[i].yaxis.set_major_formatter(yfmt)
	axs[i].ticklabel_format(axis='y', style='sci', scilimits=(0,0), useMathText=True)
	axs[i].tick_params(axis='y', labelsize=6, pad=0)
	axs[i].yaxis.offsetText.set_fontsize(6)

	axs[i].tick_params(axis='x', labelsize=6, pad=0)
	axs[i].xaxis.label.set_size(8)
	axs[i].set(xlabel=r'Time $(sec)$')
plt.show()

fig, axs = plt.subplots(1,1, figsize=figsize, dpi=dpi, sharex=True)
fig.suptitle('Single Plenum Discharge, In-space vs. On-ground', y=0.98)
fig.canvas.set_window_title('Nozzle Performance Metrics')
axs.set_title(r'({} at $T_0$={} K, $V_{{p}}=${} cm$^3$, Nozzle $\varnothing${} mm, $\lambda$={})'.format(gas_label, T_t_init, vol*10**6, d_star*1000, expansion_ratio), fontsize=9)

key, value = list(data.items())[0]
sns.lineplot(ax=axs, x=times[key], y=data[key], palette='colorblind', data=all_data, hue='flow regimes', style='gas_type', legend='full')

axs.set_ylabel(ylabels[key], color='#413839', fontsize=fontsize)

axs.set_xlim(left=0.1)
axs.set_ylim(bottom=0)
# if key == 'visc_losses':
# 	axs.set_ylim(top=105)

axs.set(xscale="log")
# if key == 'Re_star':
	# axs.set(yscale="log")
	# axs.set_ylim(bottom=100)


class ScalarFormatterForceFormat(mpl.ticker.ScalarFormatter):
		def _set_format(self):  # Override function that finds format to use.
			self.format = "%1.1f"  # Give format here

axs.legend(loc='upper right', fontsize=6, framealpha=0.9)

yfmt = ScalarFormatterForceFormat()
yfmt.set_powerlimits((0,0))
axs.yaxis.set_major_formatter(yfmt)
axs.ticklabel_format(axis='y', style='sci', scilimits=(0,0), useMathText=True)
axs.tick_params(axis='y', labelsize=6, pad=0)
axs.yaxis.offsetText.set_fontsize(6)

axs.tick_params(axis='x', labelsize=6, pad=0)
axs.xaxis.label.set_size(8)
axs.set(xlabel=r'Time $(sec)$')

plt.show()
