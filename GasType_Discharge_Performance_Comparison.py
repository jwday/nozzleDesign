# Nozzle master control
# This code is only for plotting nozzle and impulse for Chapter 2 of thesis. Nothing else.
from nozzle_code_3 import nozzle
import matplotlib.pyplot as plt
from scipy import stats
import pandas as pd
import numpy as np
import seaborn as sns
import seaborn as sns
import matplotlib as mpl
import scipy.optimize as opti
import itertools
from nozzle_helperFuncs import *
import matplotlib.ticker as ticker


## ==================================================================================
## ---- USER OPTIONS ----------------------------------------------------------------
## ==================================================================================

gas_types = ['CO2', 'R134a']				# Gas choices: R236fa, R134a, N2, CO2, H2, air
								# Gas choice will determine geometry based on desired output that was determined in the course of this project
d_upstream = 2 / 1000			# Upstream "pipe" diameter (really just the valve orifice diameter), units of m (mm / 1000)
L_upstream = 40 / 1000			# Upstream "pipe length" aka path length through valve. Just an estimate.
T_wall_init = 293					# Valve brass body wall tempertaure used to evaluate heat transfer conductivity
m_brass	= 70 / 1000				# Mass brass, kg
cp_brass = 380					# Specific heat capacity of brass, J/kg-K

figsize = (6, 5)				# Figure size (in)
dpi = 150						# Figure dpi
fudge_factor = 1
half_angle_conv = 110/2

# Choose state transition process. 'mass-energy-balance', 'isentropic', 'isenthalpic', 'isothermal'
process = 'mass-energy-balance'
# process = 'isentropic'

# Include thermal model?
thermal_model = True




# --------------------------------------------------------------------------------
# Define a plot formatter and standardize plot formatting schemes
class ScalarFormatterForceFormat(mpl.ticker.ScalarFormatter):
		def _set_format(self):  # Override function that finds format to use.
			self.format = "%1.1f"  # Give format here
sns.axes_style("white")
sns.set_style("whitegrid", {"xtick.major.size": 0, "ytick.major.size": 0, 'grid.linestyle': '--'})
sns.set_context("paper", font_scale = 1, rc={"grid.linewidth": .5})
# sns.set_palette("colorblind")
# default_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']




## ==================================================================================
## ---- BEGIN DISCHARGE CODE --------------------------------------------------------
## ==================================================================================

# Set up Pandas DataFrame to store sim parameters
all_parameters = pd.DataFrame(columns=[	'gas_type',
										'Propellant',
										'P_t_init',
										'P_amb',
										'T_t_init',
										'vol',
										'd_star',
										'half_angle',
										'expansion_ratio'])

# Set up Pandas DataFrame to store aggregate data
all_data = pd.DataFrame(columns=[	'time',
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
sat_data = pd.DataFrame(['gas type', 'real temp', 'real pres', 'interpolated temp', 'interpolated pres'])


# --------------------------------------------------------------------------------
for gas_type in gas_types:
	# Select gas properties
	if gas_type == 'CO2':
		gas_label = 'CO$_{2}$'
		P_t_init = 114.7 * 6894.76  	# Init Total Pressure, units of Pa (psia * 6894.76)
		P_amb = 14.7 * 6894.76  		# Ambient Pressure, units of Pa (psia * 6894.76)
		cutoff_cond = 0.0001			# Cutoff condition, defined by the fractional change in pressure (relative to P_t_init) per second, units of 1/sec
		T_t_init = 0 + 273.15  			# Init Total Temperature, units of K (C + 273.15)
		vol = 30 / 10**6  				# Plenum volume, units of m^3 (cm^3 / 10^6)
		d_star = 0.6 / 1000  			# Nozzle throat diameter, units of m (mm / 1000)
		half_angle = 10  				# (Conical) Nozzle expansion angle (degrees)
		expansion_ratio = 1.17			# Nozzle expansion ratio (Exit Area / Throat Area)
		right_limit = 1.4*((vol * 10**6)/30)*(d_star * 1000)/(0.6) # For X-axis time scale
		visc_loss_param = 0.39			# As determined from NASA TN D-3056 (This stays pretty constant)
		k = 1.289						# Specific heat ratio (NOT thermal conductivity)
		R = 8.314/0.04401  				# Specific gas constant (J/kg-K)
		L_lv = 574						# Enthalpy (latent heat) of vaporiation, kJ/kg
		L_sl = 184						# Enthalpy (latent heat) of fusion, kJ/kg
		L_sv = 571						# Enthalpy (latent heat) of sublimation, kJ/kg
		T_trip = 216.58  				# Triple point temperature (K)
		P_trip = 518500  				# Triple point pressure (Pa)
		T_crit = 304.18					# Critical Temperature (K)
		P_crit = 7382500				# Critical Pressure (Pa)
		Z = 0.94						# Compressibility factor (doesn't vary much at such low pressures)
		
		fluid_props = pd.read_excel('CO2_props_NIST.xlsx', sheet_name=None)
		fluid_props_vol = pd.read_excel('CO2_props_costVol_NIST.xlsx', sheet_name=None)

		# def visc_func(P,T):
		# 	f = 4.97025E-2*T + 1.337E-3
		# 	return np.array([f])


	elif gas_type == 'R134a':
		gas_label = 'R134a'
		P_t_init = 60 * 6894.76  		# Init Total Pressure, units of Pa (psia * 6894.76)
		P_amb = 0 * 6894.76  			# Ambient Pressure, units of Pa (psia * 6894.76)
		cutoff_cond = 0.005				# Cutoff condition, defined by the fractional change in pressure (relative to P_t_init) per second, units of 1/sec
		T_t_init = 20 + 273.15  		# Init Total Temperature, units of K (C + 273.15)
		vol = 11.2 / 10**6  			# Plenum volume, units of m^3 (cm^3 / 10^6)
		d_star = 0.212 / 1000  			# Nozzle throat diameter, units of m (mm / 1000)
		half_angle = 10  				# (Conical) Nozzle expansion angle (degrees)
		expansion_ratio = 30			# Nozzle expansion ratio (Exit Area / Throat Area)
		right_limit = 28*((vol * 10**6)/11.2)*(d_star * 1000)/(0.212) # For X-axis time scale
		visc_loss_param = 4.0			# As determined from NASA TN D-3056 (This varies a lot over Re numbers)
		k = 1.127  						
		R = 8.314/0.10203  				# Specific gas constant (J/kg-K)
		T_trip = 169.85  				# Triple point temperature (K)
		P_trip = 389.56  				# Triple point pressure (Pa)
		Z = 1							# Compressibility factor (unknown, so for now = 1)

		fluid_props = pd.read_excel('R134a_props_NIST.xlsx', sheet_name=None)
		fluid_props_vol = pd.read_excel('R134a_props_const_vol_NIST.xlsx', sheet_name=None)


	fg_pres_from_temp, fg_temp_from_pres, phase_data = create_phase_funcs(fluid_props, P_trip, T_trip)
	visc_func = create_visc_func_gas(fluid_props)
	

	cp_func = create_cp_func(fluid_props)
	cv_func = create_cv_func(fluid_props)
	ktc_func = create_ktc_func(fluid_props)
	h_from_PT_gas_func = create_h_from_PT_gas_func(fluid_props)		# Returns kJ/kg for (Pa, K) input
	u_from_PT_gas_func = create_u_from_PT_gas_func(fluid_props)		# Returns kJ/kg for (Pa, K) input
	r_from_PT_gas_func = create_r_from_PT_gas_func(fluid_props)		# Returns g/ml for (Pa, K) input

	# Add r u data to "phase_data" (there's gotta be a better way)
	phase_data['Internal Energy, v (kJ/kg)'] = [u_from_PT_gas_func(phase_data['Pressure (Pa)'].iloc[x], phase_data['Temperature (K)'].iloc[x])[0] for x in phase_data.index]
	phase_data['Density, v (kg/m^3)'] = [r_from_PT_gas_func(phase_data['Pressure (Pa)'].iloc[x], phase_data['Temperature (K)'].iloc[x])[0]*1E3 for x in phase_data.index]
	phase_data['Enthalpy, v (kJ/kg)'] = [h_from_PT_gas_func(phase_data['Pressure (Pa)'].iloc[x], phase_data['Temperature (K)'].iloc[x])[0] for x in phase_data.index]

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
	# Calculate full nozzle geometry, from converging section to exit
	length_inlet = 0																				# Inlet length prior to converging section, m (mm / 1000)
	length_conv = (((d_upstream/2) - d_star)/2) * np.tan(np.radians(90 - half_angle_conv))							# Converging section length, m
	length_throat = 0.0 / 1000																				# Throat length (estimated), m
	length_nozzle = ((np.sqrt(expansion_ratio*d_star**2) - d_star)/2) / np.tan(np.radians(half_angle))		# Nozzle length, m
	total_length = length_inlet + length_conv + length_throat + length_nozzle

	def area_ratio_at_pos(pos):
		if pos <= length_inlet:
			area_ratio = ((d_upstream/2)/d_star)**2

		elif (pos > length_inlet) and (pos <= (length_inlet + length_conv)):
			dist_from_throat = length_inlet + length_conv - pos
			y = dist_from_throat / np.tan(np.radians(90 - half_angle_conv))
			area_ratio = ((2*y + d_star)/d_star)**2

		elif pos > (length_inlet + length_conv) and pos <= (length_inlet + length_conv + length_throat):
			area_ratio = 1

		elif pos > (length_inlet + length_conv + length_throat) and pos <= total_length:
			dist_from_throat = pos - (length_inlet + length_conv + length_throat)
			y = dist_from_throat * np.tan(np.radians(half_angle))
			area_ratio = ((2*y + d_star)/d_star)**2

		else:
			area_ratio = NaN

		return area_ratio
	
	geometry_length = list(np.linspace(0, total_length, 5000))
	geometry_ARs = [area_ratio_at_pos(x) for x in geometry_length]
	geometry_rads = [d_star*np.sqrt(x)/2 for x in geometry_ARs]

	# geo_length_plot = geometry_length.extend([geometry_length[-1], geometry_length[0], geometry_length[0]])
	# geo_plot_rads = geometry_rads.extend([geometry_rads[0]*1.2, geometry_rads[0]*1.2, geometry_rads[0]])

		

	## ==================================================================================
	## ---- INIT DATA LISTS -------------------------------------------------------------
	## ==================================================================================
	time = [0]
	P_t_plenum		= [P_t_init] 
	T_t_plenum		= [T_t_init]
	rho_t_plenum 	= [rho_t_init]
	mu_t_plenum 	= [mu_t_init]
	u_sp_plenum 	= [u_sp_init]
	h_sp_plenum 	= [h_sp_init]
	list_of_qual_plenum = []

	Re_upstream = []
	Nu_upstream = []
	Pr_upstream = []
	visc_upstream 	= []
	P_t_inlet		= [] 
	T_t_inlet 		= []
	T_inlet 		= []
	rho_t_inlet 	= []
	mu_t_inlet 		= []
	M_inlet			= []

	T_wall 			= [T_wall_init]
	m_gas 			= [m_init]
	list_of_mdots 	= []

	list_of_P_stars = []
	list_of_T_stars = []
	list_of_rho_stars = []
	list_of_M_stars = []
	list_of_Re_stars = []

	list_of_v_stars = []

	list_of_P_exits = []
	list_of_v_exits = []
	list_of_c_exits = []
	list_of_M_exits = []

	list_of_thrusts = []
	list_of_pressure_ratios = []
	list_of_T_exits = []
	list_of_rho_exits = []
	list_of_thrust_coeffs = []
	list_of_visc_losses = []
	thrust_eff = []
	list_of_F_mdotv = []
	list_of_F_pdiff = []

	list_of_average_thrusts = []
	cumulative_impulse = []
	ISP = []
	cumulative_mass = []

	list_of_flow_regimes = [] 
	list_of_area_ratios_at_shock = []

	list_of_P_fg_exit = []
	list_of_T_fg_exit = []
	list_of_P_fg_t = []
	list_of_T_fg_t = []

	iter = []

	## ==================================================================================
	## ---- HELPER FUNCTIONS ------------------------------------------------------------
	## ==================================================================================
	
	# The common term seen in most Mach number relations, defined here as a function for ease
	def Z_func(M):
		Z = 1 + L*M**2
		return Z

	# Mach number-Temperature relation for Rayleigh flow
	def rayleigh_machTempRelation(X, *S):
		M1 = S[0]				# Mach number at station 1
		T1 = S[1]				# Static temperature at station 1
		T2 = S[2]				# Static temperature at station 2
		f = ( ((1 + k*M1**2)/(1 + k*X))*(np.sqrt(X)/M1) )**2 - (T2/T1)
		return f

	# Gnielinski correlation for estimating the Nusselt number
	def nusseltGnielinski(Re, Pr, f, d_upstream, L_upstream, T_upstream, T_wall):
		Nu = (f/8)*(Re - 1000)*Pr*(1 + (d_upstream/L_upstream)**(2/3))*((T_upstream/T_wall)**0.45) / (1 + 12.7*np.sqrt(f/8)*((Pr**(2/3)) - 1))
		return Nu

	# Isentropic flow relation
	def stream_props(mach_no, P_t, T_t):
		# mach_no = supersonic_mach_anywhere(area_ratio)
		pres = P_t * (1 + L*mach_no**2)**-W
		temp = T_t * (1 + L*mach_no**2)**-1
		c = np.sqrt(k*R*temp)
		return pres, temp, c
	
	# def subsonic_mach_from_area(X, *S):
	# Let's build a function to get subsonic mach given an area. You'd think this would be easy, but I guess it isn't. We'll use interp1d for it.



	# area = S[0]
	# m_dot = S[1]
	# P_t = S[2]
	# T_t = S[3]
	# f = m_dot - (area*P_t/np.sqrt(T_t)) * np.sqrt(k/R) * np.sqrt(X) * (1 + L*X)**(-Q/2)
	# return f



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

	# --------------------------------------------------------------------------------
	# Solve the Isentropic Mach Number-Area relation for the (critical) Mach Numbers at the exit
	# NOTE: This equation is only valid for fully subsonic OR fully supersonic flow throughout the entire nozzle
	# It is not valid for the regime which includes a normal shock within the nozzle
	# It will instead just determine the "boundaries" for the validity of the isentropic equations
	# These values will be fed into the nozzle code and it will determine what's going on inside the nozzle
	def machAreaRelation(X, S):
		f = (1 + L*X)**Q - S*X*(P**Q)  # Where 'X' is M^2 and S is AREA RATIO SQUARED (so that this same function can be used for different area ratios i.e. along the length of the nozzle)
		return f

	x0 = np.array([0.00001, 20])	# List of M^2 to solve for roots over
	sol0 = opti.fsolve(machAreaRelation, x0, args=S, maxfev=100000, full_output=False, xtol=0.000000001)	# Solve for M^2 given S as an extra argument

	M_crit_sub = np.sqrt(sol0[0])  # Subsonic critical mach no.
	M_crit_sup = np.sqrt(sol0[1])  # Supersonic critical mach no.


	# --------------------------------------------------------------------------------
	# Solve the same relation for the upstream subsonic Mach number through the valve (very approximate) so that you can get the static pressure and temperature of the upstream flow
	# I approximate the Area Ratio as the area of the valve orifice to the nozzle throat area
	A_upstream = np.pi*((d_upstream)**2)/4		# This is the CROSS-SECTIONAL AREA
	S_upstream = (d_upstream/d_star)**2
	sol1 = opti.fsolve(machAreaRelation, x0, args=S_upstream, maxfev=100000, full_output=False, xtol=0.000000001)	# Solve for M^2 given S as an extra argument

	M_sub_upstream = np.sqrt(sol1[0])  # Subsonic critical mach no.
	M_sup_upstream = np.sqrt(sol1[1])  # Supersonic critical mach no.


	# This will be used to help us estimate the viscous losses
	# I'm going to define a series of Mach vs. Area Ratio, then use some linear interpolation to determine a Mach number at a given area ratio (for viscous loss integral). Should only have to do this once.
	# Currently not being used
	range_of_subsonic_mach_nos = np.linspace(0.001, 1.0, 501)
	range_of_supersonic_mach_nos = np.linspace(1.0, 5.1, 501)

	range_of_subsonic_area_ratios = [np.sqrt( ((1 + L*x**2)**Q) / ((x**2)*(P**Q)) ) for x in range_of_subsonic_mach_nos]
	range_of_supersonic_area_ratios = [np.sqrt( ((1 + L*x**2)**Q) / ((x**2)*(P**Q)) ) for x in range_of_supersonic_mach_nos]

	subsonic_mach_anywhere = interp1d(range_of_subsonic_area_ratios, range_of_subsonic_mach_nos)
	supersonic_mach_anywhere = interp1d(range_of_supersonic_area_ratios, range_of_supersonic_mach_nos)




	## ==================================================================================
	## ---- LOOP THROUGH SUPPLY PRESSURES -----------------------------------------------
	## ==================================================================================
	# Loop control items
	i = 0
	n_max = 0
	end_loop_flag = False
	time_step = time_step_init
	delta_pres = 1

	# Start with a P_t_init and T_t_init and repeatedly run the nozzle code
	# Nozzle code will return throat and exit flow properties (along with thrust, ISP, impulse) based on supply properties
	# Results are used to determine mass flow rate, which is then used to calculate new supply pressure
	while delta_pres > cutoff_cond and P_t_plenum[-1] > P_amb:
		# --------------------------------------------------------------------------------
		if thermal_model:
			# Step 1: Get initial mass flow rate (This is for the FIRST time step only!)
			if i == 0:
				P_star, T_star, rho_star, Re_star, M_star, v_star, P_exit, T_exit, rho_exit, M_exit, v_exit, c_exit, m_dot, F, CF, flow_regime, area_ratio_at_shock, F_mdotv, F_pdiff = nozzle(k, R, M_crit_sub, M_crit_sup, P_t_init, T_t_init, rho_t_init, P_amb, d_star, expansion_ratio, half_angle, gas_type, visc_func, r_from_PT_gas_func)
				list_of_mdots.append(m_dot*1000)

			# Step 2: Calculate the upstream (i.e. "in") flow properties (temperature, pressure, viscosity)
			T_upstream = T_t_plenum[-1]/Z_func(M_sub_upstream)
			P_upstream = P_t_plenum[-1]/(Z_func(M_sub_upstream)**W)
			visc_upstream.append(visc_func(P_upstream, T_upstream)[0] / 1000000)							# Convert from uPa-s to Pa-s (equivalent to kg/m-s), and pull value out of resultant array
			
			# Step 3: Calculate the upstream Reynolds number using the aforementioned flow properties
			Re_upstream.append((list_of_mdots[-1]/1000)*d_upstream/(A_upstream*visc_upstream[-1]))
			
			# Step 4: Determine the Nusselt number. There are a few different ways but this seems to be the most legitimate method.
			f = 0.07																						# Darcy friction factor. Just an estimate from Moody chart based on brass or steel or something similar.
			# f = (0.79*np.log(Re_upstream[-1])-1.64)**-2													# Another way to estimate the Darcy friction factor. This version stays very low...like, between 0.023 and 0.035

			cp_upstream = cp_func(P_upstream, T_upstream)[0] * 1000											# Specific heat at constant pressure, J/kg-K
			ktc_upstream = ktc_func(P_upstream,  np.average([T_wall[-1], T_upstream]))[0]					# Thermal conductivity of fluid, W/m-K. (This value is typically evaluated at the average of the wall temperature and fluid temperature)
			Pr_upstream.append(visc_upstream[-1] * cp_upstream / ktc_upstream)								# Prandtl number, though this stays almost entirely between 0.770 and 0.784, so maybe just used a fixed value...say, 0.777?

			if Re_upstream[-1] >= 3000:
				Nu_upstream.append( nusseltGnielinski(Re_upstream[-1], Pr_upstream[-1], f, d_upstream, L_upstream, T_upstream, T_wall[-1]) )				# Use the Gnielinski correlation to calculate the Nusselt number
			else:
				Nu_upstream.append(Nu_upstream[-1]*(list_of_mdots[-1]/list_of_mdots[-2]))					# Make the Nusselt number a linear function of mass flow rate

			# Step 5: Determine the Heat Transfer Coefficient for convective heat transfer based on the definition of the Nusselt number
			htc_upstream = fudge_factor * Nu_upstream[-1] * ktc_upstream / L_upstream						# Heat transfer coefficient, W/m^2-K
			
			# Step 6: Calculate the heat rate
			Q_dot_upstream = htc_upstream * (np.pi*d_upstream*L_upstream) * (T_wall[-1] - T_upstream)		# Heat rate, W (J/s)

			# Step 7: Calculate the change in enthalpy of the fluid due to heating
			delh_upstream = Q_dot_upstream / (list_of_mdots[-1]/1000)										# Change in specific enthalpy, J/kg
			
			# Step 8: Calculate the change in temperature (both of the fluid and of the wall) from the change in enthalpy
			delT_through_valve = delh_upstream / cp_upstream
			delT_wall = Q_dot_upstream*time_step / (m_brass*cp_brass)
			
			# Step 8.5a: Calculate alternate T_out and T_wall based on derivation you just made
			Q_dot = htc_upstream * (np.pi*d_upstream*L_upstream) * (T_wall[-1] - T_upstream)				# Heat rate, W (J/s), using the wall temperature calculated on the previous loop
			T_out = Q_dot * (1000/(Z*R*list_of_mdots[-1])) + T_upstream										# This assumes T_gas = T_upstream and is the least conservative estimate because the Delta-Temp will always be maximum

			alpha = (htc_upstream*np.pi*d_upstream*L_upstream)/(2*Z*R*list_of_mdots[-1])					# Simplication for use in the "complex" formulation which uses the average of the inlet and outlet flow temperatures. Probably a little more accurate.
			T_out_complex = (1/(1+alpha)) * (alpha*(2*T_wall[-1] - T_upstream) + T_upstream)

			T_wall_new = (-Q_dot*time_step / (m_brass*cp_brass)) + T_wall[-1]

			# Step 9: Calculate the fluid temperature coming out of the valve and pass that to the inlet of the nozzle. Use Rayleigh flow equations to calculate new fluid properties.
			# T_inlet.append(T_upstream + delT_through_valve)
			T_inlet.append(T_out)
			T_wall.append(T_wall[-1] - delT_wall)

			x1 = np.array([0.000001, 1])	# List of M to solve for roots over
			sol2 = opti.fsolve(rayleigh_machTempRelation, x1, args=(M_sub_upstream, T_upstream, T_inlet[-1]), maxfev=100000, full_output=False, xtol=0.000000001)
			M_inlet.append(np.sqrt(sol2[0]))

			P_t_inlet.append( P_t_plenum[-1]*((1+k*M_sub_upstream**2)/(1+k*M_inlet[-1]**2))*((Z_func(M_inlet[-1]))/(Z_func(M_sub_upstream)))**W )
			T_t_inlet.append(T_inlet[-1] * Z_func(M_inlet[-1]))
			rho_t_inlet.append(r_from_PT_gas_func(P_t_inlet[-1], T_t_inlet[-1])[0]*1000)

			if i == 0:
				del(list_of_mdots[0])

		else:
			P_t_inlet.append(P_t_plenum[-1])
			T_t_inlet.append(T_t_plenum[-1])
			rho_t_inlet.append(rho_t_plenum[-1])
			Re_upstream.append(0)
			Nu_upstream.append(0)
			Pr_upstream.append(0)
			visc_upstream.append(0)	
			T_inlet.append(0)
			T_wall.append(0)
			M_inlet.append(0)


		# --------------------------------------------------------------------------------
		# Now go through and calculate nozzle stuff (with or without thermally-updated T_t)
		P_star, T_star, rho_star, Re_star, M_star, v_star, P_exit, T_exit, rho_exit, M_exit, v_exit, c_exit, m_dot, F, CF, flow_regime, area_ratio_at_shock, F_mdotv, F_pdiff = nozzle(k, R, M_crit_sub, M_crit_sup, P_t_inlet[-1], T_t_inlet[-1], rho_t_inlet[-1], P_amb, d_star, expansion_ratio, half_angle, gas_type, visc_func, r_from_PT_gas_func)


		# --------------------------------------------------------------------------------
		# Append returned items to their respective lists
		list_of_P_stars.append(P_star)									# Throat Pressure, Pa
		list_of_T_stars.append(T_star)									# Throat Temperature, K
		list_of_rho_stars.append(rho_star)								# Throat Density, kg/m^s
		list_of_Re_stars.append(Re_star)								# Throat Reynolds Number			
		list_of_M_stars.append(M_star)									# Throat Mach No.
		list_of_v_stars.append(v_star)									# Throat Velocity, m/s

		list_of_P_exits.append(P_exit)									# Exit Pressure, Pa
		list_of_T_exits.append(T_exit)									# Exit Temperature, K
		list_of_rho_exits.append(rho_exit)								# Exit Density, kg/m^3
		list_of_M_exits.append(M_exit)									# Exit Mach Number
		list_of_v_exits.append(v_exit)									# Exit Velocity, m/s
		list_of_c_exits.append(c_exit)									# Exit Speed of Sound, m/s

		list_of_mdots.append(m_dot*1000)  								# Mass Flow Rate, g/s
		list_of_thrusts.append(F)										# Thrust, N
		list_of_thrust_coeffs.append(CF)								# Thrust Coefficient
		list_of_flow_regimes.append(flow_regime)						# Exit Flow Regimes (underexpanded, perfectly expanded, weak shock outside, normal shock at exit, normal shock in nozzle, subsonic, no flow)
		list_of_area_ratios_at_shock.append(area_ratio_at_shock)		# Area ratio of normal shock in nozzle, if it exists (otherwise NaN)
		list_of_F_mdotv.append(F_mdotv)									# Thrust due to momentum exchange, N
		list_of_F_pdiff.append(F_pdiff)									# Thrust due to pressure differential, N

		list_of_average_thrusts.append(np.average(list_of_thrusts))			# Time-average cumulative thrust, N
		ISP.append( 1000*list_of_thrusts[-1]/(9.81*list_of_mdots[i]) )	# Specific Impulse, s

		if i == 0:  # If we're on the first time step...
			cumulative_impulse.append( time_step*list_of_thrusts[-1])
		else:  # If we're on any other time step...
			cumulative_impulse.append(time_step*list_of_thrusts[-1] + cumulative_impulse[-1])


		# --------------------------------------------------------------------------------
		# Determine the phase coexistence boundary at the exit
		list_of_P_fg_exit.append(fg_pres_from_temp(T_exit))
		list_of_T_fg_exit.append(fg_temp_from_pres(P_exit))


		# --------------------------------------------------------------------------------
		# Estimate viscous losses (based on NASA TN-)
		list_of_visc_losses.append(visc_loss_param/np.sqrt(Re_star*np.tan(np.deg2rad(half_angle))))
		thrust_eff.append(CF-list_of_visc_losses[-1])


	

		# --------------------------------------------------------------------------------
		# Calculate these properties in preparation for the next loop. Loop will end if the newly calculated pressure drops below the ambient pressure.
		m_gas.append(m_gas[-1] - m_dot*time_step) 
		rho_t_plenum.append(m_gas[-1]/(vol))
		mu_t_plenum.append(1/rho_t_plenum[-1])

		if process == 'mass-energy-balance':
			# From mass and energy balance, and assuming the flow is isentropic (is it? Yes I believe it is. This model assumes it is NOT isentropic across the c.v., but it is isentropic afterwards, hence the "upstream" relations made here).
			process_label = 'Conservation of Energy'
			P_upstream = P_t_plenum[-1]*Z_func(M_sub_upstream)**-W									# Pa
			T_upstream = T_t_plenum[-1]*Z_func(M_sub_upstream)**-1									# K
			rho_upstream = rho_t_plenum[-1]*Z_func(M_sub_upstream)**-(W/k)

			h_upstream = h_from_PT_gas_func(P_upstream, T_upstream)[0]*1000							# J/kg
			v_upstream = M_sub_upstream*np.sqrt(k*R*T_upstream)										# m/s

			# dE_dt = m_dot*(h_upstream + (v_upstream**2)/2)										# J/s
			dE_dt = m_dot*(h_sp_plenum[-1])				# J/s, what if we did this instead? Basically assume that the enthalpy is the same for the gas leaving the volume? It makes almost no difference because the velocity component of enthalpy is extremely small (it's like Mach 0.15, so...)
			delta_E = dE_dt*time_step																

			u_sp_plenum.append( ((u_sp_plenum[-1]*m_gas[-2]) - delta_E)/m_gas[-1] )

			P_from_ru = P_from_ru_func(rho_t_plenum[-1], u_sp_plenum[-1]/1000)[0]*1000000
			P_t_plenum.append(P_from_ru)

			# It is at this point that you must check the phase of the propellent in the plenum
			plenum_phase = check_ru_phase(rho_t_plenum[-1], u_sp_plenum[-1]/1000, phase_data, P_from_ru_func, T_from_ru_func, P_trip, T_trip)
			if plenum_phase == 'vapor':
				T_from_ru = T_from_ru_func(rho_t_plenum[-1], u_sp_plenum[-1]/1000)[0]
				T_t_plenum.append(T_from_ru)
				list_of_qual_plenum.append(1)
				h_sp_plenum.append( h_from_PT_gas_func(P_t_plenum[-1], T_t_plenum[-1])[0]*1000)					# J/s, calculate new specific enthalpy 

				P_fg_t = fg_pres_from_temp(T_from_ru)															# Track the phase change pressure at given temperature
				T_fg_t = fg_temp_from_pres(P_from_ru)															# Track the phase change temperature at given pressure

			else:
				T_t_plenum.append(fg_temp_from_pres(P_t_plenum[-1]))
				h_sp_plenum.append( h_from_PT_gas_func(P_t_plenum[-1], T_t_plenum[-1])[0]*1000)

				# Now let's determine the quality
				# This will get you the Enthalpy of Sublimation at a given temperature (determined empircally using NIST data)
				enth_of_sub = y = -0.2985*T_t_plenum[-1] + 644.19
				
				# Used for determine the quality
				rho_temp = rho_t_plenum[-1]

				# Get the enthalpy for the saturated vapor phase enthalpy at this given density
				idx_of_closest_under = abs(phase_data[phase_data['Density, v (kg/m^3)'] < rho_temp]['Density, v (kg/m^3)'] - rho_temp).idxmin()
				idx_of_closest_over = abs(phase_data[phase_data['Density, v (kg/m^3)'] > rho_temp]['Density, v (kg/m^3)'] - rho_temp).idxmin()

				r_x1 = phase_data['Density, v (kg/m^3)'].iloc[idx_of_closest_under]
				r_x2 = phase_data['Density, v (kg/m^3)'].iloc[idx_of_closest_over]
				h_y1 = phase_data['Enthalpy, v (kJ/kg)'].iloc[idx_of_closest_under]
				h_y2 = phase_data['Enthalpy, v (kJ/kg)'].iloc[idx_of_closest_over]
				
				m = (h_y2-h_y1)/(r_x2-r_x1)
				h_v_sat = m*(rho_temp - r_x1) + h_y1

				list_of_qual_plenum.append(1 - ((h_v_sat - h_sp_plenum[-1]/1000) / enth_of_sub))

				P_fg_t = fg_pres_from_temp(T_t_plenum[-1])															# Track the phase change pressure at given temperature
				T_fg_t = fg_temp_from_pres(P_t_plenum[-1])															# Track the phase change temperature at given pressure



		elif process == 'isothermal':
			# Isothermal process in plenum
			process_label = 'Isothermal'
			T_t_plenum.append( T_t_init )
			P_t_plenum.append( P_from_rT_func(rho_t_plenum[-1], T_t_plenum[-1])[0]*1000000 )
			u_sp_plenum.append( u_from_PT_gas_func(P_t_plenum[-1], T_t_plenum[-1])[0]*1000 )
			h_sp_plenum.append( h_from_PT_gas_func(P_t_plenum[-1], T_t_plenum[-1])[0]*1000 )
			P_fg_t = fg_pres_from_temp(T_t_plenum[-1])															# Track the phase change pressure at given temperature
			T_fg_t = fg_temp_from_pres(P_t_plenum[-1])															# Track the phase change temperature at given pressure
			list_of_qual_plenum.append(1)

		elif process == 'isentropic':
			# Isentropic process in plenum
			process_label = 'Isentropic'
			P_t_plenum.append( P_t_plenum[-1]*(rho_t_plenum[-1]/rho_t_plenum[-2])**k )
			T_t_plenum.append( T_t_plenum[-1]*(rho_t_plenum[-1]/rho_t_plenum[-2])**(k-1) )
			u_sp_plenum.append( u_from_PT_gas_func(P_t_plenum[-1], T_t_plenum[-1])[0]*1000 )
			h_sp_plenum.append( h_from_PT_gas_func(P_t_plenum[-1], T_t_plenum[-1])[0]*1000 )
			P_fg_t = fg_pres_from_temp(T_t_plenum[-1])															# Track the phase change pressure at given temperature
			T_fg_t = fg_temp_from_pres(P_t_plenum[-1])															# Track the phase change temperature at given pressure
			list_of_qual_plenum.append(1)

		elif process == 'isenthalpic':
			# Isenthalpic process in plenum
			process_label = 'Isenthalpic'
			h_sp_plenum.append( h_sp_init )
			P_from_rh = P_from_rh_func(rho_t_plenum[-1], h_sp_plenum[-1]/1000)[0]*1000000
			T_from_rh = T_from_rh_func(rho_t_plenum[-1], h_sp_plenum[-1]/1000)[0]
			P_t_plenum.append(P_from_rh)
			T_t_plenum.append(T_from_rh)
			u_sp_plenum.append( u_from_PT_gas_func(P_t_plenum[-1], T_t_plenum[-1])[0]*1000 )
			P_fg_t = fg_pres_from_temp(T_t_plenum[-1])															# Track the phase change pressure at given temperature
			T_fg_t = fg_temp_from_pres(P_t_plenum[-1])															# Track the phase change temperature at given pressure
			list_of_qual_plenum.append(1)
		
		else:
			print('No state transition process defined!')
			break
		
		list_of_P_fg_t.append(P_fg_t)
		list_of_T_fg_t.append(T_fg_t)

		time.append((i+1)*time_step)  # The first iteration is at t=0, so the first time[] entry will be 0.

		# --------------------------------------------------------------------------------
		# Loop counter to check cutoff condition
		if i>=2:
			delta_pres = np.absolute((P_t_plenum[-2] - P_t_plenum[-1])/P_t_init)/time_step
		print('Time step: ' + str(round(time_step, 6)) + ' sec, P_t: ' + str(round(P_t_plenum[-1]/ 6894.76, 1)) + ' psia, Change in pres: ' + str(round(delta_pres*100, 3)) + '%/sec ' + str(round(time[-1], 4)) + ' sec', end='\r', flush=True)
		iter.append(i)
		i+=1
	print('\n')

	# By the nature of this loop, anything that has an init value will end up with one extra element in its list
	# So we must manually remove the last element once all is said and done in order to make all the array lengths the same
	del m_gas[-1], rho_t_plenum[-1], T_t_plenum[-1], P_t_plenum[-1], time[-1]



	## ==================================================================================
	## ---- ALL DATA SAVED TO PD DATAFRAME ----------------------------------------------
	## ==================================================================================
	single_sim_params = pd.DataFrame([[	gas_type,
										gas_label,
										P_t_init,
										P_amb,
										T_t_init,
										vol,
										d_star,
										half_angle,
										expansion_ratio]],
							columns=[	'gas_type',
										'Propellant',
										'P_t_init',
										'P_amb',
										'T_t_init',
										'vol',
										'd_star',
										'half_angle',
										'expansion_ratio'])

	all_parameters = all_parameters.append(single_sim_params)

	current_data = pd.DataFrame(zip(time,
									iter,
									P_t_plenum, 
									T_t_plenum, 
									rho_t_plenum,
									mu_t_plenum,
									h_sp_plenum,
									u_sp_plenum,
									list_of_qual_plenum,

									Re_upstream,
									Nu_upstream,
									Pr_upstream,
									visc_upstream,

									P_t_inlet, 
									T_t_inlet,
									T_inlet,
									T_wall,
									M_inlet,

									list_of_P_stars, 
									list_of_T_stars, 
									list_of_rho_stars, 
									list_of_Re_stars,
									list_of_M_stars,
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
									thrust_eff,
									list_of_average_thrusts,
									[x*1000 for x in cumulative_impulse],
									ISP,
									list_of_flow_regimes,
									list_of_area_ratios_at_shock,

									list_of_P_fg_t,
									list_of_T_fg_t,
									list_of_P_fg_exit,
									list_of_T_fg_exit),
						columns = [	'time',
									'iter',
									'P_t',
									'T_t',
									'rho_t',
									'mu_t',
									'h_sp',
									'u_sp',
									'X',

									'Re_up',
									'Nu_up',
									'Pr_up',
									'visc_up',

									'P_t_in',
									'T_t_in',
									'T_in',
									'T_wall',
									'M_in',

									'P_star',
									'T_star',
									'rho_star',
									'Re_star',
									'M_star',
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
									'visc_losses',
									'thrust_eff',
									'avg_thrust',
									'cum_impulse',
									'ISP',
									'flow regimes',
									'area ratios at shock',
									
									'P_fg_t',
									'T_fg_t',
									'P_fg_exit',
									'T_fg_exit'])
	current_data['gas_type'] = gas_type
	current_data['Propellant'] = gas_label
	current_data['P_trip'] = P_trip
	current_data['T_trip'] = T_trip

	if not current_data[current_data['flow regimes'] == 'No Flow'].empty:
		idx_drop = current_data[current_data['flow regimes'] == 'No Flow'].index[0]
		current_data.drop(idx_drop, inplace=True)

	all_data = all_data.append(current_data, ignore_index=True)



	# --------------------------------------------------------------------------------
	# Create a plot of phase transition boundaries and overlay P-T path of propellant during discharge
	fig, ax = plt.subplots(1, 1, figsize=(6,3), dpi=dpi)
	ax.set_title(r'Phase Diagram and Plenum P-T Path ({})'.format(process_label))
	# ax.set_title(r'({} at $P_0$={} kPa, $T_0$={} K, $V_{{p}}=${} cm$^3$, Nozzle $\varnothing${} mm, $\lambda$={})'.format(gas_label, 791, 273, vol*10**6, d_star*1000, expansion_ratio), fontsize=7)
	ax.set_xlabel(r'Temperature, $K$')
	ax.set_ylabel(r'Pressure, $Pa$')
	ax.set(yscale="log")

	# Plot phase data
	sns.scatterplot(ax=ax, x='Temperature (K)', y='Pressure (Pa)', palette='colorblind', data=phase_data[phase_data['Dataset']=='NIST Data'], hue='Dataset', zorder=10)
	sns.lineplot(ax=ax, x='Temperature (K)', y='Pressure (Pa)', palette='colorblind', data=phase_data[phase_data['Dataset']=='Extrapolated'], hue='Phase')

	# Plot Plenum and Nozzle Exit P-T path
	sim_flow = current_data[current_data['flow regimes'].isin(['Underexpanded', 'Overexpanded', 'Normal Shock', 'Normal Shock in Nozzle', 'Subsonic'])][['P_t', 'T_t', 'P_exit', 'T_exit']]
	ax.plot(sim_flow['T_t'], sim_flow['P_t'], label='Plenum P-T Path', linestyle='--', color='orange')
	# ax.plot(sim_flow['T_exit'], sim_flow['P_exit'], label='Exit P-T Path', linestyle='-.', color='crimson')

	ax.plot(sim_flow['T_t'][0], sim_flow['P_t'][0], 'o', fillstyle='none', label='Start',markeredgecolor='green')
	ax.plot(sim_flow['T_t'][-1:], sim_flow['P_t'][-1:], 'x', label='Finish', zorder=20, markeredgecolor='red')

	# ax.plot(sim_flow['T_exit'][0], sim_flow['P_exit'][0], 'o', fillstyle='none', markeredgecolor='green')
	# ax.plot(sim_flow['T_exit'][-1:], sim_flow['P_exit'][-1:], 'x', zorder=20, markeredgecolor='red')
	

	if gas_type == 'CO2':
		ax.set_xlim([160, 280])
		ax.set_ylim([80000, 1000000])

		ax.text(185, 4.5E5, 'Solid Phase', style='italic', fontsize=7,
			bbox={'facecolor': 'red', 'alpha': 0.2, 'pad': 5})
		ax.text(215, 1.2E5, 'Vapor Phase', style='italic', fontsize=7,
			bbox={'facecolor': 'red', 'alpha': 0.2, 'pad': 5})

	if gas_type == 'R134a':
		ax.set_xlim([150, 300])
		ax.set_ylim([100, 1000000])

	# Manually place the grid ticks
	ax.set_yticks([1E5, 2E5, 4E5, 6E5, 8E5, 1E6], minor=False)
	ax.yaxis.grid(True, which='major')

	# Change tick formatting to make it look nicer
	yfmt = ScalarFormatterForceFormat()
	yfmt.set_powerlimits((0,0))
	ax.xaxis.label.set_size(8)
	ax.tick_params(axis='x', labelsize=7, pad=0)
	ax.tick_params(axis='y', labelsize=7, pad=0)
	ax.yaxis.set_major_formatter(yfmt)
	ax.yaxis.offsetText.set_fontsize(6)
	ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0), useMathText=True)

	# Scale the y-axis so it's in more reasonable units. This only adjusts the labels, not the actual data.
	# scale_y = 1E-3
	# ticks_y = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x/scale_y))
	# ax.yaxis.set_major_formatter(ticks_y)
	
	# Manually adjust the legend labels
	legend, handles = ax.get_legend_handles_labels()								# Get the legend and hadles so we can modify them
	del legend[0], handles[0], legend[-2], handles[-2]								# Remove some of the dumb legend stuff that gets added by sns.scatterplot and sns.lineplot
	legend = legend[-1:] + legend[:-1]												# Move the last element to the beginning
	handles = handles[-1:] + handles[:-1]											# Move the last handle to the beginning
	ax.legend(legend, handles, loc='lower right', fontsize=7, framealpha=0.9)		# Make a new legend with the modified handles

	plt.tight_layout()
	plt.show()




	# --------------------------------------------------------------------------------
	# Create a plot of u vs. rho to identify point of phase transition
	fig, ax = plt.subplots(1, 1, figsize=(6,3), dpi=dpi)
	ax.set_title(r'Plenum $\rho$-u Path ({})'.format(gas_label, process_label))
	# ax.set(yscale="log")

	# Plot phase data
	sns.scatterplot(ax=ax, x='Internal Energy, v (kJ/kg)', y='Density, v (kg/m^3)', palette='colorblind', data=phase_data[phase_data['Dataset']=='NIST Data'], hue='Dataset', zorder=10)
	sns.lineplot(ax=ax, x='Internal Energy, v (kJ/kg)', y='Density, v (kg/m^3)', palette='colorblind', data=phase_data[phase_data['Dataset']=='Extrapolated'], hue='Phase')

	# Plot r-u path
	sim_flow = current_data[current_data['flow regimes'].isin(['Underexpanded', 'Overexpanded', 'Normal Shock', 'Normal Shock in Nozzle', 'Subsonic'])][['u_sp', 'rho_t']]
	ax.plot([x/1000 for x in sim_flow['u_sp']], sim_flow['rho_t'], label=r'Plenum, Total $\rho$-u Path', linestyle='--')
	ax.plot([x/1000 for x in sim_flow['u_sp']][0], sim_flow['rho_t'][0], 'o', fillstyle='none', label='Start')
	ax.plot([x/1000 for x in sim_flow['u_sp']][-1:], sim_flow['rho_t'][-1:], 'x', label='Finish')
	
	# Change tick formatting to make it look nicer
	ax.xaxis.label.set_size(8)
	ax.tick_params(axis='x', labelsize=6, pad=0)
	ax.tick_params(axis='y', labelsize=6, pad=0)
	yfmt = ScalarFormatterForceFormat()
	yfmt.set_powerlimits((0,0))
	ax.yaxis.set_major_formatter(yfmt)
	ax.yaxis.offsetText.set_fontsize(6)
	ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0), useMathText=True)
	ax.set_xlabel(r'Internal Energy, $kJ/kg$')				# Gotta do it here because sns.lineplot takes over the formatting if you do it before. I think.
	ax.set_ylabel(r'Density, $kg/m^{3}$')
	
	legend, handles = ax.get_legend_handles_labels()		# Get the legend and hadles so we can modify them

	if gas_type == 'CO2':
		ax.set_xlim([350, 430])
		ax.set_ylim([2, 17])

		handles[1] = 'Saturated Vapor-Solid Boundary'
		handles[2] = 'Saturated Liquid Boundary'
		handles[3] = 'Saturated Vapor Boundary'
		del legend[0], handles[0], legend[-2], handles[-2]		# Remove some of the dumb legend stuff that gets added by sns.scatterplot and sns.lineplot
		legend = legend[-1:] + legend[:-1]						# Move the last element to the beginning
		handles = handles[-1:] + handles[:-1]					# Move the last handle to the beginning
		ax.legend(legend, handles, loc='upper left')			# Make a new legend with the modified handles

		ax.text(408, 6.5, 'Vapor Phase', style='italic',
			bbox={'facecolor': 'red', 'alpha': 0.2, 'pad': 7})
		ax.text(354, 5.5, 'Two-Phase Solid/Vapor', style='italic',
			bbox={'facecolor': 'red', 'alpha': 0.2, 'pad': 7})

	if gas_type == 'R134a':
		ax.set_xlim([300, 400])
		ax.set_ylim([-1, 23])

		handles[1] = 'Saturated Vapor-Solid Boundary'
		handles[2] = 'Saturated Liquid Boundary'
		handles[3] = 'Saturated Vapor Boundary'
		del legend[0], handles[0], legend[-2], handles[-2]		# Remove some of the dumb legend stuff that gets added by sns.scatterplot and sns.lineplot
		legend = legend[-1:] + legend[:-1]						# Move the last element to the beginning
		handles = handles[-1:] + handles[:-1]					# Move the last handle to the beginning
		ax.legend(legend, handles, loc='upper left')			# Make a new legend with the modified handles

		ax.text(328, 8, 'Two-Phase Liquid/Vapor', style='italic',
			bbox={'facecolor': 'red', 'alpha': 0.2, 'pad': 7})
		ax.text(375, 2, 'Vapor Phase', style='italic',
			bbox={'facecolor': 'red', 'alpha': 0.2, 'pad': 7})

	plt.tight_layout()
	plt.show()








all_parameters = all_parameters.set_index('gas_type')


## ==================================================================================
## ---- PLOT ------------------------------------------------------------------------
## ==================================================================================
# ---- Plot 2x3 [Thrust, Impulse, ISP, Re, Ma, Density] all vs. Inlet + Throat + Exit Pressure
linewidth = 2
fontsize = 8

data = 	{ 
			'P_t': 				all_data['P_t'],
			# 'T_t': 				all_data['T_t'],
			# 'rho_t':			all_data['rho_t'],
			# 'mu_t':				all_data['mu_t'],
			# 'h_sp':				all_data['h_sp'],
			# 'u_sp':				all_data['u_sp'],
			'X':				all_data['X'],

			# 'Re_up':			all_data['Re_up'],
			# 'Nu_up':			all_data['Nu_up'],
			# 'Pr_up':			all_data['Pr_up'],
			# 'visc_up':			all_data['visc_up'],

			# 'P_t_in': 			all_data['P_t_in'],
			# 'T_t_in': 			all_data['T_t_in'],
			# 'T_in': 			all_data['T_in'],
			# 'T_wall':			all_data['T_wall'],
			# 'M_in':				all_data['M_in'],

			# 'P_star': 		all_data['P_star'],
			# 'T_star': 		all_data['T_star'],
			# 'rho_star': 		all_data['rho_star'],
			# 'Re_star': 		all_data['Re_star'],
			# 'M_star': 			all_data['M_star'],
			# 'v_star': 			all_data['v_star'],

			# 'P_exit': 		all_data['P_exit'],
			# 'T_exit': 		all_data['T_exit'],
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
			# 'visc_losses': 		all_data['visc_loss'],
			# 'thrust_eff': 		all_data['thrust_eff'],
			# 'avg_thrust': 	all_data['avg_thrust'],
			# 'cum_impulse': 	all_data['cum_impulse'],
			# 'ISP': 			all_data['ISP'],
			# 'ARs at shock':	all_data['area ratios at shock']

			# 'P_fg_exit':		all_data['P_fg_t'],
			# 'T_fg_exit':		all_data['T_fg_t'],
			# 'P_fg_exit':		all_data['P_fg_exit'],
			# 'T_fg_exit':		all_data['T_fg_exit'],
			  }

figname = { 
			'P_t': 				'Plenum Total Pressure',
			'T_t': 				'Plenum Total Temperature', 							
			'rho_t': 			'Plenum Total Density',
			'mu_t': 			'Plenum Total Specific Volume',
			'h_sp':				'Specific Enthalpy',
			'u_sp':				'Specific Internal Energy',
			'X':				'Plenum Vapor Quality',

			'Re_up':			'Upstream Reynolds No.',
			'Nu_up':			'Upstream Nusselt No.',
			'Pr_up':			'Upstream Prandtl No.',
			'visc_up':			'Upstream Viscosity',

			'P_t_in':			'Inlet Total Pressure',
			'T_t_in':			'Inlet Total Temperature', 
			'T_in':				'Inlet Static Temperature', 
			'T_wall':			'Wall Temperature',
			'M_in':				'Inlet Mach No.',						

			'P_star': 			'Throat Pressure',
			'T_star': 			'Throat Temperature',
			'rho_star': 		'Throat Density', 		
			'Re_star': 			'Throat Reynold\'s No.',
			'M_star': 			'Throat Mach No.',
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
			'thrust_eff':		'Effective Thrust Coefficient',
			'avg_thrust': 		'Time Average Thrust',
			'cum_impulse': 		'Net Impulse',
			'ISP': 				'$I_{SP}$',
			'ARs at shock':		'Area Ratio at Shock',

			'P_fg_t':			'Saturated Vapor Pressure at Current',
			'T_fg_t':			'Saturated Vapor Temperature at Current',
			'P_fg_exit':		'Saturated Vapor Pressure at Exit',
			'T_fg_exit':		'Saturated Vapor Temperature at Exit'
		  }

times = {
			'P_t': 				all_data['time'],
			'T_t': 				all_data['time'],
			'rho_t':			all_data['time'],
			'mu_t':				all_data['time'],
			'h_sp':				all_data['time'],
			'u_sp':				all_data['time'],
			'X':				all_data['time'],

			'Re_up':			all_data['time'],
			'Nu_up':			all_data['time'],
			'Pr_up':			all_data['time'],
			'visc_up':			all_data['time'],

			'P_t_in':			all_data['time'],
			'T_t_in':			all_data['time'],
			'T_in':				all_data['time'],
			'T_wall':			all_data['time'],
			'M_in':				all_data['time'],

			'P_star': 			all_data['time'],
			'T_star': 			all_data['time'],
			'rho_star': 		all_data['time'],
			'Re_star': 			all_data['time'],
			'M_star': 			all_data['time'],
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
			'thrust_eff': 		all_data['time'],
			'avg_thrust': 		all_data['time'],
			'cum_impulse': 		all_data['time'],
			'ISP': 				all_data['time'],
			'ARs at shock': 	all_data['time'],

			'P_fg_t':	 	all_data['time'],
			'T_fg_t':	 	all_data['time'],
			'P_fg_exit': 	all_data['time'],
			'T_fg_exit': 	all_data['time']
		}

ylabels = {
			'P_t': 				'Plenum Total Pres, $Pa$',
			'T_t': 				'Plenum Total Temp, $K$', 						
			'rho_t': 			'Plenum Total Dens, $kg/m^3$',
			'mu_t': 			'Plenum Total Sp. Vol, $m^3/kg$',
			'h_sp':				'Specific Enthalpy, $kJ/kg$',
			'u_sp':				'Specific Internal Energy, $kJ/kg$',
			'X':				'Plenum Quality',

			'Re_up':			'Upstream Reynolds No.',
			'Nu_up':			'Upstream Nusselt No.',
			'Pr_up':			'Upstream Prandtl No.',
			'visc_up':			'Upstream Viscosity, $Pa-s$',

			'P_t_in':			'Inlet Total Pres, $Pa$',
			'T_t_in':			'Inlet Total Temp, $K$', 						
			'T_in':				'Inlet Static Temp, $K$', 						
			'T_wall':			'Wall Temp, $K$',
			'M_in':				'Inlet Mach No.',

			'P_star': 			'Throat Pressure, $Pa$',
			'T_star': 			'Throat Temperature, $K$',
			'rho_star': 		'Throat Density, $kg/m^3$', 		
			'Re_star': 			'Throat Reynold\'s No.',
			'M_star': 			'Throat Mach No.', 				
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
			'thrust_coeff':	   r'Thrust Coefficient, $C_f$',
			'visc_losses': 	   r'Viscous Losses, $C_{f_v}$',
			'thrust_eff': 	   r'Effective Thrust Coefficient',
			'avg_thrust':		'Time Average Thrust, $mN$',
			'cum_impulse': 		'Impulse, $mN-s$',
			'ISP': 				'$I_{SP}$, $s$', 		
			'ARs at shock':	   r'Area Ratio, $\lambda$',

			'P_fg_t':			'Saturated Pressure, $Pa$',
			'T_fg_t':			'Saturated Temperature, $K$',
			'P_fg_exit':		'Saturated Pressure, $Pa$',
			'T_fg_exit':		'Saturated Temperature, $K$',
		   }

legend_locs = {
			'P_t': 				'upper right',
			'T_t': 				'upper right',
			'rho_t':			'upper right',
			'mu_t':				'upper right',
			'h_sp':				'upper right',
			'u_sp':				'upper right',
			'X':				'upper right',

			'Re_up':			'upper right',
			'Nu_up':			'upper right',
			'Pr_up':			'upper right',
			'visc_up':			'upper right',

			'P_t_in':			'upper right',
			'T_t_in':			'upper right',
			'T_in':				'upper right',
			'T_wall':			'upper right',
			'M_in':				'upper right',

			'P_star': 			'upper right',
			'T_star': 			'upper right',
			'rho_star': 		'upper right',
			'Re_star': 			'upper right',
			'M_star': 			'upper right',
			'v_star': 			'upper right',

			'P_exit': 			'upper right',
			'T_exit': 			'upper right',
			'rho_exit': 		'upper right',
			'M_exit': 			'upper right',
			'v_exit': 			'upper right',
			'c_exit': 			'upper right',

			'm_gas': 			'upper right',
			'mdot': 			'upper right',
			'F_mdotv': 			'upper right',
			'F_pdiff': 			'upper right',
			'thrust': 			'upper right',
			'thrust_coeff': 	'upper right',
			'visc_losses': 		'upper left',
			'thrust_eff': 		'upper right',
			'avg_thrust': 		'upper right',
			'cum_impulse': 		'lower right',
			'ISP': 				'upper right',
			'ARs at shock': 	'upper right',

			'P_fg_t':	 	'upper right',
			'T_fg_t':	 	'upper right',
			'P_fg_exit': 	'upper right',
			'T_fg_exit': 	'upper right'
		}

titles = { 	
			'CO2': 'On-Ground Single Plenum Discharge Reynolds Number',
			'R134a': 'In-Space Single Plenum Discharge Reynolds Number'
		 }



# --------------------------------------------------------------------------------
# Let's plot Re of both gasses on one figure with two subplots of differing time scales and make sure you include ax.set_titles
fig, ax = plt.subplots(2,1, figsize=(6, 4.5), dpi=300)
fig.suptitle(r'Throat Reynolds Number for CO$_2$ and R134a', y=0.95)
for i, gas in enumerate(gas_types):
	ax[i].set_title(r'{} at $T_0$={} K, $V_{{p}}=${} cm$^3$, Nozzle $\varnothing${} mm, $\lambda$={}'.format(all_parameters.loc[gas]['Propellant'], 273.15, all_parameters.loc[gas]['vol']*10**6, all_parameters.loc[gas]['d_star']*1000, all_parameters.loc[gas]['expansion_ratio']), fontsize=7)
	ax[i].plot(all_data[all_data['gas_type']==gas]['time'], all_data[all_data['gas_type']==gas]['Re_star'])
	ax[i].set_ylabel(r'Reynolds No. ({})'.format(all_parameters.loc[gas]['Propellant']), color='#413839', fontsize=fontsize)
	ax[i].set_xlim(left=0)

	yfmt = ScalarFormatterForceFormat()
	yfmt.set_powerlimits((0,0))
	ax[i].yaxis.set_major_formatter(yfmt)
	ax[i].ticklabel_format(axis='y', style='sci', scilimits=(0,0), useMathText=True)
	ax[i].tick_params(axis='y', labelsize=6, pad=0)
	ax[i].yaxis.offsetText.set_fontsize(6)

	ax[i].tick_params(axis='x', labelsize=6, pad=0)

ax[1].xaxis.label.set_size(8)
ax[1].set(xlabel=r'Time $(sec)$')

plt.tight_layout()
plt.subplots_adjust(top=0.85)
plt.show()



# --------------------------------------------------------------------------------
# Only looking at temperatures on same graph
for gas in gas_types:
	fig, ax = plt.subplots(1,1, figsize=(6,3), dpi=300)
	fig.suptitle('Valve Upstream and Downstream Flow Total Temperatures', y=0.98)
	ax.set_title(r'({} at $T_0$={} K, $V_{{p}}=${} cm$^3$, Nozzle $\varnothing${} mm, $\lambda$={})'.format(all_parameters.loc[gas]['Propellant'], 273.15, all_parameters.loc[gas]['vol']*10**6, all_parameters.loc[gas]['d_star']*1000, all_parameters.loc[gas]['expansion_ratio']), fontsize=7)

	temp_styles=['-', '--']
	temp_labels=['Upstream flow', 'Downstream flow']
	for i, key in enumerate(data):
		ax.plot(all_data[all_data['gas_type']==gas]['time'], all_data[all_data['gas_type']==gas][key], label=temp_labels[i], linestyle=temp_styles[i])
	
	ax.set_ylabel(r'Temperature, $K$', color='#413839', fontsize=fontsize)
	ax.set_xlim(left=0)
	ax.legend()

	# yfmt = ScalarFormatterForceFormat()
	# yfmt.set_powerlimits((0,0))
	# ax.yaxis.set_major_formatter(yfmt)
	# ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0), useMathText=True)
	ax.tick_params(axis='y', labelsize=6, pad=0)
	ax.yaxis.offsetText.set_fontsize(6)

	ax.tick_params(axis='x', labelsize=6, pad=0)
	ax.xaxis.label.set_size(8)
	ax.set(xlabel=r'Time $(sec)$')

	plt.tight_layout()
	plt.subplots_adjust(top=0.85)
plt.show()



# --------------------------------------------------------------------------------
# Let's see if we can plot exit pres & sat pres at exit on same plot, and also temp on another
for gas in gas_types:
	fig_sat, axs = plt.subplots(2,1, figsize=figsize, dpi=dpi, sharex=True)
	fig_sat.suptitle('{} and {}'.format(figname[list(data.keys())[0]], figname[list(data.keys())[1]]), y=0.98)
	axs[0].set_title(r'({} at $T_0$={} K, $V_{{p}}=${} cm$^3$, Nozzle $\varnothing${} mm, $\lambda$={})'.format(all_parameters.loc[gas]['Propellant'], 273.15, all_parameters.loc[gas]['vol']*10**6, all_parameters.loc[gas]['d_star']*1000, all_parameters.loc[gas]['expansion_ratio']), fontsize=7)
	fig_sat.canvas.set_window_title('Saturated Pressure Stuff')

	for i, key in enumerate(data):
		sns.lineplot(ax=axs[i], x='time', y=key, palette='colorblind', data=all_data[all_data['gas_type']==gas], hue='flow regimes', legend='full')

		# if key=='P_t':
		# 	# sns.lineplot(ax=axs[i], x=times['P_fg_t'], y=all_data['P_fg_t'], palette='colorblind', data=all_data[all_data['gas_type']==gas], legend='full')
		# 	axs[i].plot(all_data[all_data['gas_type']==gas]['time'], all_data[all_data['gas_type']==gas]['P_fg_t'], color='red', label='Phase Change Pres at Plenum Temp')
		# 	# axs[i].plot(all_data[all_data['gas_type']==gas]['time'], all_data[all_data['gas_type']==gas]['P_trip'], color='green', linestyle='--', label='Triple Point')
		# if key=='T_t':
		# 	# sns.lineplot(ax=axs[i], x=times['T_fg_t'], y=all_data['T_fg_t'], palette='colorblind', data=all_data[all_data['gas_type']==gas], legend='full')
		# 	axs[i].plot(all_data[all_data['gas_type']==gas]['time'], all_data[all_data['gas_type']==gas]['T_fg_t'], color='red', label='Phase Change Temp at Plenum Pres')
		# 	# axs[i].plot(all_data[all_data['gas_type']==gas]['time'], all_data[all_data['gas_type']==gas]['T_trip'], color='green', linestyle='--', label='Triple Point')

		# if key=='P_t':
		# 	axs[i].plot(all_data[all_data['gas_type']==gas]['time'], [fg_pres_from_temp(x) for x in all_data[all_data['gas_type']==gas]['T_t'].values], color='red', label='Phase Change Pres at Exit Temp')

		# if key=='T_t':
		# 	axs[i].plot(all_data[all_data['gas_type']==gas]['time'], [fg_temp_from_pres(x) for x in all_data[all_data['gas_type']==gas]['P_t'].values], color='red', label='Phase Change Temp at Exit Pres')
			

		axs[i].set_ylabel(ylabels[key], color='#413839', fontsize=fontsize)
		axs[i].set_xlim(left=0)
		# axs[i].legend(loc=legend_locs[key], fontsize=6, framealpha=0.9)
		axs[i].legend().texts[0].set_text('Flow Regimes')

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




# --------------------------------------------------------------------------------
# Compare the two gas types together on one plot and examine the flow regimes
fig, axs = plt.subplots(1,1, figsize=(6,3), dpi=300, sharex=True)
fig.canvas.set_window_title('Nozzle Performance Metrics')
fig.suptitle('Throat Reynolds No. for Single Plenum Discharge', y=1)
axs.set_title('In-space vs. On-ground')

key, value = list(data.items())[0]
sns.lineplot(ax=axs, x=times[key], y=data[key], palette='colorblind', data=all_data, style='Propellant', legend='full')

axs.set_ylabel(ylabels[key], color='#413839', fontsize=fontsize)

axs.set(xscale="log")
axs.set_xlim(left=0.1)

# axs.set(yscale="log")
# axs.set_ylim(bottom=0.1)

# if key == 'visc_losses':
# 	axs.set_ylim(top=105)
# if key == 'Re_star':
	# axs.set(yscale="log")
	# axs.set_ylim(bottom=100)


axs.legend(loc='bottom right', fontsize=6, framealpha=0.9)


axs.yaxis.set_major_formatter(yfmt)
axs.ticklabel_format(axis='y', style='sci', scilimits=(0,0), useMathText=True)
axs.tick_params(axis='y', labelsize=6, pad=0)
axs.yaxis.offsetText.set_fontsize(6)

axs.tick_params(axis='x', labelsize=6, pad=0)
axs.xaxis.label.set_size(8)
axs.set(xlabel=r'Time $(sec)$')

plt.show()




# --------------------------------------------------------------------------------
# Plot properties along nozzle length
fig_noz, axs = plt.subplots(2,2, figsize=(8, 4.5), dpi=dpi, sharex=True)
fig_noz.suptitle('Flow Properties over Nozzle Length', y=0.98)
fig_noz.canvas.set_window_title('Flow Properties over Nozzle Length')
# axs[0, 0].set_title(r'({} at $T_0$={} K, $V_{{p}}=${} cm$^3$, Nozzle $\varnothing${} mm, $\lambda$={})'.format(gas_label, T_t_init, vol*10**6, d_star*1000, expansion_ratio), fontsize=9)

# Iterate across length of nozzle to get spacial distribution of properties
# Init 2D arrays
mach_no_2D = []
P_2D = []
T_2D = []
sound = []
flow_regimes = []

time_idxs = len(time)
time_idxs_to_plot = [0, 780, 860]

for i, idx in enumerate(time_idxs_to_plot):
	mach_no_2D.append(list(np.zeros(len(geometry_length))))
	P_2D.append(list(np.zeros(len(geometry_length))))
	T_2D.append(list(np.zeros(len(geometry_length))))
	sound.append(list(np.zeros(len(geometry_length))))

	if list_of_flow_regimes[idx] == 'Underexpanded':
		print('underexpanded')
		for j, x in enumerate(geometry_length):
			if x <= (length_inlet + length_conv):																							# Subsonic up to beginning of throat
				mach_no_2D[i][j] = subsonic_mach_anywhere(geometry_ARs[j])[()]
			elif x > (length_inlet + length_conv) and x <= (length_inlet + length_conv + length_throat):									# Sonic at and through the throat
				mach_no_2D[i][j] = 1
			elif x > (length_inlet + length_conv + length_throat):																			# Supersonic in diverging section
				mach_no_2D[i][j] = supersonic_mach_anywhere(geometry_ARs[j])[()]
			P_2D[i][j], T_2D[i][j], sound[i][j] = stream_props(mach_no_2D[i][j], P_t_inlet[idx], T_t_inlet[idx])

	if list_of_flow_regimes[idx] == 'Normal Shock':
		print('weak shock outside')
		for j, x in enumerate(geometry_length):
			if x <= (length_inlet + length_conv):																							# Subsonic up to beginning of throat
				mach_no_2D[i][j] = subsonic_mach_anywhere(geometry_ARs[j])[()]
			elif x > (length_inlet + length_conv) and x <= (length_inlet + length_conv + length_throat):									# Sonic at and through the throat
				mach_no_2D[i][j] = 1
			elif x > (length_inlet + length_conv + length_throat):																			# Supersonic in diverging section
				mach_no_2D[i][j] = supersonic_mach_anywhere(geometry_ARs[j])[()]
			P_2D[i][j], T_2D[i][j], sound[i][j] = stream_props(mach_no_2D[i][j], P_t_inlet[idx], T_t_inlet[idx])

	if list_of_flow_regimes[idx] == 'Normal Shock in Nozzle':
		print('normal shock in nozzle')
		shock_pos = length_inlet + length_conv + length_throat + (d_star*np.sqrt(list_of_area_ratios_at_shock[idx]) - d_star) / (2*np.tan(np.radians(half_angle)))		# Position of shock relative to beginning of whole nozzle
		shock_pos_diff = [abs(x-shock_pos) for x in geometry_length]																		# Differential between discrete nozzle locations and calculated location of shock
		shock_pos_discrete_idx = shock_pos_diff.index(min(shock_pos_diff))																	# Index of mininmum, i.e. index corresponding to discrete shock location nearest to the calculated shock location
		for j, x in enumerate(geometry_length):
			if x <= (length_inlet + length_conv):																							# Subsonic up to beginning of throat
				mach_no_2D[i][j] = subsonic_mach_anywhere(geometry_ARs[j])[()]
			elif x > (length_inlet + length_conv) and x <= (length_inlet + length_conv + length_throat):									# Sonic at and through the throat
				mach_no_2D[i][j] = 1
			elif x > (length_inlet + length_conv + length_throat) and x < geometry_length[shock_pos_discrete_idx]:							# Supersonic in diverging section up to shock
				mach_no_2D[i][j] = supersonic_mach_anywhere(geometry_ARs[j])[()]
			elif x >= geometry_length[shock_pos_discrete_idx]:																				# Subsonic in diverging section after shock
				mach_no_2D[i][j] = subsonic_mach_anywhere(geometry_ARs[j])[()]
			P_2D[i][j], T_2D[i][j], sound[i][j] = stream_props(mach_no_2D[i][j], P_t_inlet[idx], T_t_inlet[idx])

	if list_of_flow_regimes[idx] == 'Subsonic':
		print('subsonic')
		area = []
		machs = list(np.linspace(0, 1, 100))
		for mach in machs:
			area.append((list_of_mdots[idx]/1000) / ( (P_t_inlet[idx]/np.sqrt(T_t_inlet[idx])) * np.sqrt(k/R) * mach * (1 + L*(mach**2))**(-Q/2) ))
		subsonic_mach_from_area = interp1d(area, machs)
		for j, x in enumerate(geometry_length):
			area = geometry_ARs[j]*A_star
			mach_no_2D[i][j] = subsonic_mach_from_area(area)
			P_2D[i][j], T_2D[i][j], sound[i][j] = stream_props(mach_no_2D[i][j], P_t_inlet[idx], T_t_inlet[idx])

i = 0
axs[0, 0].plot([x*1000 for x in geometry_length], [y*1000 for y in geometry_rads], color='green')
axs[0, 0].plot([x*1000 for x in geometry_length], [-y*1000 for y in geometry_rads], color='green')
for Ma, c, P, T in zip(mach_no_2D, sound, P_2D, T_2D):
	label = r'{} sec ({})'.format(round(time[time_idxs_to_plot[i]],2), list_of_flow_regimes[time_idxs_to_plot[i]])
	axs[1, 0].plot([x*1000 for x in geometry_length], Ma, label=label)
	# axs[1, 0].plot([x*1000 for x in geometry_length], np.multiply(Ma, c), label=label)
	axs[0, 1].plot([x*1000 for x in geometry_length], [x/1000 for x in P], label=label)
	axs[1, 1].plot([x*1000 for x in geometry_length], T, label=label)
	i += 1

axs[0, 0].set_ylabel(r'Nozzle Radius, $mm$', color='#413839', fontsize=fontsize)
axs[1, 0].set_ylabel('Mach No.', color='#413839', fontsize=fontsize)
# axs[1, 0].set_ylabel(r'Velocity, $m/s$', color='#413839', fontsize=fontsize)
axs[0, 1].set_ylabel(r'Pressure, $kPa$', color='#413839', fontsize=fontsize)
axs[1, 1].set_ylabel(r'Temperature, $K$', color='#413839', fontsize=fontsize)

# axs[0, 0].legend(loc='upper left', fontsize=6, framealpha=0.9)
axs[1, 0].legend(loc='upper left', fontsize=6, framealpha=0.9)
axs[0, 1].legend(loc='upper left', fontsize=6, framealpha=0.9)
axs[1, 1].legend(loc='upper left', fontsize=6, framealpha=0.9)

yfmt = ScalarFormatterForceFormat()
yfmt.set_powerlimits((0,0))

axs[0, 0].yaxis.set_major_formatter(yfmt)
axs[1, 0].yaxis.set_major_formatter(yfmt)
axs[0, 1].yaxis.set_major_formatter(yfmt)
axs[1, 1].yaxis.set_major_formatter(yfmt)

axs[0, 0].ticklabel_format(axis='y', style='sci', scilimits=(0,0), useMathText=True)
axs[1, 0].ticklabel_format(axis='y', style='sci', scilimits=(0,0), useMathText=True)
axs[0, 1].ticklabel_format(axis='y', style='sci', scilimits=(0,0), useMathText=True)
axs[1, 1].ticklabel_format(axis='y', style='sci', scilimits=(0,0), useMathText=True)

axs[0, 0].tick_params(axis='y', labelsize=6, pad=0)
axs[1, 0].tick_params(axis='y', labelsize=6, pad=0)
axs[0, 1].tick_params(axis='y', labelsize=6, pad=0)
axs[1, 1].tick_params(axis='y', labelsize=6, pad=0)

axs[0, 0].yaxis.offsetText.set_fontsize(6)
axs[1, 0].yaxis.offsetText.set_fontsize(6)
axs[0, 1].yaxis.offsetText.set_fontsize(6)
axs[1, 1].yaxis.offsetText.set_fontsize(6)

axs[0, 0].tick_params(axis='x', labelsize=6, pad=0)
axs[1, 0].tick_params(axis='x', labelsize=6, pad=0)
axs[0, 1].tick_params(axis='x', labelsize=6, pad=0)
axs[1, 1].tick_params(axis='x', labelsize=6, pad=0)

axs[0, 0].xaxis.label.set_size(8)
axs[1, 0].xaxis.label.set_size(8)
axs[0, 1].xaxis.label.set_size(8)
axs[1, 1].xaxis.label.set_size(8)

axs[1, 0].set(xlabel=r'Nozzle Length, $mm$')
axs[1, 1].set(xlabel=r'Nozzle Length, $mm$')

axs[0, 0].set_ylim(bottom=-0.55, top=0.55)

plt.show()



# We're going to try estimating viscous losses now
# mu_w = visc_func(list_of_P_stars[-1], T_wall[-1])				# Pseudo-wall viscosity (using throat Pressure and valve wall temperature)
# mu_star = visc_func(list_of_P_stars[-1], list_of_T_stars[-1])	# Throat viscosity

# f0 = gamma * np.sqrt( (mu_w/m_star) * (list_of_rho_stars[-1]/rho_t_inlet[-1]) * np.sqrt(list_of_T_stars[-1]/T_t_inlet[-1]) ) / np.sqrt( 2*(T_wall[-1]/T_t_inlet[-1]) )

# for pos in np.range(0, length_nozzle, 100):
# 	AR = area_ratio_at_pos(pos)
# 	mach_no = supersonic_mach_anywhere(AR)
# 	temp = T_t * (1 + L*mach_no**2)**-1
# 	f1 = 1/np.sqrt(np.sqrt(AR)-1)