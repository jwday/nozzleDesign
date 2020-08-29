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


## ==================================================================================
## ---- USER OPTIONS ----------------------------------------------------------------
## ==================================================================================

gas_types = ['CO2']	# Gas choices: R236fa, R134a, N2, CO2, H2, air
								# Gas choice will determine geometry based on desired output that was determined in the course of this project
d_upstream = 2 / 1000			# Upstream "pipe" diameter (really just the valve orifice diameter), units of m (mm / 1000)
L_upstream = 40 / 1000			# Upstream "pipe length" aka path length through valve. Just an estimate.
T_wall = 293					# Valve brass body wall tempertaure used to evaluate heat transfer conductivity
m_brass	= 60 / 1000				# Mass brass, kg
cp_brass = 380					# Specific heat capacity of brass, J/kg-K
cutoff_cond = 0.0001			# Cutoff condition, defined by the fractional change in pressure (relative to P_t_init) per second, units of 1/sec
figsize = (7.5, 4)				# Figure size (in)
dpi = 150						# Figure dpi
fudge_factor = 1

# Choose state transition process. 'mass-energy-balance', 'isentropic', 'isenthalpic', 'isothermal'
process = 'mass-energy-balance'

# Include thermal model?
thermal_model = False




# --------------------------------------------------------------------------------
# Define a plot formatter and standardize plot formatting schemes
class ScalarFormatterForceFormat(mpl.ticker.ScalarFormatter):
		def _set_format(self):  # Override function that finds format to use.
			self.format = "%1.1f"  # Give format here
yfmt = ScalarFormatterForceFormat()
yfmt.set_powerlimits((0,0))
# sns.set()
sns.axes_style("white")
sns.set_style("whitegrid", {"xtick.major.size": 0, "ytick.major.size": 0, 'grid.linestyle': '--'})
sns.set_context("paper", font_scale = 1, rc={"grid.linewidth": .5})
# sns.set_palette("colorblind")
# default_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']



def stream_props(area_ratio, P_t, T_t):
	mach_no = supersonic_mach_anywhere(area_ratio)
	pres = P_t * (1 + L*mach_no**2)**-W
	temp = T_t * (1 + L*mach_no**2)**-1
	return mach_no, pres, temp



## ==================================================================================
## ---- BEGIN DISCHARGE CODE --------------------------------------------------------
## ==================================================================================

# Set up Pandas DataFrame to store sim parameters
all_parameters = pd.DataFrame(columns=[	'gas_type',
										'gas_label',
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

		fg_pres_from_temp, fg_temp_from_pres, phase_data = create_phase_funcs(fluid_props, P_trip, T_trip)
		# visc_func = create_visc_func_gas(fluid_props)
		
		def visc_func(P,T):
			f = 4.97025E-2*T + 1.337E-3
			return np.array([f])

		cp_func = create_cp_func(fluid_props)
		cv_func = create_cv_func(fluid_props)
		ktc_func = create_ktc_func(fluid_props)
		h_from_PT_gas_func = create_h_from_PT_gas_func(fluid_props)
		u_from_PT_gas_func = create_u_from_PT_gas_func(fluid_props)
		r_from_PT_gas_func = create_r_from_PT_gas_func(fluid_props)

		P_from_ru_func, T_from_ru_func = create_PT_from_ru_gas_func(fluid_props_vol)
		P_from_rh_func, T_from_rh_func = create_PT_from_rh_gas_func(fluid_props_vol)
		P_from_rT_func = create_P_from_rT_gas_func(fluid_props_vol)
		


	elif gas_type == 'R134a':
		gas_label = 'R134a'
		P_t_init = 82.9 * 6894.76  		# Init Total Pressure, units of Pa (psia * 6894.76)
		P_amb = 0 * 6894.76  			# Ambient Pressure, units of Pa (psia * 6894.76)
		T_t_init = 20 + 273.15  		# Init Total Temperature, units of K (C + 273.15)
		vol = 11.2 / 10**6  			# Plenum volume, units of m^3 (cm^3 / 10^6)
		d_star = 0.212 / 1000  			# Nozzle throat diameter, units of m (mm / 1000)
		half_angle = 10  				# (Conical) Nozzle expansion angle (degrees)
		expansion_ratio = 4			# Nozzle expansion ratio (Exit Area / Throat Area)
		right_limit = 28*((vol * 10**6)/11.2)*(d_star * 1000)/(0.212) # For X-axis time scale
		visc_loss_param = 4.0			# As determined from NASA TN D-3056 (This varies a lot over Re numbers)
		k = 1.127  						
		R = 8.314/0.10203  				# Specific gas constant (J/kg-K)
		T_trip = 169.85  				# Triple point temperature (K)
		P_trip = 389.56  				# Triple point pressure (Pa)
		Z = 1							# Compressibility factor (unknown, so for now = 1)

		fluid_props = pd.read_excel('R134a_props_NIST.xlsx', sheet_name=None)

		fg_pres_from_temp, fg_temp_from_pres, phase_data = create_phase_funcs(fluid_props, P_trip, T_trip)
		visc_func = create_visc_func_gas(fluid_props)
		cp_func = create_cp_func(fluid_props)
		ktc_func = create_ktc_func(fluid_props)


	# --------------------------------------------------------------------------------
	# Calculate initial properties
	# rho_t_init = P_t_init/(Z*R*(T_t_init))  								# Initial density, kg/m^3 (or g/l), using GENERALIZED COMPRESSABILITY
	rho_t_init = r_from_PT_gas_func(P_t_init, T_t_init)[0]*1000				# Initial density, kg/m^3 (or g/l), using REAL DATA
	mu_t_init = 1/rho_t_init  												# Initial specific volume, m^3/kg (or l/g)
	m_init = rho_t_init*vol  												# Initial propellant mass, kg
	u_sp_init = u_from_PT_gas_func(P_t_init, T_t_init)[0]*1000				# Initial specific internal energy, J/kg
	h_sp_init = h_from_PT_gas_func(P_t_init, T_t_init)[0]*1000				# Initial specific enthalpy, J/kg
	dia = 2*(vol*(3/4)/np.pi)**(1/3)  										# Plenum diameter, m
	time_step_init = 15.28*vol*(P_t_init-P_amb)/((P_t_init*d_star)**2)		# Time step, s

	# Recalculate P_init based on rho_t_init and h_sp_init. I know it shouldn't be different, but it is, based on whichever function is being used (h_from_PT vs P_from_rh)
	# P_t_init = P_from_rh_func(rho_t_init, h_sp_init/1000)[0]*1000000

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

	T_wall 			= [T_wall]
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
	## ---- SOLVE MACH-AREA RELATION ----------------------------------------------------
	## ==================================================================================
	# Init control items
	i = 0
	n_max = 0
	end_loop_flag = False
	time_step = time_step_init
	delta_pres = 1

	# --------------------------------------------------------------------------------
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
	
	def Z_func(M):
		Z = 1 + L*M**2
		return Z

	def rayleigh_machTempRelation(X, *S):
		M1 = S[0]				# Mach number at station 1
		T1 = S[1]				# Static temperature at station 1
		T2 = S[2]				# Static temperature at station 2
		f = ( ((1 + k*M1**2)/(1 + k*X))*(np.sqrt(X)/M1) )**2 - (T2/T1)
		return f

	def nusseltGnielinski(Re, Pr, f):
		Nu = (f/8)*(Re - 1000)*Pr / (1 + 12.7*np.sqrt(f/8)*(Pr**(2/3) - 1))
		return Nu

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
	# I'm going to define a series of Mach vs. Area Ratio, then use some linear interpolation to determine a Mach number at a given area ratio (for viscous loss integral). Should only have to do this one.
	# Currently not being used
	range_of_subsonic_mach_nos = np.linspace(0.1, 1.1, 501)
	range_of_supersonic_mach_nos = np.linspace(1.1, 5.1, 501)

	range_of_subsonic_area_ratios = [np.sqrt( ((1 + L*x**2)**Q) / ((x**2)*(P**Q)) ) for x in range_of_subsonic_mach_nos]
	range_of_supersonic_area_ratios = [np.sqrt( ((1 + L*x**2)**Q) / ((x**2)*(P**Q)) ) for x in range_of_supersonic_mach_nos]

	subsonic_mach_anywhere = interp1d(range_of_subsonic_area_ratios, range_of_subsonic_mach_nos)
	supersonic_mach_anywhere = interp1d(range_of_supersonic_area_ratios, range_of_supersonic_mach_nos)




	## ==================================================================================
	## ---- LOOP THROUGH SUPPLY PRESSURES -----------------------------------------------
	## ==================================================================================
	# Start with a P_t_init and T_t_init and repeatedly run the nozzle code
	# Nozzle code will return throat and exit flow properties (along with thrust, ISP, impulse) based on supply properties
	# Results are used to determine mass flow rate, which is then used to calculate new supply pressure
	while delta_pres > cutoff_cond and P_t_plenum[-1] > P_amb:
		if thermal_model:
			# Put the CONSTANT PRESSURE temperature change here
			# First calculate Nusselt number using the Gnielinski correlation
			# ...to do that you need a Reynolds number, Prandtl number, and Darcy friction factor
			# ...to get the Reynolds number you need mass flow rate (nozzle function) and viscosity
			# ...viscosity and mass flow rate both depend on temperature
			# Realistically, the temperature will depend on the HTC and temperature gradient during this time period

			# Step 0: Solve machAreaRelation at nozzle INLET area to get SUBSONIC mach number, Ma
			# Step 1: Run nozzle() at isentropic T_t, P_t, rho_t to return m_dot
			# Step 2: Use isentropic static temperature to get viscosity, mu
			# Step 3: Use isentropic m_dot and mu to calculate Re
			# Step 4: Use Re, Pr, and f to calculate Nu using the Gnielinski correlation
			# Step 5: Use Nu to calculate heat transfer coefficient, htc
			# Step 6: Use h, A, delta_T (T_wall - T_static), and m_dot to get Q (total specific heat, J/g) from Q_dot (J/s)
			# Step 7: Use Q = Cp*delT to estimate the delT of the fluid.
			# Step 8: Add delT to T_static, then calculate new T_t

			# Step 1:
			if i == 0:
				P_star, T_star, rho_star, Re_star, M_star, v_star, P_exit, T_exit, rho_exit, M_exit, v_exit, c_exit, m_dot, F, CF, flow_regime, area_ratio_at_shock, F_mdotv, F_pdiff = nozzle(k, R, M_crit_sub, M_crit_sup, P_t_init, T_t_init, rho_t_init, P_amb, d_star, expansion_ratio, half_angle, gas_type, visc_func, r_from_PT_gas_func)
				list_of_mdots.append(m_dot*1000)

			# Step 2:
			T_upstream = T_t_plenum[-1]/Z_func(M_sub_upstream)
			P_upstream = P_t_plenum[-1]/(Z_func(M_sub_upstream)**W)
			visc_upstream.append(visc_func(P_upstream, T_upstream)[0] / 1000000)									# Convert from uPa-s to Pa-s (equivalent to kg/m-s), and pull value out of resultant array
			
			# Step 3:
			Re_upstream.append((list_of_mdots[-1]/1000)*d_upstream/(A_upstream*visc_upstream[-1]))
			
			# Step 4:
			f = 0.07																						# Darcy friction factor. Just an estimate from Moody chart based on brass or steel or something similar.
			# f_alt = (0.79*np.log(Re_upstream[-1])-1.64)**-2												# This stays very low...like, between 0.023 and 0.035

			cp_upstream = cp_func(P_upstream, T_upstream)[0] * 1000											# Convert from J/g-K to J/kg-K
			ktc_upstream = ktc_func(P_upstream,  np.average([T_wall[-1], T_upstream]))[0]					# Thermal conductivity, W/m-K. This value is typically evaluated at the average of the wall temperature and fluid temperature
			Pr_upstream.append(visc_upstream[-1] * cp_upstream / ktc_upstream)								# Calculate Prandtl number, though this stays almost entirely between 0.784 and 0.77, so maybe just used a fixed value...say, 0.78?

			if Re_upstream[-1] >= 3000:
				Nu_upstream.append( nusseltGnielinski(Re_upstream[-1], Pr_upstream[-1], f) )
				Nu_last = Nu_upstream[-1]
				m_dot_last = m_dot
			else:
				Nu_upstream.append(m_dot*(Nu_last/m_dot_last))
				# Nu_upstream.append(4.36)

			# Step 5:
			htc_upstream = Nu_upstream[-1] * ktc_upstream / L_upstream										# Heat transfer coefficient, W/m^2-K
			
			# Step 6:
			Q_dot_upstream = htc_upstream * (np.pi*d_upstream*L_upstream) * (T_wall[-1] - T_upstream)		# Heat rate, W (J/s)
			delh_upstream = Q_dot_upstream / (list_of_mdots[-1]/1000)										# Change in specific enthalpy, J/kg
			
			# Step 7:
			delT_through_valve = delh_upstream / cp_upstream
			delT_wall = Q_dot_upstream*time_step / (m_brass*cp_brass)
			
			# Step 8:
			T_inlet.append(T_upstream + delT_through_valve)
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
			



		# Now go through and calculate nozzle stuff with this new T_t
		P_star, T_star, rho_star, Re_star, M_star, v_star, P_exit, T_exit, rho_exit, M_exit, v_exit, c_exit, m_dot, F, CF, flow_regime, area_ratio_at_shock, F_mdotv, F_pdiff = nozzle(k, R, M_crit_sub, M_crit_sup, P_t_inlet[-1], T_t_inlet[-1], rho_t_inlet[-1], P_amb, d_star, expansion_ratio, half_angle, gas_type, visc_func, r_from_PT_gas_func)

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
		list_of_flow_regimes.append(flow_regime)						# Exit Flow Regimes (underexpanded, overexpanded, shock in nozzle, subsonic)
		list_of_area_ratios_at_shock.append(area_ratio_at_shock)		# Area ratio of normal shock in nozzle, if it exists (otherwise NaN)
		list_of_F_mdotv.append(F_mdotv)									# Thrust due to momentum exchange, N
		list_of_F_pdiff.append(F_pdiff)									# Thrust due to pressure differential, N

		list_of_average_thrusts.append(np.average(list_of_thrusts))			# Time-average cumulative thrust, N
		ISP.append( 1000*list_of_thrusts[-1]/(9.81*list_of_mdots[i]) )	# Specific Impulse, s

		if i == 0:  # If we're on the first time step...
			cumulative_impulse.append( time_step*list_of_thrusts[-1])
		else:  # If we're on any other time step...
			cumulative_impulse.append(time_step*list_of_thrusts[-1] + cumulative_impulse[-1])

		## --------------------------------------------------------------------------------
		# Now we should try to see if we can determine the actual PHASE mixture of the flow AT THE EXIT ONLY (for now)
		# First, let's try to verify that the entropy is actually constant maybe? Well no, because the equations you're using assume an isentropic process (adiabatic and internally reversable), and the entropy change for an adiabatic process is 0.
		# Maybe instead we should just try to plot the phase change lines along with the actual pressure.
		# That means writing a function that goes through every sheet in the R134a (and later CO2) data and returning the (P,T) at which the phase change takes place
		# Done
		list_of_P_fg_exit.append(fg_pres_from_temp(T_exit))
		list_of_T_fg_exit.append(fg_temp_from_pres(P_exit))

		## --------------------------------------------------------------------------------
		# Now that we've established the boundaries for phase change, maybe NOW we can start looking at what the phase mixture is.
		# AFAIK, the isentropic trend will just continue into the mixed phase region.
		# So...once you determine the P,T of the mixture, see if it falls below that line.
		# If it does, first figure out what the phase bounds are at that state.
		# Then, use the known entropy (which will be affected by temperature, so...take it from T_t AFTER the point where you apply your thermal model) and determine a QUALITY 


		## --------------------------------------------------------------------------------
		# Now let's calculate viscous losses
		# This is my temporary function which is only really valid for CO2 but I'll use it for both for the time being
		list_of_visc_losses.append(visc_loss_param/np.sqrt(Re_star*np.tan(np.deg2rad(half_angle))))
		thrust_eff.append(100*(CF-list_of_visc_losses[-1])/CF)

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
			dE_dt = m_dot*(h_sp_plenum[-1])				# J/s, what if we did this instead? Basically assume that the enthalpy is the same for the gas leaving the volume...?
			delta_E = dE_dt*time_step																

			u_sp_plenum.append( ((u_sp_plenum[-1]*m_gas[-2]) - delta_E)/m_gas[-1] )
			
			P_from_ru = P_from_ru_func(rho_t_plenum[-1], u_sp_plenum[-1]/1000)[0]*1000000
			T_from_ru = T_from_ru_func(rho_t_plenum[-1], u_sp_plenum[-1]/1000)[0]

			P_t_plenum.append(P_from_ru)
			T_t_plenum.append(T_from_ru)

			h_sp_plenum.append( h_from_PT_gas_func(P_from_ru, T_from_ru)[0]*1000)					# J/s, calculate new specific enthalpy 

			list_of_P_fg_t.append(fg_pres_from_temp(T_from_ru))		# Track the phase change pressure at given temperature
			list_of_T_fg_t.append(fg_temp_from_pres(P_from_ru))		# Track the phase change temperature at given pressure

			# I'd like a third function to determine quality based on the aforementioned parameters.

		elif process == 'isothermal':
			# Isothermal process in plenum
			process_label = 'Isothermal'
			T_t_plenum.append( T_t_init )
			P_t_plenum.append( P_from_rT_func(rho_t_plenum[-1], T_t_plenum[-1])[0]*1000000 )
			u_sp_plenum.append( u_from_PT_gas_func(P_t_plenum[-1], T_t_plenum[-1])[0]*1000 )
			list_of_h_sp.append( h_from_PT_gas_func(P_t_plenum[-1], T_t_plenum[-1])[0]*1000 )

		elif process == 'isentropic':
			# Isentropic process in plenum
			process_label = 'Isentropic'
			P_t_plenum.append( P_t_plenum[-1]*(rho_t_plenum[-1]/rho_t_plenum[-2])**k )
			T_t_plenum.append( T_t_plenum[-1]*(rho_t_plenum[-1]/rho_t_plenum[-2])**(k-1) )
			u_sp_plenum.append( u_from_PT_gas_func(P_t_plenum[-1], T_t_plenum[-1])[0]*1000 )
			list_of_h_sp.append( h_from_PT_gas_func(P_t_plenum[-1], T_t_plenum[-1])[0]*1000 )

		elif process == 'isenthalpic':
			# Isenthalpic process in plenum
			process_label = 'Isenthalpic'
			list_of_h_sp.append( h_sp_init )
			P_from_rh = P_from_rh_func(rho_t_plenum[-1], list_of_h_sp[-1]/1000)[0]*1000000
			T_from_rh = T_from_rh_func(rho_t_plenum[-1], list_of_h_sp[-1]/1000)[0]
			P_t_plenum.append(P_from_rh)
			T_t_plenum.append(T_from_rh)
			u_sp_plenum.append( u_from_PT_gas_func(P_t_plenum[-1], T_t_plenum[-1])[0]*1000 )
		
		else:
			print('No state transition process defined!')
			break


		time.append((i+1)*time_step)  # The first iteration is at t=0, so the first time[] entry will be 0.


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
										'gas_label',
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
									cumulative_impulse,
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
	current_data['P_trip'] = P_trip
	current_data['T_trip'] = T_trip

	all_data = all_data.append(current_data, ignore_index=True)


	# --------------------------------------------------------------------------------
	# Create a plot of phase transition boundaries and overlay P-T path of propellant during discharge
	fig, ax = plt.subplots(1, 1, figsize=figsize, dpi=dpi)
	# Plot phase data
	sns.scatterplot(ax=ax, x='Temperature (K)', y='Pressure (Pa)', palette='colorblind', data=phase_data[phase_data['Dataset']=='NIST Data'], hue='Dataset', zorder=10)
	sns.lineplot(ax=ax, x='Temperature (K)', y='Pressure (Pa)', palette='colorblind', data=phase_data[phase_data['Dataset']=='Extrapolated'], hue='Phase')

	# Plot P-T path
	sim_flow = current_data[current_data['flow regimes'].isin(['underexpanded', 'weak shock outside', 'normal shock at exit', 'normal shock in nozzle', 'subsonic'])][['P_t', 'T_t']]
	ax.plot(sim_flow['T_t'], sim_flow['P_t'], label='Plenum, Total P-T Path', linestyle='--')
	ax.plot(sim_flow['T_t'][0], sim_flow['P_t'][0], 'o', fillstyle='none', label='Start')
	ax.plot(sim_flow['T_t'][-1:], sim_flow['P_t'][-1:], 'x', label='Finish')
	
	ax.set_title(r'{} Phase Diagram and Plenum P-T Path ({})'.format(gas_label, process_label))
	ax.set_xlabel(r'Temperature, $K$')
	ax.set_ylabel(r'Pressure, $Pa$')
	ax.set(yscale="log")
	if gas_type == 'CO2':
		ax.set_xlim([135, 280])
		ax.set_ylim([3000, 2000000])

	# Change tick formatting to make it look nicer
	ax.xaxis.label.set_size(8)
	ax.tick_params(axis='x', labelsize=6, pad=0)
	ax.tick_params(axis='y', labelsize=6, pad=0)
	ax.yaxis.set_major_formatter(yfmt)
	ax.yaxis.offsetText.set_fontsize(6)
	ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0), useMathText=True)

	
	legend, handles = ax.get_legend_handles_labels()		# Get the legend and hadles so we can modify them
	del legend[0], legend[-2], handles[0], handles[-2]		# Remove some of the dumb legend stuff that gets added by sns.scatterplot and sns.lineplot
	legend = legend[-1:] + legend[:-1]						# Move the last element to the beginning
	handles = handles[-1:] + handles[:-1]					# Move the last handle to the beginning
	ax.legend(legend, handles, loc='lower right')			# Make a new legend with the modified handles

	plt.show()


	# --------------------------------------------------------------------------------
	# Create a plot of u vs. rho to begin identifying phase transitions
	fig, ax = plt.subplots(1, 1, figsize=figsize, dpi=dpi)
	# Plot phase data
	u_fg_gas = [u_from_PT_gas_func(x,y)[0] for x,y in all_data[all_data['gas_type']==gas_type][['P_fg_t', 'T_fg_t']].values]
	r_fg_gas = [r_from_PT_gas_func(x,y)[0]*1000 for x,y in all_data[all_data['gas_type']==gas_type][['P_fg_t', 'T_fg_t']].values]
	u_r_phase_data = pd.DataFrame({'Internal Energy (kJ/kg)': u_fg_gas, 'Density (kg/m^3)': r_fg_gas})
	# sns.scatterplot(ax=ax, x='Temperature (K)', y='Pressure (Pa)', palette='colorblind', data=phase_data[phase_data['Dataset']=='NIST Data'], hue='Dataset', zorder=10)
	# sns.lineplot(ax=ax, x='Internal Energy (kJ/kg)', y='Density (kg/m^3)', palette='colorblind', data=u_r_phase_data)
	ax.plot(u_r_phase_data['Internal Energy (kJ/kg)'], u_r_phase_data['Density (kg/m^3)'], label='Phase Change', linestyle='-')

	# Plot h-p path
	sim_flow = current_data[current_data['flow regimes'].isin(['underexpanded', 'weak shock outside', 'normal shock at exit', 'normal shock in nozzle', 'subsonic'])][['u_sp', 'rho_t']]
	ax.plot([x/1000 for x in sim_flow['u_sp']], sim_flow['rho_t'], label='Plenum, Total h-r Path', linestyle='--')
	ax.plot([x/1000 for x in sim_flow['u_sp']][0], sim_flow['rho_t'][0], 'o', fillstyle='none', label='Start')
	ax.plot([x/1000 for x in sim_flow['u_sp']][-1:], sim_flow['rho_t'][-1:], 'x', label='Finish')
	
	ax.set_title(r'Plenum $\rho$-h Path ({})'.format(gas_label, process_label))
	ax.set_xlabel(r'Internal Energy, $kJ/kg$')
	ax.set_ylabel(r'Density, $kg/m^3$')
	# ax.set(yscale="log")
	if gas_type == 'CO2':
		ax.set_xlim([360, 430])
		ax.set_ylim([3, 17])

	# Change tick formatting to make it look nicer
	ax.xaxis.label.set_size(8)
	ax.tick_params(axis='x', labelsize=6, pad=0)
	ax.tick_params(axis='y', labelsize=6, pad=0)
	ax.yaxis.set_major_formatter(yfmt)
	ax.yaxis.offsetText.set_fontsize(6)
	ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0), useMathText=True)

	
	legend, handles = ax.get_legend_handles_labels()		# Get the legend and hadles so we can modify them
	del legend[0], legend[-2], handles[0], handles[-2]		# Remove some of the dumb legend stuff that gets added by sns.scatterplot and sns.lineplot
	legend = legend[-1:] + legend[:-1]						# Move the last element to the beginning
	handles = handles[-1:] + handles[:-1]					# Move the last handle to the beginning
	ax.legend(legend, handles, loc='lower right')			# Make a new legend with the modified handles

	plt.show()








all_parameters = all_parameters.set_index('gas_type')


## ==================================================================================
## ---- PLOT ------------------------------------------------------------------------
## ==================================================================================
# ---- Plot 2x3 [Thrust, Impulse, ISP, Re, Ma, Density] all vs. Inlet + Throat + Exit Pressure
linewidth = 2
fontsize = 10

data = 	{ 
			'P_t': 				all_data['P_t'],
			'T_t': 				all_data['T_t'],
			# 'rho_t':			all_data['rho_t'],
			# 'mu_t':				all_data['mu_t'],
			# 'h_sp':				all_data['h_sp'],
			# 'u_sp':				all_data['u_sp'],

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
			'Re_star': 		all_data['Re_star'],
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
			'thrust_eff': 		all_data['thrust_eff'],
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
			'thrust_eff':		'Thrust Efficiency',
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
			'thrust_eff': 	   r'Thrust Efficiency, $\%$',
			'avg_thrust':		'Time Average Thrust, $mN$',
			'cum_impulse': 		'Impulse, $mN-s$',
			'ISP': 				'$I_{SP}$, $s$', 		
			'ARs at shock':	   r'Area Ratio, $\lambda$',

			'P_fg_t':			'Saturated Pressure, $Pa$',
			'T_fg_t':			'Saturated Temperature, $K$',
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

for gas in gas_types:
	fig_sat, axs = plt.subplots(4,1, figsize=figsize, dpi=dpi, sharex=True)
	fig_sat.suptitle('Plenum Pressure & Temperature w/ Saturation Data', y=0.98)
	axs[0].set_title(r'({} at $T_0$={} K, $V_{{p}}=${} cm$^3$, Nozzle $\varnothing${} mm, $\lambda$={})'.format(all_parameters.loc[gas]['gas_label'], all_parameters.loc[gas]['T_t_init'], all_parameters.loc[gas]['vol']*10**6, all_parameters.loc[gas]['d_star']*1000, all_parameters.loc[gas]['expansion_ratio']), fontsize=9)
	fig_sat.canvas.set_window_title('Saturated Pressure Stuff')

	for i, key in enumerate(data):
		sns.lineplot(ax=axs[i], x='time', y=key, palette='colorblind', data=all_data[all_data['gas_type']==gas], hue='flow regimes', legend='full')

		if key=='P_t':
			# sns.lineplot(ax=axs[i], x=times['P_fg_t'], y=all_data['P_fg_t'], palette='colorblind', data=all_data[all_data['gas_type']==gas], legend='full')
			axs[i].plot(all_data[all_data['gas_type']==gas]['time'], all_data[all_data['gas_type']==gas]['P_fg_t'], color='red', label='Phase Change Pres at Plenum Temp')
			# axs[i].plot(all_data[all_data['gas_type']==gas]['time'], all_data[all_data['gas_type']==gas]['P_trip'], color='green', linestyle='--', label='Triple Point')
		if key=='T_t':
			# sns.lineplot(ax=axs[i], x=times['T_fg_t'], y=all_data['T_fg_t'], palette='colorblind', data=all_data[all_data['gas_type']==gas], legend='full')
			axs[i].plot(all_data[all_data['gas_type']==gas]['time'], all_data[all_data['gas_type']==gas]['T_fg_t'], color='red', label='Phase Change Temp at Plenum Pres')
			# axs[i].plot(all_data[all_data['gas_type']==gas]['time'], all_data[all_data['gas_type']==gas]['T_trip'], color='green', linestyle='--', label='Triple Point')

		# if key=='P_t':
		# 	axs[i].plot(all_data[all_data['gas_type']==gas]['time'], [fg_pres_from_temp(x) for x in all_data[all_data['gas_type']==gas]['T_t'].values], color='red', label='Phase Change Pres at Exit Temp')

		# if key=='T_t':
		# 	axs[i].plot(all_data[all_data['gas_type']==gas]['time'], [fg_temp_from_pres(x) for x in all_data[all_data['gas_type']==gas]['P_t'].values], color='red', label='Phase Change Temp at Exit Pres')
			

		axs[i].set_ylabel(ylabels[key], color='#413839', fontsize=fontsize)
		axs[i].set_xlim(left=0)
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


# --------------------------------------------------------------------------------
# Compare the two gas types together on one plot and examine the flow regimes
fig, axs = plt.subplots(1,1, figsize=figsize, dpi=dpi, sharex=True)
fig.suptitle('Single Plenum Discharge, In-space vs. On-ground', y=0.98)
fig.canvas.set_window_title('Nozzle Performance Metrics')
axs.set_title(r'({} at $T_0$={} K, $V_{{p}}=${} cm$^3$, Nozzle $\varnothing${} mm, $\lambda$={})'.format(gas_label, T_t_init, vol*10**6, d_star*1000, expansion_ratio), fontsize=9)

key, value = list(data.items())[0]
sns.lineplot(ax=axs, x=times[key], y=data[key], palette='colorblind', data=all_data, hue='flow regimes', style='gas_type', legend='full')

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


axs.legend(loc='upper right', fontsize=6, framealpha=0.9)


axs.yaxis.set_major_formatter(yfmt)
axs.ticklabel_format(axis='y', style='sci', scilimits=(0,0), useMathText=True)
axs.tick_params(axis='y', labelsize=6, pad=0)
axs.yaxis.offsetText.set_fontsize(6)

axs.tick_params(axis='x', labelsize=6, pad=0)
axs.xaxis.label.set_size(8)
axs.set(xlabel=r'Time $(sec)$')

plt.show()
