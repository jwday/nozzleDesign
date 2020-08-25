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
import scipy.optimize as opti
import itertools
from helperFuncs import *


## ==================================================================================
## ---- USER OPTIONS ----------------------------------------------------------------
## ==================================================================================

gas_types = ['CO2']	# Gas choices: R236fa, R134a, N2, CO2, H2, air
								# Gas choice will determine geometry based on desired output that was determined in the course of this project
d_upstream = 2 / 1000			# Upstream "pipe" diameter (really just the valve orifice diameter), units of m (mm / 1000)
L_upstream = 20 / 1000			# Upstream "pipe length" aka path length through valve. Just an estimate.
T_wall = 293					# Valve brass body wall tempertaure used to evaluate heat transfer conductivity
m_brass	= 80 / 1000				# Mass brass, kg
cp_brass = 380					# Specific heat capacity of brass, J/kg-K
cutoff_cond = 0.0001			# Cutoff condition, defined by the fractional change in pressure (relative to P_t_init) per second, units of 1/sec
figsize = (7.5, 4)				# Figure size (in)
dpi = 150						# Figure dpi


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
		visc_func = create_visc_func_gas(fluid_props)
		cp_func = create_cp_func(fluid_props)
		cv_func = create_cv_func(fluid_props)
		ktc_func = create_ktc_func(fluid_props)
		h_from_PT_gas_func = create_h_from_PT_gas_func(fluid_props)
		u_from_PT_gas_func = create_u_from_PT_gas_func(fluid_props)

		P_from_rh_func, T_from_rh_func = create_PT_from_rh_gas_func(fluid_props_vol)
		P_from_ru_func, T_from_ru_func = create_PT_from_ru_gas_func(fluid_props_vol)


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
	rho_t_init = P_t_init/(Z*R*(T_t_init))  								# Initial density, kg/m^3 (or g/l), using GENERALIZED COMPRESSABILITY
	mu_t_init = 1/rho_t_init  												# Initial specific volume, m^3/kg (or l/g)
	m_init = rho_t_init*vol  												# Initial propellant mass, kg
	u_sp_init = u_from_PT_gas_func(P_t_init, T_t_init)[0]*1000				# Initial specific internal energy, J/kg
	h_sp_init = h_from_PT_gas_func(P_t_init, T_t_init)[0]*1000				# Initial specific enthalpy, J/kg
	dia = 2*(vol*(3/4)/np.pi)**(1/3)  										# Plenum diameter, m
	time_step_init = 15.28*vol*(P_t_init-P_amb)/((P_t_init*d_star)**2)		# Time step, s

	# Recalculate P_init based on rho_t_init and h_sp_init. I know it shouldn't be different, but it is, based on whichever function is being used (h_from_PT vs P_from_rh)
	P_t_init = P_from_rh_func(rho_t_init, h_sp_init/1000)[0]*1000000

	## ==================================================================================
	## ---- INIT DATA LISTS -------------------------------------------------------------
	## ==================================================================================
	time = [0]
	list_of_P_ts = [P_t_init] 
	list_of_T_ts = [T_t_init]
	list_of_rho_ts = [rho_t_init]
	list_of_mu_ts = [mu_t_init]
	list_of_u_sp = [u_sp_init]
	list_of_h_sp = [h_sp_init]
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

	list_of_average_thrusts = []
	cumulative_impulse = []
	ISP = []
	cumulative_mass = []

	list_of_flow_regimes = [] 
	list_of_area_ratios_at_shock = []

	list_of_P_fg_exit = []
	list_of_T_fg_exit = []


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
	while delta_pres > cutoff_cond and list_of_P_ts[-1] > P_amb:
		if thermal_model:
			# Put the CONSTANT PRESSURE temperature change here
			# First calculate Nusselt number using the Gnielinski correlation
			# ...to do that you need a Reynolds number, Prandtl number, and Darcy friction factor
			# ...to get the Reynolds number you need mass flow rate (nozzle function) and viscosity
			# ...viscosity and mass flow rate both depend on temperature
			# Realistically, the temperature will depend on the HTC and temperature gradient during this time period
			# For the sake of brevity, let's linearize it about the Isentropic Temperature and determine all the parameters at THAT temperature

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
			P_star, T_star, rho_star, Re_star, v_star, P_exit, T_exit, rho_exit, M_exit, v_exit, c_exit, m_dot, F, CF, flow_regime, area_ratio_at_shock, F_mdotv, F_pdiff = nozzle(k, R, M_crit_sub, M_crit_sup, list_of_P_ts[-1], list_of_T_ts[-1], list_of_rho_ts[-1], P_amb, d_star, expansion_ratio, half_angle, gas_type, visc_func)
			# Step 2:
			T_upstream = list_of_T_ts[-1]/Z_func(M_sub_upstream)
			P_upstream = list_of_P_ts[-1]/(Z_func(M_sub_upstream)**W)
			mu_upstream = visc_func(P_upstream, T_upstream)[0] / 1000000									# Convert from uPa-s to Pa-s (equivalent to kg/m-s), and pull value out of resultant array
			# Step 3:
			Re_upstream = m_dot*d_upstream/(A_upstream*mu_upstream)
			# Step 4:
			f = 0.07																						# Darcy friction factor. Just an estimate from Moody chart based on brass or steel or something similar.
			f_alt = (0.79*np.log(Re_upstream)-1.64)**-2
			cp_upstream = cp_func(P_upstream, T_upstream)[0] * 1000											# Convert from J/g-K to J/kg-K
			ktc_upstream = ktc_func(P_upstream,  np.average([T_wall, T_upstream]))[0]						# Thermal conductivity, W/m-K. This value is typically evaluated at the average of the wall temperature and fluid temperature
			# ktc_upstream_alt = ktc_func(P_upstream,  T_upstream)[0]
			Pr_upstream = mu_upstream * cp_upstream / ktc_upstream											# Calculate Prandtl number
			# Nu_upstream = (f/8)*(Re_upstream - 1000)*Pr_upstream / ( 1 + 12.7*np.sqrt(f/8)*(Pr_upstream**(2/3) - 1) )				# Gnielinsky correlation
			Nu_upstream_alt = (f_alt/8)*(Re_upstream - 1000)*Pr_upstream / ( 1 + 12.7*np.sqrt(f_alt/8)*(Pr_upstream**(2/3) - 1) )
			# Step 5:
			htc_upstream = Nu_upstream_alt * ktc_upstream / L_upstream										# Heat transfer coefficient, W/m^2-K
			# htc_upstream_alt = 0.023 * (ktc_upstream_alt / d_upstream) * (m_dot*d_upstream/(A_upstream*mu_upstream))**0.8 * (mu_upstream*cp_upstream/ktc_upstream_alt)**0.4		# Dittus-Boelter equation
			# Step 6:
			Q_dot_upstream = htc_upstream * (np.pi*d_upstream*L_upstream) * (T_wall - T_upstream)			# Heat rate, W (J/s)
			delh_upstream = Q_dot_upstream / m_dot															# Change in specific enthalpy, J/kg
			delT_wall = -Q_dot_upstream*time_step/m_brass*cp_brass
			T_wall += delT_wall
			# Step 7:
			delT_upstream = delh_upstream / cp_upstream
			# Step 8:
			T_upstream = T_upstream + delT_upstream
			list_of_T_ts[-1] = T_upstream * Z_func(M_sub_upstream)
			list_of_rho_ts[-1] = list_of_P_ts[-1] / (R * list_of_T_ts[-1])


		# Separate T_t from T_t_plenum (same with P_t)


		# Now go through and calculate nozzle stuff with this new T_t
		P_star, T_star, rho_star, Re_star, v_star, P_exit, T_exit, rho_exit, M_exit, v_exit, c_exit, m_dot, F, CF, flow_regime, area_ratio_at_shock, F_mdotv, F_pdiff = nozzle(k, R, M_crit_sub, M_crit_sup, list_of_P_ts[-1], list_of_T_ts[-1], list_of_rho_ts[-1], P_amb, d_star, expansion_ratio, half_angle, gas_type, visc_func)

		# Append returned items to their respective lists
		list_of_P_stars.append(P_star)									# Throat Pressure, Pa
		list_of_T_stars.append(T_star)									# Throat Temperature, K
		list_of_rho_stars.append(rho_star)								# Throat Density, kg/m^s
		list_of_Re_stars.append(Re_star)								# Throat Reynolds Number
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
		list_of_rho_ts.append(m_gas[-1]/(vol))
		list_of_mu_ts.append(1/list_of_rho_ts[-1])

		# From mass and energy balance, and assuming the flow is isentropic (is it? Yes I believe it is. This model assumes it is NOT isentropic across the c.v., but it is isentropic afterwards, hence the "upstream" relations made here).
		P_upstream = list_of_P_ts[-1]*Z_func(M_sub_upstream)**-W			# Pa
		T_upstream = list_of_T_ts[-1]*Z_func(M_sub_upstream)**-1			# K
		h_upstream = h_from_PT_gas_func(P_upstream, T_upstream)[0]*1000		# J/kg
		v_upstream = M_sub_upstream*np.sqrt(k*R*T_upstream)					# m/s
		dE_dt = m_dot*(h_upstream + (v_upstream**2)/2)						# J/s
		delta_E = dE_dt*time_step

		list_of_h_sp.append( ((list_of_h_sp[-1]*m_gas[-2]) - delta_E)/m_gas[-1] )		# This value won't be calculated correctly if you assume an isentropic process in the plenum. This should go down, but it in fact goes up.
		list_of_u_sp.append( ((list_of_u_sp[-1]*m_gas[-2]) - delta_E)/m_gas[-1] )		

		P_from_rh = P_from_rh_func(list_of_rho_ts[-1], list_of_h_sp[-1]/1000)[0]*1000000
		T_from_rh = T_from_rh_func(list_of_rho_ts[-1], list_of_h_sp[-1]/1000)[0]

		list_of_P_ts.append(P_from_rh)
		list_of_T_ts.append(T_from_rh)

		# I'd like another function to determine the phase based on (u,v)
		# I'd like a third function to determine quality based on the aforementioned parameters.

		# Isothermal process in plenum
		# list_of_T_ts.append(T_t_init)

		# Isentropic process in plenum
		# list_of_P_ts.append( list_of_P_ts[-1]*(list_of_rho_ts[-1]/list_of_rho_ts[-2])**k )
		# list_of_T_ts.append( list_of_T_ts[-1]*(list_of_rho_ts[-1]/list_of_rho_ts[-2])**(k-1) )					

		# Isenthalpic process in plenum
		# cv = cv_func(list_of_P_ts[-1], list_of_T_ts[-1])[0]*1000		# J/kg-K
		# list_of_P_ts.append( list_of_P_ts[-1]*(m_gas[-1]/m_gas[-2])**((cv-R)/(cv+R)) )
		# list_of_T_ts.append( list_of_T_ts[-1]*(m_gas[-1]/m_gas[-2])**(2*R/(cv+R)) )

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
									list_of_P_ts, 
									list_of_T_ts, 
									list_of_rho_ts,
									list_of_mu_ts,
									list_of_h_sp,
									list_of_u_sp,

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
									list_of_average_thrusts,
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
									'mu_t',
									'h_sp',
									'u_sp',

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
									'visc_losses',
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
	
	ax.set_title(r'{} Phase Diagram and Corresponding Plenum P-T Path'.format(gas_label))
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


all_parameters = all_parameters.set_index('gas_type')


## ==================================================================================
## ---- PLOT ------------------------------------------------------------------------
## ==================================================================================
# ---- Plot 2x3 [Thrust, Impulse, ISP, Re, Ma, Density] all vs. Inlet + Throat + Exit Pressure
linewidth = 2
fontsize = 12

data = 	{ 
			'P_t': 				all_data['P_t'],
			'T_t': 				all_data['T_t'],
			'rho_t':			all_data['rho_t'],
			# 'mu_t':				all_data['mu_t'],
			'h_sp':				all_data['h_sp'],
			# 'u_sp':				all_data['u_sp'],

			# 'P_star': 		all_data['P_star'],
			# 'T_star': 		all_data['T_star'],
			# 'rho_star': 		all_data['rho_star'],
			# 'Re_star': 		all_data['Re_star'],
			# 'v_star': 		all_data['v_star'],

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
			'mu_t': 			'Total Specific Volume',
			'h_sp':				'Specific Enthalpy',
			'u_sp':				'Specific Internal Energy',

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
			'mu_t':				all_data['time'],
			'h_sp':				all_data['time'],
			'u_sp':				all_data['time'],

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
			'T_t': 				'Total Temperature, $K$', 						
			'rho_t': 			'Total Density, $kg/m^3$',
			'mu_t': 			'Total Specific Volume, $m^3/kg$',
			'h_sp':				'Specific Enthalpy, $kJ/kg$',
			'u_sp':				'Specific Internal Energy, $kJ/kg$',

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

for gas in gas_types:
	fig_sat, axs = plt.subplots(4,1, figsize=figsize, dpi=dpi, sharex=True)
	fig_sat.suptitle('Exit Pressure & Temperature w/ Saturation Data', y=0.98)
	axs[0].set_title(r'({} at $T_0$={} K, $V_{{p}}=${} cm$^3$, Nozzle $\varnothing${} mm, $\lambda$={})'.format(all_parameters.loc[gas]['gas_label'], all_parameters.loc[gas]['T_t_init'], all_parameters.loc[gas]['vol']*10**6, all_parameters.loc[gas]['d_star']*1000, all_parameters.loc[gas]['expansion_ratio']), fontsize=9)
	fig_sat.canvas.set_window_title('Saturated Pressure Stuff')

	for i, key in enumerate(data):
		sns.lineplot(ax=axs[i], x='time', y=key, palette='colorblind', data=all_data[all_data['gas_type']==gas], hue='flow regimes', legend='full')

		if key=='P_exit':
			# sns.lineplot(ax=axs[i], x=times['P_fg_exit'], y=all_data['P_fg_exit'], palette='colorblind', data=all_data[all_data['gas_type']==gas], legend='full')
			axs[i].plot(all_data[all_data['gas_type']==gas]['time'], all_data[all_data['gas_type']==gas]['P_fg_exit'], color='red', label='Phase Change Pres at Exit Temp')
			axs[i].plot(all_data[all_data['gas_type']==gas]['time'], all_data[all_data['gas_type']==gas]['P_trip'], color='green', linestyle='--', label='Triple Point')
		if key=='T_exit':
			# sns.lineplot(ax=axs[i], x=times['T_fg_exit'], y=all_data['T_fg_exit'], palette='colorblind', data=all_data[all_data['gas_type']==gas], legend='full')
			axs[i].plot(all_data[all_data['gas_type']==gas]['time'], all_data[all_data['gas_type']==gas]['T_fg_exit'], color='red', label='Phase Change Temp at Exit Pres')
			axs[i].plot(all_data[all_data['gas_type']==gas]['time'], all_data[all_data['gas_type']==gas]['T_trip'], color='green', linestyle='--', label='Triple Point')

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
