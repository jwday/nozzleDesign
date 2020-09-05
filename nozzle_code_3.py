# 1D_Isentropic_Nozzle_Flow Code
# This code has been modified as of 4/13/2020 to allow for P_amb = 0 (in-space conditions) by inverting the conditions and functions
import numpy as np
import sys
import math
import pandas as pd
import scipy.optimize as opti
from scipy.interpolate import interp1d
from scipy.interpolate import interp2d

def nozzle(k, R, M_crit_sub, M_crit_sup, P_t, T_t, rho_t, P_amb, d_star, expansion_ratio, half_angle, gas_type, visc_func, r_from_PT):
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
	def area_ratio(M):
		f = (P**(-Q/2)) * (Z_func(M)**(Q/2))/M
		return f

	# Calculate relevant pressure ratios
	PR_crit_sub = 1/(( 1 + L*M_crit_sub**2 )**W)  	# Entirely subsonic if (P_amb/P_t)_ea < P_amb/P_t < 1
	PR_crit_sup = 1/(( 1 + L*M_crit_sup**2 )**W)  	# Entirely supersonic if P_amb/P_t < (P_amb/P_t)_eb
	def pres_ratio(M):
		f = (1 + L*M**2)**(-W)  					# P/P_t, Isentropic Mach-Pressure relation.
		return f

	# Calculate compression ratio
	Z = P_t / (rho_t * R * T_t)


	## ==================================================================================
	## ---- SOLVE FLOW MACH NOS. BASED ON PR --------------------------------------------
	## ==================================================================================
	shock_in_nozzle_flag = False		# Default
	flow_is_supersonic_flag = True		# Default
	area_ratio_at_shock = None			# Default

	# Case 1: Entirely subsonic. SONIC CONDITIONS DO NOT APPLY AND THEREFORE THROAT STATE MUST BE RECALCULATED.
	if P_amb/P_t < 1 and P_amb/P_t >= PR_crit_sub:
		M_exit = np.sqrt( (1/L)*((P_amb/P_t)**(-1/W) - 1) )
		P_exit = P_amb
		
		# Since it's entirely subsonic, you gotta go about things a little differently...
		# In order to calculate the subsonic throat Mach number, you need mdot.
		# But to get m_dot we need to use the exit conditions since we haven't calculated the throat conditions yet.
		T_exit = T_t/Z_func(M_exit)  				# Exit temperature (K)
		c_exit = np.sqrt(k*R*T_exit)  				# Exit speed of sound (m/s)
		v_exit = M_exit*c_exit  					# Exit velocity (m/s)
		rho_exit = (P_exit*(rho_t**k)/P_t)**(1/k) 	# Units of kg/m^3
		m_dot = Z*rho_exit*A_exit*v_exit  			# Units of kg/s

		area = []
		machs = list(np.linspace(0, 1, 100))
		for mach in machs:
			area.append(m_dot / ( (P_t/np.sqrt(T_t)) * np.sqrt(k/R) * mach * (1 + L*(mach**2))**(-Q/2) ))
		subsonic_mach_from_area = interp1d(area, machs)

		# def objective2(X):
		# 	f = m_dot - (A_star*P_t/math.sqrt(T_t)) * math.sqrt(k/R) * np.sqrt(X) * (1 + L*X)**(-Q/2)
		# 	return f
		# x0 = np.array([0.00001])
		# sol0 = opti.fsolve(objective2, x0, maxfev=100000, full_output=False, xtol=0.000000001)
		# M_star = np.sqrt(sol0[0])
		M_star = subsonic_mach_from_area(A_star)
		flow_is_supersonic_flag = False
		flow_regime = 'Subsonic'
		

	## ----------------------------------------------------------------------------------
	# Case 2: Shock somewhere in the flow (Overexpansion)
	elif P_amb/P_t < PR_crit_sub and P_amb/P_t > PR_crit_sup:
		# Flow is isentropic before and after the shock, but not across it.
		# The total pressure decreases across the shock.
		# Must solve (P_t)(A_star) = (P_t_exit*)(A_exit)f(M) to find Mach number at exit.
		M_star = 1
		M_exit_behindshock = np.sqrt( (1 + L*M_crit_sup**2) / (k*M_crit_sup**2 - L) )  					# Mach-shock relation (relates mach numbers just before and just after normal shock). This is the mach number just after the shock if the shock is at the exit.
		PR_exit_shock = 1/(np.sqrt( (S*(P**Q)) * (M_exit_behindshock**2) * Z_func(M_exit_behindshock) ))  	# Pressure ratio required to have a shock sit exactly at the nozzle exit, (P_amb/P_t)_exit_shock

		# Case 2a: Shock occurs inside the nozzle. Exit pressure = ambient and exit flow is subsonic.
		if P_amb/P_t > PR_exit_shock:
			def shock_function(X):
				f = ((P_t/P_amb)**2)/S - (P**Q)*X*(1 + L*X)  # Where 'X' is M^2
				return f
			x1 = np.array([0.00001, 20])
			sol1 = opti.fsolve(shock_function, x1, maxfev=100000, full_output=False, xtol=0.000000001)
			M_exit = np.sqrt(sol1[0])
			P_exit = P_amb

			# This block computes the location (in terms of area ratio) of the shock that presumably exists inside the nozzle
			P_t_0 = P_t
			P_t_1 = P_exit * (1 + L*(M_exit**2))**W
			def normal_shock_pressure_ratio(X):
				f = ((( (P*X)/(1 + L*X) )**W) * ((P/(k*X - L))**(W/k))) - P_t_1/P_t_0	# Where 'X' is M^2
				return f
			x2 = np.array([5])
			sol2 = opti.fsolve(normal_shock_pressure_ratio, x2, maxfev=100000, full_output=False, xtol=0.000000001)
			M_just_before_shock = np.sqrt(sol2[0])
			area_ratio_at_shock = area_ratio(M_just_before_shock)
			shock_in_nozzle_flag = True
			flow_regime = 'Normal Shock in Nozzle'

		# Case 2b: Shock occurs exactly at nozzle exit. Exit pressure = ambient and exit flow is subsonic.
		elif P_amb/P_t == PR_exit_shock:
			M_exit = M_exit_behindshock
			P_exit = P_amb
			flow_regime = 'Normal Shock'

		# Case 2c: Shock occurs just outside the nozzle. Exit pressure < ambient and exit flow is supersonic.
		elif P_amb/P_t < PR_exit_shock:
			M_exit = M_crit_sup
			P_exit = pres_ratio(M_exit)*P_t
			flow_regime = 'Overexpanded'

		# Catch-all. Probably logically redundant. 
		else:
			sys.exit("Some sort of error occured when trying to figure out the exit conditions for a nozzle shock")
	

	## ----------------------------------------------------------------------------------
	# Case 3: No shock in the flow (Perfect expansion)
	elif P_amb/P_t == PR_crit_sup:
		# Exit pressure = Ambient pressure and supersonic.
		M_star = 1
		M_exit = M_crit_sup
		P_exit = P_amb
		flow_regime = 'Perfectly Expanded'


	## ----------------------------------------------------------------------------------
	# Case 4: No shock in the flow (Underexpansion)
	elif P_amb/P_t < PR_crit_sup:
		# Exit pressure > Ambient pressure.
		# If ambient is vaccuum then this will always be true
		M_star = 1
		M_exit = M_crit_sup
		P_exit = pres_ratio(M_exit)*P_t
		flow_regime = 'Underexpanded'


	## ----------------------------------------------------------------------------------
	# Case 5: No flow at all (Catch-all)
	else:
		P_exit = P_amb
		M_star = 0
		M_exit = 0
		m_dot = 0
		T_exit = T_t
		c_exit = np.sqrt(k*R*T_exit)
		v_exit = 0
		rho_exit = rho_t
		flow_is_supersonic_flag = False
		flow_regime = 'No Flow'




	## ==================================================================================
	## ---- SOLVE FOR THROAT CONDITIONS -------------------------------------------------
	## ==================================================================================
	T_star = T_t*(2/(k+1))							# Throat Static Temperature, K
	P_star = P_t*(2/(k+1))**(k/(k-1)) 				# Throat Static Pressure, Pa
	rho_star = rho_t*(2/(k+1))**(1/(k-1)) 			# kg/m^3, P_t, T_t, rho_t all supplied from input

	c_star = np.sqrt(k*R*T_star) 					# m/s
	v_star = M_star*c_star							# m/s
	visc_star = visc_func(P_star, T_star)[0] / 1000000
	Re_star = rho_star*v_star*d_star/visc_star

	if flow_is_supersonic_flag:
		m_dot = Z*rho_star*A_star*v_star					# kg/s
		# m_dot = (A_star*P_t/math.sqrt(T_t)) * math.sqrt(k/R) * M_star * (Z_func(M_star))**(-Q/2)

		## ==================================================================================
		## ---- SOLVE FOR EXIT CONDITIONS ---------------------------------------------------
		## ==================================================================================
		T_exit = T_t/Z_func(M_exit)  				# Exit Static Temperature (K)
		c_exit = np.sqrt(k*R*T_exit)  				# Exit speed of sound (m/s)
		v_exit = M_exit*c_exit  					# Exit velocity (m/s)
		rho_exit = (P_exit*(rho_t**k)/P_t)**(1/k) 	# Units of kg/m^3

	else:
		pass





	# Determine thrust
	correction = (1 + np.cos(np.radians(half_angle)))/2
	F_mdotv = correction*m_dot*v_exit
	F_pdiff = (P_exit-P_amb)*A_exit
	F = F_mdotv + F_pdiff  # Thrust (N)
	CF = F/(P_t * A_star)


	# return a whole lot of stuff
	return P_star, T_star, rho_star, Re_star, M_star, v_star, P_exit, T_exit, rho_exit, M_exit, v_exit, c_exit, m_dot, F, CF, flow_regime, area_ratio_at_shock, F_mdotv, F_pdiff