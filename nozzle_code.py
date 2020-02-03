# 1D_Isentropic_Nozzle_Flow Code
import numpy as np
import sys
import math
import matplotlib.pyplot as plt
import sympy as sm
import scipy.optimize as opti

def nozzle(P_t, T_t, P_amb, d_throat, expansion_ratio, half_angle, gas_type):

	if gas_type == 'R236fa':
		k = 1.083  # Heat capacity ratio (Cp/Cv)
		R = 8.314/0.152039  # Specific gas constant (J/kg-K)

	if gas_type == 'R134a':
		k = 1.127  # Override value to compare results to MST paper
		R = 8.314/0.10203  # Specific gas constant (J/kg-K)

	if gas_type == 'N2':
		k = 1.039
		R = 8.314/0.028014  # Specific gas constant (J/kg-K)
		
	if gas_type == 'CO2':
		k = 1.289
		R = 8.314/0.04401  # Specific gas constant (J/kg-K)
	
	if gas_type == 'air':
		k = 1.401
		R = 8.314/0.0289645  # Specific gas constant (J/kg-K)


	## --------------------------------------------------------------------------------
	rho_t = P_t/(R*T_t)  # Total density (kg/m^3)
	c_0 = np.sqrt(k*R*T_t)

	# Critical (Sonic) Conditions (Isentropic)
	A_throat = np.pi*(d_throat**2)/4  # Throat area
	T_throat = T_t*(2/(k+1))
	P_throat = P_t*(2/(k+1))**(k/(k-1)) # Pa
	rho_throat = rho_t*(2/(k+1))**(1/(k-1)) # kg/m^3
	c_throat = np.sqrt(k*R*T_throat) # m/s
	mu_throat = (0.0492*T_throat + 0.3276)*(10**-6) # Pa-s, Emperical formula derived from NIST data. Viscosity is highly invariable with pressure between 0.1 and 0.8 MPa, and only somewhat variable with temperature, and linearly at that. Use 13 uPa-s in a pinch.
	Re_throat = rho_throat*c_throat*d_throat/mu_throat

	# Exit Conditions
	A_exit = A_throat*expansion_ratio  # Exit area
	d_exit = np.sqrt(4*A_exit/np.pi)  # Exit diameter (m)


	## --------------------------------------------------------------------------------
	# Establish thermodynamic relationships in more compact terms
	P = (k+1)/2
	L = (k-1)/2
	W = k/(k-1)
	Q = P/L  # aka (k+1)/(k-1)
	S = (A_exit/A_throat)**2
	Y = np.sqrt(k/R)
	
	def pres_ratio(M):
		f = (1 + L*M**2)**(-W)  # P/P_t
		return f


	## --------------------------------------------------------------------------------
	# Solve the Isentropic Area-Mach Number relation for the critical Mach Numbers at the exit
	# NOTE: This equation is only valid for fully subsonic OR fully supersonic flow throughout the entire nozzle
	# It is not valid for the regime which includes a normal shock within the nozzle
	def objective(X):
		f = (1 + L*X)**Q - S*X*(P**Q)  # Where 'X' is M^2
		return f

	x0 = np.array([0.00001, 20])
	sol0 = opti.fsolve(objective, x0, maxfev=100000, full_output=False, xtol=0.000000001)

	M_crit_sub = np.sqrt(sol0[0])  # Subsonic critical mach no.
	M_crit_sup = np.sqrt(sol0[1])  # Supersonic critical mach no.

	# Use the above solutions to determine the corresponding pressure ratios that produce them
	PR_crit_sub = ( 1 + L*M_crit_sub**2 )**W  # (P_t/P_amb)_ea
	PR_crit_sup = ( 1 + L*M_crit_sup**2 )**W  # (P_t/P_amb)_eb

	M_exit_behindshock = np.sqrt( (1 + L*M_crit_sup**2)/(k*M_crit_sup**2 - L) )  # Should be < 1
	PR_exit_shock = np.sqrt( S*(P**Q)*M_exit_behindshock*(1 + L*M_exit_behindshock**2) )  # (P_t/P_amb)_exit_shock  

	## --------------------------------------------------------------------------------
	# Entirely subsonic. SONIC CONDITIONS DO NOT APPLY AND THEREFORE THROAT STATE MUST BE RECALCULATED.
	if P_t/P_amb > 1 and P_t/P_amb <= PR_crit_sub:
		M_exit = np.sqrt( (1/L)*((P_t/P_amb)**(1/W) - 1) )  # Should be < 1
		# Exit pressure = Ambient pressure. Any and all thrust is purely from momentum exchange.
		P_exit = P_amb




	# Shock somewhere in the flow (Overexpansion)
	elif P_t/P_amb > PR_crit_sub and P_t/P_amb < PR_crit_sup:
		# The flow is overexpanded and a shock will occur, but it's anyone's guess as to whether the shock occurs inside the nozzle or outside
		# We would need to know where the shock occurs based on the area ratio at the location of the shock
		# The only thing we do know is that the flow is isentropic before and after the shock, but not between, and the exit pressure = P_ambient
		# The total pressure decreases across the shock

		# Shock! But where? Start by checking the extreme case of a shock at the nozzle exit.
		# look above

		if P_t/P_amb < PR_exit_shock:
			# Shock occurs inside the nozzle
			def shock_function(X):
				f = ((P_t/P_amb)**2)/S - (P**Q)*X*(1 + L*X)  # Where 'X' is M^2
				return f

			x1 = np.array([0.00001, 20])
			sol1 = opti.fsolve(shock_function, x1, maxfev=100000, full_output=False, xtol=0.000000001)
			
			M_exit = np.sqrt(sol1[0])
			# Exit pressure = Ambient pressure and subsonic. Any and all thrust is purely from momentum exchange.
			P_exit = P_amb
		
		elif P_t/P_amb == PR_exit_shock:
			M_exit = M_exit_behindshock
			# Exit pressure = Ambient pressure and subsonic. Any and all thrust is purely from momentum exchange.
			P_exit = P_amb

		elif P_t/P_amb > PR_exit_shock:
			M_exit = M_crit_sup
			# Exit pressure < Ambient pressure but not enough to cause it to shock before the exit. All thrust is from momentum exchange and pressure differential.
			P_exit = pres_ratio(M_exit)*P_t

		else:
			sys.exit("Some sort of error occured when trying to figure out the exit conditions for a nozzle shock")
	
	# No shock in the flow (Perfect expansion)
	elif P_t/P_amb == PR_crit_sup:
		M_exit = M_crit_sup
		# Exit pressure = Ambient pressure and supersonic. All thrust is from momentum exchange.
		P_exit = P_amb

	# No shock in the flow (Underexpansion)
	elif P_t/P_amb > PR_crit_sup:
		M_exit = M_crit_sup
		# Exit pressure > Ambient pressure. All thrust is from momentum exchange and pressure differential.
		P_exit = pres_ratio(M_exit)*P_t

	# No flow at all
	else:
		M_exit = 0  # Because the only other possible condition left is if P_t = P_amb
		P_exit = P_amb


	## --------------------------------------------------------------------------------    
	def Z_func(M):
		Z = 1 + L*M**2
		return Z
	
	T_exit = T_t/Z_func(M_exit)  # Exit temperature (K)
	a_exit = np.sqrt(k*R*T_exit)  # Exit speed of sound (m/s)
	v_exit = M_exit*a_exit  # Exit velocity (m/s)
	rho_exit = (P_exit*(rho_t**k)/P_t)**(1/k) # Units of kg/m^3
	m_dot = rho_exit*A_exit*v_exit  # Units of kg/s

	if P_t/P_amb > 1 and P_t/P_amb <= PR_crit_sub:
		def objective2(X):
			f = m_dot - (A_throat*P_t/math.sqrt(T_t)) * math.sqrt(k/R) * X * (1 + L*X**2)**(-Q/2)
			return f

		x0 = np.array([0.3])
		sol0 = opti.fsolve(objective2, x0, maxfev=100000, full_output=False, xtol=0.000000001)
		M_throat = sol0[0]

		P_throat = P_t*Z_func(M_throat)**-W
		rho_throat = rho_t*Z_func(M_throat)**-(W/k)
		T_throat = T_t*Z_func(M_throat)**(-1)
		c_throat = math.sqrt(k*R*T_throat)
		v_throat = M_throat*c_throat
		mu_throat = (0.0492*T_throat + 0.3276)*(10**-6) # Pa-s
		Re_throat = rho_throat*v_throat*d_throat/mu_throat

	# Q_scfm = (m_dot/1.98)*35.3147*60  # Volumetric flow rate (m^3/s)
	# CO2 is 1.98 kg/m^3 at STP
	# 35.3147 ft^3 / m^3

	# Determine sub and supersonic exit conditions
	# T_exit = []
	# P_exit = []
	# v_exit = []
	# c_exit = []

	# for i in [M_crit_sub, M_crit_sup]:
	#     Z = 1 + L*(i**2)
	#     T_exit.append(T_t/Z)  # Exit temperature (K)
	#     v_exit.append( i*np.sqrt(k*R*(T_t/Z)) )  # Exit velocity (m/s)
	#     c_exit.append(np.sqrt(k*R*(T_t/Z)))  # Exit speed of sound (m/s)
	#     P_exit.append(P_t*Z**(-W))

		# if P_exit[i] > P_amb:
		#     P_exit.append(P_t*Z**(-W))  # Exit pressure (Pa)
		# else:
		#     P_exit.append(P_amb)
		#     P_total_shock = P_amb*( (1 + L*M_crit_sub**2)**(k/(k-1)) )



	# Determine thrust
	correction = (1 + np.cos(np.radians(half_angle)))/2
	F_mdotv = correction*m_dot*v_exit
	F_pdiff = (P_exit-P_amb)*A_exit
	F = F_mdotv + F_pdiff  # Thrust (N)


	# Reynold's Number at the throat
	# Re_throat = rho_throat*c_throat*d_throat/mu_throat



	# return m_dot, P_throat, T_throat, c_throat, Re_throat, M_crit_sub, M_crit_sup, P_exit, T_exit, v_exit, F, F_mdotv, F_pdiff
	return m_dot, M_crit_sub, M_crit_sup, PR_crit_sub, PR_crit_sup, PR_exit_shock, M_exit_behindshock, M_exit, P_exit, v_exit, F, P_throat, T_throat, rho_throat, Re_throat, T_exit, rho_exit