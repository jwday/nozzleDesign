# 1D_Isentropic_Nozzle_Flow Code
import numpy as np
import sys
import math
import matplotlib.pyplot as plt
import sympy as sm
import scipy.optimize as opti

def nozzle(P_t, T_t, P_amb, d_star, expansion_ratio, half_angle, gas_type):

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


    # Unit Conversion
    d_star /= 1000  # Throat Diameter (m)
    P_t *= (137.9/20)*1000  # Total initial pressure (Pa)
    # T_t += 273.15  # Total temperature (K)
    rho_t = P_t/(R*T_t)  # Total density (kg/m^3)
    P_amb *= (137.9/20)*1000  # Ambient pressure (Pa)
    a_0 = np.sqrt(k*R*T_t)

    # Critical (Sonic) Conditions (Isentropic)
    A_star = np.pi*(d_star**2)/4  # Throat area
    # T_star = T_t/(1 + (k-1)/2)
    T_star = T_t*(2/(k+1))
    # P_star = P_t*(T_t/T_star)**(-k/(k-1))
    P_star = P_t*(2/(k+1))**(k/(k-1))
    rho_star = P_star/(R*T_star)
    a_star = np.sqrt(k*R*T_star)
    mu_star = (2*(10**-9)*(T_star**3) - 5*(10**-6)*(T_star**2) + 0.007*T_star + 0.0613)*(10**-5)

    # Exit Conditions
    A_exit = A_star*expansion_ratio  # Exit area
    d_exit = np.sqrt(4*A_exit/np.pi)  # Exit diameter (m)


    ## --------------------------------------------------------------------------------
    # Establish thermodynamic relationships in more compact terms
    P = (k+1)/2
    L = (k-1)/2
    W = k/(k-1)
    Q = P/L  # aka (k+1)/(k-1)
    S = (A_exit/A_star)**2
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
    # Find exit mach number for entirely subsonic
    if P_t/P_amb > 1 and P_t/P_amb <= PR_crit_sub:
        # Entirely subsonic.
        M_exit = np.sqrt( (1/L)*((P_t/P_amb)**(1/W) - 1) )  # Should be < 1
        # Exit pressure = Ambient pressure. Any and all thrust is purely from momentum exchange.
        P_exit = P_amb



    # Find exit mach number for a shock somewhere within the flow
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
        
    elif P_t/P_amb == PR_crit_sup:
        M_exit = M_crit_sup
        # Exit pressure = Ambient pressure and supersonic. All thrust is from momentum exchange.
        P_exit = P_amb

    elif P_t/P_amb > PR_crit_sup:
        M_exit = M_crit_sup
        # Exit pressure > Ambient pressure. All thrust is from momentum exchange and pressure differential.
        P_exit = pres_ratio(M_exit)*P_t

    else:
        M_exit = 0  # Because the only other possible condition left is if P_t = P_amb
        P_exit = P_amb


    ## --------------------------------------------------------------------------------    
    Z = 1 + L*M_exit**2
    T_exit = T_t/Z  # Exit temperature (K)
    c_exit = np.sqrt(k*R*T_exit)  # Exit speed of sound (m/s)
    v_exit = M_exit*c_exit  # Exit velocity (m/s)
    rho_exit = P_exit/(R*T_exit)  # Units of kg/m^3
    m_dot = rho_exit*A_exit*v_exit  # Units of kg/s


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
    # Re_star = rho_star*a_star*d_star/mu_star



    # return m_dot, P_star, T_star, a_star, Re_star, M_crit_sub, M_crit_sup, P_exit, T_exit, v_exit, F, F_mdotv, F_pdiff
    return m_dot, M_crit_sub, M_crit_sup, PR_crit_sub, PR_crit_sup, PR_exit_shock, M_exit_behindshock, M_exit, P_exit/6894.76, v_exit, F, P_star/6894.76, T_star, rho_star, T_exit, rho_exit