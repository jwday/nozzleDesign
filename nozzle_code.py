# 1D_Isentropic_Nozzle_Flow Code
import numpy as np
import math
import matplotlib.pyplot as plt
import sympy as sm
import scipy.optimize as opti

def nozzle(P_0, T_0, P_amb, d_star, expansion_ratio, half_angle, gas_type):

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
        R = 8.314/0.04401


    ## --------------------------------------------------------------------------------


    # Unit Conversion
    d_star /= 1000  # Throat Diameter (m)
    P_0 *= (137.9/20)*1000  # Total initial pressure (Pa)
    T_0 += 273.15  # Total temperature (K)
    rho_0 = P_0/(R*T_0)  # Total density (kg/m^3)
    P_amb *= (137.9/20)*1000  # Ambient pressure (Pa)
    a_0 = np.sqrt(k*R*T_0)

    # Critical (Sonic) Conditions (Isentropic)
    A_star = np.pi*(d_star**2)/4  # Throat area
    # T_star = T_0/(1 + (k-1)/2)
    T_star = T_0*(2/(k+1))
    # P_star = P_0*(T_0/T_star)**(-k/(k-1))
    P_star = P_0*(2/(k+1))**(k/(k-1))
    rho_star = P_star/(R*T_star)
    a_star = np.sqrt(k*R*T_star)
    mu_star = (2*(10**-9)*(T_star**3) - 5*(10**-6)*(T_star**2) + 0.007*T_star + 0.0613)*(10**-5)

    # Exit Conditions
    A_exit = A_star*expansion_ratio  # Exit area
    d_exit = np.sqrt(4*A_exit/np.pi)  # Exit diameter (m)


    ## --------------------------------------------------------------------------------

    # Use the all-important Isentropic Area-Mach Number relation to solve for the Mach Number at the exit
    P = (k+1)/2
    L = (k-1)/2
    W = k/(k-1)
    Q = P/L  # aka (k+1)/(k-1)
    S = (A_exit/A_star)**2

    def objective(X):
        f = (1 + L*X)**Q - S*X*(P**Q)  # Where 'X' is M^2
        return f

    x0 = np.array([0.00001, 20])
    sol = opti.fsolve(objective, x0, maxfev=100000, full_output=False, xtol=0.000000001)

    M_exit_sub = np.sqrt(sol[0])  # Subsonic mach no.
    M_exit_sup = np.sqrt(sol[1])  # Supersonic mach no.


    ## --------------------------------------------------------------------------------
    

    Y = np.sqrt(k/R)  # Just another component in the mass flow rate equation
    m_dot = ( A_star*P_0/np.sqrt(T_0) )*Y*( P )**(-Q/2)  # Mass flow rate (kg/s)
    Q_scfm = (m_dot/1.98)*35.3147*60  # Volumetric flow rate (m^3/s)
    # CO2 is 1.98 kg/m^3 at STP
    # 35.3147 ft^3 / m^3

    # Determine sub and supersonic exit conditions
    T_exit = []
    P_exit = []
    v_exit = []
    a_exit = []

    for i in [M_exit_sub, M_exit_sup]:
        Z = 1 + L*(i**2)
        T_exit.append(T_0/Z)  # Exit temperature (K)
        v_exit.append( i*np.sqrt(k*R*(T_0/Z)) )  # Exit velocity (m/s)
        a_exit.append(np.sqrt(k*R*(T_0/Z)))  # Exit speed of sound (m/s)
        P_exit.append(P_0*Z**(-W))

        # if P_exit[i] > P_amb:
        #     P_exit.append(P_0*Z**(-W))  # Exit pressure (Pa)
        # else:
        #     P_exit.append(P_amb)
        #     P_total_shock = P_amb*( (1 + L*M_exit_sub**2)**(k/(k-1)) )



    # Determine thrust
    correction = (1 + np.cos(np.radians(half_angle)))/2
    F = correction*(m_dot*v_exit[1]) + (P_exit[1]-P_amb)*A_exit  # Thrust (N)
    F_mdotv = correction*m_dot*v_exit[1]
    F_pdiff = (P_exit[1]-P_amb)*A_exit

    # Reynold's Number at the throat
    Re_star = rho_star*a_star*d_star/mu_star



    return m_dot, P_star, T_star, a_star, Re_star, M_exit_sub, M_exit_sup, P_exit, T_exit, v_exit, a_exit, F, F_mdotv, F_pdiff