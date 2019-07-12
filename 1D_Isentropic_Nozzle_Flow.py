# 1D_Isentropic_Nozzle_Flow
import numpy as np
import math
import matplotlib.pyplot as plt
import sympy as sm
import scipy.optimize as opti
P_0 = 64.7*(137.9/20)  # Total/Stagnation/Chamber Pressure (kPa)
T_0 = 10  # Total/Stagnation/Chamber Temperature (C)
P_amb = 14.7*(137.9/20)  # Ambient Pressure (kPa)

d_star = 0.5  # Throat diameter (mm) (1/64" = 0.397mm)
expansion_ratio = 1.3225  # 1.8048 for ideal expansion at 114.7 psi supply, 2.2447 for 164.7, 1.3225 for 0.2mm and 64.7 psi, 1.1235 for 0.3 and 44.7 psi
half_angle = 10  # Conical Nozzle expansion angle (degrees)

## ---- OPTIONS ----------------------------------

# Choices: R236fa, R134a, N2, CO2
gas_type = 'CO2'

plot_relationships = False
attempt_to_solve = False


## ------------------------------

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

# --------------------------------------------------------------------------------

def exit(): 
    raise Exception("Found exit()")
    

# Unit Conversion
d_star /= 1000  # Throat Diameter (m)
P_0 *= 1000  # Total pressure (Pa)
T_0 += 273.15  # Total temperature (K)
rho_0 = P_0/(R*T_0)  # Total density (kg/m^3)
P_amb *= 1000  # Ambient pressure (Pa)
a_0 = np.sqrt(k*R*T_0)



# Critical (Sonic) Conditions
A_star = np.pi*(d_star**2)/4  # Throat area
T_star = T_0/(1 + (k-1)/2)
P_star = P_0*(T_0/T_star)**(-k/(k-1))
rho_star = P_star/(R*T_star)
a_star = np.sqrt(k*R*T_star)
mu_star = (2*(10**-9)*(T_star**3) - 5*(10**-6)*(T_star**2) + 0.007*T_star + 0.0613)*(10**-5)


# Exit Conditions
A_exit = A_star*expansion_ratio  # Exit area
d_exit = np.sqrt(4*A_exit/np.pi)  # Exit diameter (m)

print('')


## ---- ISENTROPIC RELATIONSHIPS -----------------

P = (k+1)/2
L = (k-1)/2
Q = P/L
S = (A_exit/A_star)**2

def objective(X):
    f = (1 + L*X)**Q - S*X*(P**Q)
    return f

x0 = np.array([0.00001, 20])
sol = opti.fsolve(objective, x0, maxfev=100000, full_output=False, xtol=0.000000001)

M_exit_sub = np.sqrt(sol[0])  # Subsonic exit mach no.
M_exit_sup = np.sqrt(sol[1])  # Supersonic exit mach no.

M = np.arange(0.001, 5, 0.001)
E = M**2

area_ratio = np.sqrt( (E*P**Q)/((1 + L*E)**Q) )  # A*/A
diam_ratio = np.sqrt(1/area_ratio)
temp_ratio = 1/(1 + ((k-1)/2)*M**2)  # T/T_0
pres_ratio = (1 + ((k-1)/2)*M**2)**(-k/(k-1))  # P/P_0
a_ratio = np.sqrt(temp_ratio)  # a/a_0

def objective2(X, area):
    S_inst = (area/A_star)**2
    f = (1 + L*X)**Q - S_inst*X*(P**Q)
    return f

def mach_from_area(area):
    x0 = np.array([0.00001, 20])  # Starting estimate for the roots of objective2(X)
    sol = opti.fsolve(objective2, x0, args=(area), maxfev=100000, full_output=False, xtol=0.000000001)

    M_sub = np.sqrt(sol[0])  # Subsonic solution
    M_sup = np.sqrt(sol[1])  # Supersonic solution
    return M_sub, M_sup

def flow_properties(Mach):
    temp = T_0/(1 + ((k-1)/2)*Mach**2)
    pres = P_0*(1 + ((k-1)/2)*Mach**2)**(-k/(k-1))
    a = np.sqrt(k*R*temp)
    vel = Mach*a
    return temp, pres, a, vel

exit_area = (1/expansion_ratio)*np.ones(4999)  # A_exit/A*
freestream_pres = (P_amb/P_0)*np.ones(4999)  # P_amb/P*

## -----------------------------------------------

Y = np.sqrt(k/R)  # Just another component in the mass flow rate equation
m_dot = ( A_star*P_0/np.sqrt(T_0) )*Y*( P )**(-Q/2)  # Mass flow rate (kg/s)
Q_scfm = (m_dot/1.98)*35.3147*60
# CO2 is 1.98 kg/m^3 at STP
# 35.3147 ft^3 / m^3


Z = 1 + L*(M_exit_sup**2)

# Determine exit conditions
T_exit = (T_0/Z)  # Exit temperature (K)
P_exit = (P_0*Z**(-k/(k-1)))  # Exit pressure (Pa)
v_exit = ( M_exit_sup*np.sqrt(k*R*T_exit) )  # Exit velocity (m/s)
a_exit = np.sqrt(k*R*T_exit)  # Exit speed of sound (m/s)

# Determine thrust
correction = (1 + np.cos(np.radians(half_angle)))/2
F = correction*(m_dot*v_exit) + (P_exit-P_amb)*A_exit  # Thrust (N)
F_mdotv = correction*m_dot*v_exit
F_pdiff = correction*(P_exit-P_amb)*A_exit

# Reynold's Number at the throat
Re_star = rho_star*a_star*d_star/mu_star

# --------------------------------------------------------------------------------




# ----- I Want to plot Critical/Instantaneous Ratios -----------------------------

pres_crit_ratio = P_star/(P_0*pres_ratio)  # P*/P
temp_crit_ratio = T_star/(T_0*temp_ratio)  # T*/T
exit_area_crit_ratio = (expansion_ratio)*np.ones(4999)  # A_exit/A*
freestream_pres_crit_ratio = (P_star/P_amb)*np.ones(4999)  # P*/P_amb
v_exit_ratio = M*a_ratio*a_0/a_star  # v/v*

# --------------------------------------------------------------------------------




# ----- Nozzle Geometry Plotting -------------------------------------------------

nozzle_plot = np.sqrt( (A_star/area_ratio)*(4/np.pi) )/2
nozzle_length = ((d_exit - d_star)/2)/np.tan(np.deg2rad(half_angle))

# --------------------------------------------------------------------------------


print('')
print('Chamber Pressure: ' + str(round(P_0/1000, 0)) + ' kPa  (' + str(round(P_0/6894.76, 0)) + ' psi)')
print('Chamber Temperature: ', T_0, ' K')
print('Mass Flow Rate: ' +  str(round(m_dot*1000, 4)) + ' g/s')
print('Volumetric Flow Rate: ' + str(round(Q_scfm, 3)) + ' scfm')
print('')
print('Throat Diameter: ' + str(round(d_star*1000, 3)) + ' mm  (' + str(round(d_star/0.0254, 3)) + ' in)')
print('Throat Pressure: ' + str(round(P_star/1000, 2)) + ' kPa  (' + str(round(P_star/6894.76, 0)) + ' psi)')
print('Throat Temperature: ' + str(round(T_star, 1)) + ' K')
print('Throat Velocity: ' + str(round(a_star, 1)) + ' m/s')
print('Throat Reynold\'s Number (Diameter): ' + str(round(Re_star, 0)))
print('')
print('Exit Diameter: ' + str(round(d_exit*1000, 3)) + ' mm  (' + str(round(d_exit/0.0254, 3)) + ' in)')
print('Exit Temperature: ' + str(round(T_exit, 1)) + ' K')
print('Exit Pressure: ' + str(round(P_exit/1000, 3)) + ' kPa  (' + str(round(P_exit/6894.76, 3)) + ' psi)')
# print('Entrance Mach Number: ' + str(round(M_exit_sub, 3)))
print('Exit Mach Number: ' + str(round(M_exit_sup, 3)))
print('Exit Velocity: ' + str(round(v_exit, 1)) + ' m/s')
print('Exit Speed of Sound: ' + str(round(a_exit, 1)) + ' m/s')
print('')
if F > 1:
	print('Thrust: ', str(round(F, 3)) + ' N')
	print('Thrust from momentum exchange: ', str(round(F_mdotv, 3)) + ' N')
	print('Thrust from pressure difference: ', str(round(F_pdiff, 3)) + ' N')
else:
	print('Thrust: ', str(round(F*1000, 3)) + ' mN')
	print('Thrust from momentum exchange: ', str(round(F_mdotv*1000, 3)) + ' mN')
	print('Thrust from pressure difference: ', str(round(F_pdiff*1000, 3)) + ' mN')
print('')
print('Nozzle Length: ', str(round(nozzle_length*1000, 3)) + 'mm')


#if plot_relationships == True:

plt.figure(1)
plt.plot(M, area_ratio, '#1f77b4')
plt.plot(M, temp_ratio, '#ff7f0e')
plt.plot(M, pres_ratio, '#2ca02c')
plt.plot(M, a_ratio, '#ff7f0e', linestyle='--')
plt.plot(M, v_exit_ratio, '#c28285')
plt.plot(M, freestream_pres, '#2ca02c', linestyle='--')
plt.legend(['$A*/A$', '$T/T_0$', '$P/P_0$', '$a/a_0$', 'v/v_${exit}$'])
plt.xlabel('Mach')
#plt.ylim(0, 1)
plt.grid(True)

plt.title('Isentropic Relationships for ' + gas_type)
#plt.show(block=True)


plt.figure(2)
plt.subplot(211)
plt.plot(M, 1/area_ratio, '#1f77b4')
plt.plot(M, temp_crit_ratio, '#ff7f0e')
plt.plot(M, pres_crit_ratio, '#2ca02c')
plt.plot(M, exit_area_crit_ratio, '#1f77b4', linestyle='--')
plt.plot(M, freestream_pres_crit_ratio, '#2ca02c', linestyle='--')
plt.legend(['$A/A*$', '$T*/T$', '$P*/P$'])
plt.xlabel('Mach')
#plt.ylim(0, 1)
plt.grid(True)

plt.title('Critical Ratio Relationships for ' + gas_type)


#plt.subplot(212)
#lt.plot(M, nozzle_plot)
#plt.show(block=True)
    

#exit()