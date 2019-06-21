# Nozzle maser control
from nozzle_code import nozzle
import pandas as pd
import numpy as np

## ---- OPTIONS ----------------------------------

P_0 = 114.7  # Total/Stagnation/Chamber Pressure (psi)
T_0 = 10  # Total/Stagnation/Chamber Temperature (C)
P_amb = 14.7  # Ambient Pressure (psi)

d_star = 0.3  # Throat diameter (mm) (1/64" = 0.397mm)
expansion_ratio = 1.1235  # 1.8048 for ideal expansion at 114.7 psi supply, 2.2447 for 164.7, 1.3225 for 0.2mm and 64.7 psi, 1.1235 for 0.3 and 44.7 psi
half_angle = 10  # Conical Nozzle expansion angle (degrees)

# Choices: R236fa, R134a, N2, CO2
gas_type = 'CO2'


## --------------------------------------------------------------------------------


list_of_P0 = list(np.linspace (44.7, 44.7, 1))
list_of_mdots = []
list_of_P_exit_subs = []
list_of_P_exit_supers = []
list_of_M_exit_sub = []
list_of_M_exit_supers = []
list_of_P_stars = []
list_of_thrusts = []


for i in range(len(list_of_P0)):
    P = list_of_P0[i]

    m_dot, P_star, T_star, a_star, Re_star, M_exit_sub, M_exit_sup, P_exit, T_exit, v_exit, a_exit, F, F_mdotv, F_pdiff = nozzle(P, T_0, P_amb, d_star, expansion_ratio, half_angle, gas_type)
    m_dot *= 1000
    P_exit_sub = P_exit[0]/6894.76
    P_exit_super = P_exit[1]/6894.76
    P_star /= 6894.76

    list_of_mdots.append(m_dot)
    list_of_P_exit_subs.append(P_exit_sub)
    list_of_P_exit_supers.append(P_exit_super)
    list_of_M_exit_sub.append(M_exit_sub)
    list_of_M_exit_supers.append(M_exit_sup)
    list_of_P_stars.append(P_star)
    list_of_thrusts.append(F)

# out =  [('Init Pressure (psi)', list_of_P0),
#         ('Mass Flow Rate (g/s)', list_of_mdots),
#         ('Exit Pressure (sub, psi)', list_of_P_exit_subs),
#         ('Exit Pressure (sup, psi)', list_of_P_exit_supers)]

out =  {'Init Pressure (psi)': list_of_P0,
        'Mass Flow Rate (g/s)': list_of_mdots,
        'Throat Pressure (psi)': list_of_P_stars,
        'Exit Pressure (psi)': list_of_P_exit_supers,
        'Thrust (N)': F}

df = pd.DataFrame.from_dict(out)
# df.columns = ['Mass Flow Rate (g/s)', 'Exit Pressure (sub, psi)', 'Exit Pressure (sup, psi)']
df.to_csv("nozzle_pressure.csv", index=False)


