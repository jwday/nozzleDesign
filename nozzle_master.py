# Nozzle maser control
from nozzle_code import nozzle
import pandas as pd
import numpy as np

## ---- OPTIONS ----------------------------------

# P_t = 44.7  # Total/Stagnation/Chamber Pressure (psi)
T_t = 10  # Total/Stagnation/Chamber Temperature (C)
P_amb = 14.7  # Ambient Pressure (psi)

d_star = 0.3  # Throat diameter (mm) (1/64" = 0.397mm)
expansion_ratio = 1.1235  # 1.8048 for ideal expansion at 114.7 psi supply, 2.2447 for 164.7, 1.3225 for 0.2mm and 64.7 psi, 1.1235 for 0.3 and 44.7 psi
half_angle = 10  # Conical Nozzle expansion angle (degrees)

# Choices: R236fa, R134a, N2, CO2
gas_type = 'CO2'
k = 1.289


## --------------------------------------------------------------------------------


list_of_P_ts = list(np.linspace (14.7, 64.7, 30))
list_of_mdots = []
list_of_P_exits = []
list_of_v_exits = []
list_of_M_exits = []
list_of_M_crit_subs = []
list_of_M_crit_sups = []
list_of_PR_crit_subs = []
list_of_PR_crit_sups = []
list_of_thrusts = []
list_of_pressure_ratios = []
list_of_exit_shock_PRs = []
list_of_M_exit_behindshock = []

for P_t in list_of_P_ts:
    m_dot, M_crit_sub, M_crit_sup, PR_crit_sub, PR_crit_sup, PR_exit_shock, M_exit_behindshock, M_exit, P_exit, v_exit, F, F_mdotv, F_pdiff = nozzle(P_t, T_t, P_amb, d_star, expansion_ratio, half_angle, gas_type)
    
    m_dot *= 1000
    P_exit /= 6894.76

    list_of_mdots.append(m_dot)
    list_of_P_exits.append(P_exit)
    list_of_v_exits.append(v_exit)
    list_of_M_exits.append(M_exit)
    list_of_thrusts.append(F)
    list_of_pressure_ratios.append(P_exit/P_t)
    list_of_M_crit_subs.append(M_crit_sub)
    list_of_M_crit_sups.append(M_crit_sup)
    list_of_PR_crit_subs.append(1/PR_crit_sub)
    list_of_PR_crit_sups.append(1/PR_crit_sup)
    list_of_exit_shock_PRs.append(1/PR_exit_shock)
    list_of_M_exit_behindshock.append(M_exit_behindshock)

# out =  [('Init Pressure (psi)', list_of_P_ts),
#         ('Mass Flow Rate (g/s)', list_of_mdots),
#         ('Exit Pressure (sub, psi)', list_of_P_exit_subs),
#         ('Exit Pressure (sup, psi)', list_of_P_exit_supers)]

out =  {'Init Pressure (psi)': list_of_P_ts,
        'Mass Flow Rate (g/s)': list_of_mdots,
        'Exit Pressure (psi)': list_of_P_exits,
        'Pressure Ratio': list_of_pressure_ratios,
        'Exit Mach Number': list_of_M_exits,
        'Exit Velocity (m/s)': list_of_v_exits,
        'Thrust (N)': list_of_thrusts,
        '': [np.NaN for i in range(len(list_of_P_ts))],
        'Critical Mach, Sub': list_of_M_crit_subs,
        'Critical Mach, Sup': list_of_M_crit_sups,
        'Critical PR, Sub': list_of_PR_crit_subs,
        'Critical PR, Sup': list_of_PR_crit_sups,
        '': [np.NaN for i in range(len(list_of_P_ts))],
        'PR for Shock at Exit': list_of_exit_shock_PRs,
        'Mach Behind Exit Shock': list_of_M_exit_behindshock}

df = pd.DataFrame.from_dict(out)
# df.columns = ['Mass Flow Rate (g/s)', 'Exit Pressure (sub, psi)', 'Exit Pressure (sup, psi)']
df.to_csv("nozzle_pressure.csv", index=False)


