# Nozzle master control
# Given a nozzle geometry and initial conditions, this code will sweep through a range of stagnation pressures and output the exit conditions
# It's interesting to note that the critical conditions are dependent ONLY on geometry and not stagnation conditions
from nozzle_code import nozzle
import pandas as pd
import numpy as np
import sys


## ---- OPTIONS --------------------------------------------------------------------

# Gas initial conditions
P_t_max = 64.7  # Max Total Pressure (psi)
P_amb = 14.7  # Ambient Pressure (psi)
T_t = 10  # Total Temperature (C)
vol = 20  # cu. cm

no_of_points = 30

# Nozzle geometry
d_star = 0.3  # Throat diameter (mm) (1/64" = 0.397mm)
expansion_ratio = 1.1235  # 1.8048 for ideal expansion at 114.7 psi supply, 2.2447 for 164.7, 1.3225 for 0.2mm and 64.7 psi, 1.1235 for 0.3 and 44.7 psi
half_angle = 10  # Conical Nozzle expansion angle (degrees)

# Gas Choices: R236fa, R134a, N2, CO2
gas_type = 'CO2'
k = 1.289
R = 8.314/0.04401



## ---- DO THE THING ----------------------------------------------------------------

list_of_P_ts = list(np.linspace (P_t_max, P_amb, no_of_points))

list_of_mdots = []
list_of_P_exits = []
list_of_v_exits = []
list_of_M_exits = []
list_of_thrusts = []
list_of_pressure_ratios = []
list_of_exit_shock_PRs = []

# for P_t in list_of_P_ts:
while m > 
    m_dot, M_crit_sub, M_crit_sup, PR_crit_sub, PR_crit_sup, PR_exit_shock, M_exit_behindshock, M_exit, P_exit, v_exit, F, F_mdotv, F_pdiff = nozzle(P_t, T_t, P_amb, d_star, expansion_ratio, half_angle, gas_type)
    list_of_mdots.append(m_dot*1000)
    list_of_P_exits.append(P_exit)
    list_of_v_exits.append(v_exit)
    list_of_M_exits.append(M_exit)
    list_of_thrusts.append(F)
    list_of_pressure_ratios.append(P_exit/P_t)

## ---- UPDATE THE THING ------------------------------------------------------------
# Search through the list of Pressure Ratios to find where the two critical conditions and the nozzle shock condition fit in
# Then at each of those Pressure Ratios, determine the Total Pressure and run it through nozzle() to get the exit conditions
# Then update the lists by inserting the exit conditions at the appropriate location

for i in range(len(list_of_pressure_ratios)):
    if ((1/PR_crit_sub) > list_of_pressure_ratios[i] and (1/PR_crit_sub) < list_of_pressure_ratios[i+1]) or ((1/PR_crit_sub) < list_of_pressure_ratios[i] and (1/PR_crit_sub) > list_of_pressure_ratios[i+1]):
        P_t = PR_crit_sub*P_amb
        m_dot, M_crit_sub, M_crit_sup, PR_crit_sub, PR_crit_sup, PR_exit_shock, M_exit_behindshock, M_exit, P_exit, v_exit, F, F_mdotv, F_pdiff = nozzle(P_t, T_t, P_amb, d_star, expansion_ratio, half_angle, gas_type)
        list_of_P_ts.insert(i+1, P_t)
        list_of_mdots.insert(i+1, m_dot*1000)
        list_of_P_exits.insert(i+1, P_exit)
        list_of_pressure_ratios.insert(i+1, P_exit/P_t)
        list_of_M_exits.insert(i+1, M_exit)
        list_of_v_exits.insert(i+1, v_exit)
        list_of_thrusts.insert(i+1, F) 
        print("Match for subsonic critical pressure ratio!")
        break
    else:
        pass

for i in range(len(list_of_pressure_ratios)):
    if (1/PR_exit_shock > list_of_pressure_ratios[i] and 1/PR_exit_shock < list_of_pressure_ratios[i+1]) or (1/PR_exit_shock < list_of_pressure_ratios[i] and 1/PR_exit_shock > list_of_pressure_ratios[i+1]):
        P_t = PR_exit_shock*P_amb
        m_dot, M_crit_sub, M_crit_sup, PR_crit_sub, PR_crit_sup, PR_exit_shock, M_exit_behindshock, M_exit, P_exit, v_exit, F, F_mdotv, F_pdiff = nozzle(P_t, T_t, P_amb, d_star, expansion_ratio, half_angle, gas_type)
        list_of_P_ts.insert(i+1, P_t)
        list_of_mdots.insert(i+1, m_dot*1000)
        list_of_P_exits.insert(i+1, P_exit)
        list_of_pressure_ratios.insert(i+1, P_exit/P_t)
        list_of_M_exits.insert(i+1, M_exit)
        list_of_v_exits.insert(i+1, v_exit)
        list_of_thrusts.insert(i+1, F) 
        print("Match for shock-at-exit pressure ratio!")
        break
    else:
        pass

for i in range(len(list_of_P_exits)):
    if (list_of_P_exits[i] > P_amb and list_of_P_exits[i+1] < P_amb) or (list_of_P_exits[i] < P_amb and list_of_P_exits[i+1] > P_amb):
        P_t = PR_crit_sup*P_amb
        m_dot, M_crit_sub, M_crit_sup, PR_crit_sub, PR_crit_sup, PR_exit_shock, M_exit_behindshock, M_exit, P_exit, v_exit, F, F_mdotv, F_pdiff = nozzle(P_t, T_t, P_amb, d_star, expansion_ratio, half_angle, gas_type)
        list_of_P_ts.insert(i+1, P_t)    
        list_of_mdots.insert(i+1, m_dot*1000)
        list_of_P_exits.insert(i+1, P_exit)
        list_of_pressure_ratios.insert(i+1, P_exit/P_t)
        list_of_M_exits.insert(i+1, M_exit)
        list_of_v_exits.insert(i+1, v_exit)
        list_of_thrusts.insert(i+1, F)
        print("Match for supersonic critical pressure ratio!")
        break
    else:
        pass





## ---- OUTPUT THE STUFF -----------------------------------------------------------

list_of_NaNs = [np.NaN for i in range(len(list_of_P_ts)-1)]

list_of_M_crit_subs = [M_crit_sub]
list_of_M_crit_sups = [M_crit_sup]
list_of_PR_crit_subs = [1/PR_crit_sub]
list_of_PR_crit_sups = [1/PR_crit_sup]
list_of_PR_exit_shocks = [1/PR_exit_shock]
list_of_M_exit_behindshock = [M_exit_behindshock]

for list in [list_of_M_crit_subs, list_of_M_crit_sups, list_of_PR_crit_subs, list_of_PR_crit_sups, list_of_PR_exit_shocks, list_of_M_exit_behindshock]:
        list.extend(list_of_NaNs)

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
        'Critical PR, Sup': list_of_PR_crit_subs,
        '': [np.NaN for i in range(len(list_of_P_ts))],
        'PR for Shock at Exit': list_of_PR_exit_shocks,
        'Mach Behind Exit Shock': list_of_M_exit_behindshock}

df = pd.DataFrame.from_dict(out)
df.to_csv("nozzle_pressure.csv", index=False)