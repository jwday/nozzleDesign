# plot_mach_area_relation.py
# Simple exercise to plot the Mach Number-Area relation for three gasses (CO2, R134a, Air)

import numpy as np
import sys
import math
import matplotlib.pyplot as plt
import sympy as sm
import scipy.optimize as opti
import pandas as pd
import seaborn as sns

k_gas_names = {'R236fa': 1.083, 'Air': 1.401, 'CO2': 1.289}

def objective(M):
	X = M**2
	S = ((1 + L*X)**Q) / (X*(P**Q))
	f = np.sqrt(S)
	return f

machs = np.arange(0.01, 4.01, 0.01)
mach_area_relations = pd.DataFrame(columns=['Mach'] + list(k_gas_names.keys()))
mach_area_relations['Mach'] = machs

for gas_name in list(k_gas_names.keys()):
	k = k_gas_names[gas_name]
	P = (k+1)/2
	L = (k-1)/2
	W = k/(k-1)
	Q = P/L

	for i, M in enumerate(machs):
		mach_area_relations[gas_name][i] = objective(M)

fig, ax = plt.subplots(figsize=(6, 3), dpi=150)
sns.set_style('whitegrid')
mach_area_relations.plot(ax=ax, x='Mach')
# for gas_name in list(k_gas_names.keys()):
# sns.lineplot(x='Mach',
# 				y=list(k_gas_names.keys()),
# 				data=mach_area_relations,
			#  hue='Trial', 
			#  style=gas_name,  # Show each trial individually instead of an aggregate
			#  estimator=np.mean,  # Show each trial individually instead of an aggregate
			#  marker=data_marker
				# )
ax.grid(which='major', axis='both', linestyle='--')

ax.set_xlim(0, 3.5)
ax.set_ylim(0, 10)
ax.set_xlabel("Mach Number, $M$")
ax.set_ylabel(r"Area Ratio, $\frac{A}{A^{*}}$")
plt.title("Mach Number - Area Relation")
plt.subplots_adjust(bottom=0.15)
plt.subplots_adjust(top=0.925)
plt.legend(loc='best')

plt.show()