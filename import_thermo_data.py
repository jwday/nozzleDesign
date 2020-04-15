import pandas as pd
import numpy as np
import os

file = 'R134a_Sat_Props_tabDelim.txt'

df = pd.read_csv(file, delimiter='\t')
# Goal: Import P-T-h data for saturated substance, then use it to generate or plot data points on nozzle P-T diagrams
# Calculate enthalpy of gas flow through nozzle
# i.e. Supply enthalpy data (ex. exit enthalpy), determine corresponding saturation Temperature and Pressure, plot on graphs