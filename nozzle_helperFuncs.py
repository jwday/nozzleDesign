# nozzle_helperFuncs.py
# Helper functions to interpolate/extrapolate real gas data from NIST. Used only when not using isentropic model in simulation, or when heat transfer is being estimated.

from scipy import stats
import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
from scipy.interpolate import interp2d
import scipy.optimize as opti
import itertools


## ==================================================================================
## ---- HELPER FUNCTIONS ------------------------------------------------------------
## ==================================================================================

# --------------------------------------------------------------------------------
# Create a 2D function to return viscosity of a fluid at a given (P,T) based on NIST data and using linear interpolation
def create_visc_func_gas(fluid_props):
	visc_data = pd.DataFrame(columns=['Temperature (K)'])
	pressures = []		
	for pres in list(fluid_props.keys()):
		f = fluid_props[pres][['Temperature (K)', 'Viscosity (uPa*s)', 'Phase']].rename(columns={'Viscosity (uPa*s)': pres})	# Pull out Temperature, Viscosity, and Phase data for each Pressure, and rename column to numercial pressure (in MPa)
		f = f.drop(f.index[f['Phase'] != 'vapor']).drop('Phase', axis=1) 														# Drop all rows that are not pertaining to vapor, then drop the 'Phase' column entirely so we're just left with Viscosity
		if not f.empty:
			visc_data = pd.merge(visc_data, f, how='outer', on='Temperature (K)')												# Merge on 'Temperature (K)'. Non-overlapping values will have NaNs.
			visc_data = visc_data.sort_values('Temperature (K)').reset_index(drop=True)											# Sort by temperature and reset index
			visc_data.interpolate(inplace=True)
			pressures.append(fluid_props[pres]['Pressure (MPa)'][0]*1E6)														# Interpolate across the non-overlapping values to fill in any NaNs that lie BETWEEN data points.
		else:
			pass
	visc_data = visc_data.fillna(method='bfill', axis=1)																		# We can't interpolate over NaN's, and it throws a NaN for all the pressures where there is no viscosity data (or when a phase change occurs), so...just backfill them for the time being
	visc_data = visc_data.dropna(axis=0, how='any', thresh=2)																	# Drop any rows that are completely NaN'ed out
	visc_data = visc_data.dropna(axis=1, how='all')																				# Drop the remaining all-Nan columns
	visc_func = interp2d(pressures[:len(visc_data.columns)-1], visc_data['Temperature (K)'].values, visc_data.iloc[:,1:].values)
	return visc_func



# --------------------------------------------------------------------------------
# Create two 1D functions to return the P_fg given T, and another to return the T_fg given P. Also determine the phase boundaries and return an ax object for future plotting
def create_phase_funcs(fluid_props, P_trip, T_trip):
	saturated_temp = [T_trip]
	saturated_pres = [P_trip]
	for pres in list(fluid_props.keys()):
		fluid_props[pres]['Phase Change'] = fluid_props[pres]['Phase'].ne(fluid_props[pres]['Phase'].shift().bfill()).astype(bool)		# Add a column called 'Phase Change' to identify when a phase change takes place
		if fluid_props[pres]['Phase Change'].any():																						# If any of the values in ['Phase Change'] are True:
			saturated_temp.append(fluid_props[pres][fluid_props[pres]['Phase Change'] == True]['Temperature (K)'].values[0])				# Returns the temperature of phase change at specified pressure
			saturated_pres.append(fluid_props[pres]['Pressure (MPa)'][0]*1E6)																# Returns the specified pressure (in Pa)
		else:
			# saturated_temp.append(np.nan)
			# saturated_pres.append(np.nan)
			pass
	log_saturated_pres = [np.log(x) for x in saturated_pres]																			# Behavior can be modeled and interpolated logarithmically (T vs log(P) is very linear)
	f = interp1d(pd.DataFrame(saturated_temp).dropna()[0].values, pd.DataFrame(log_saturated_pres).dropna()[0].values, kind='linear', fill_value='extrapolate')		# Linearlly interpolate over T vs log(P), return a function f. Need to convert to df to drop na's to use quadratic fit
	g = interp1d(pd.DataFrame(log_saturated_pres).dropna()[0].values, pd.DataFrame(saturated_temp).dropna()[0].values, kind='linear', fill_value='extrapolate')		# Linearlly interpolate over T vs log(P), return a function g. Need to convert to df to drop na's to use quadratic fit

	def saturated_pres_from_temp(temp):
		h = np.e**f(temp)																												# When a temp is specified, use the function f to return a log(P), then exponentiate it to return P
		return h

	def saturated_temp_from_pres(pres):
		h = float(g(np.log(pres)))																										# When a temp is specified, supply the function g with a log(P) to return a T. g returns a 0D array -- convert to float
		return h

	# --------------------------------------------------------------------------------
	# Function to set up phase data to plot, returns DataFrame containing all phase data
	# OR USE CLAUSIUS-CLAPEYRON EQUATION? Need latent heat of sublimation for R134a, don't have it, can't find it.
	phase_data = pd.DataFrame(columns=['Temperature (K)', 'Pressure (Pa)', 'Phase', 'Dataset'])
	phase_data = pd.concat([phase_data, pd.DataFrame({'Temperature (K)': saturated_temp, 'Pressure (Pa)':saturated_pres, 'Phase':'', 'Dataset':'NIST Data'})], ignore_index=True)

	# Liquid-Vapor Phase Change Line
	xnew_lv = np.arange(T_trip, max(saturated_temp)+40, 1)
	ynew_lv = [saturated_pres_from_temp(x) for x in xnew_lv]
	# ynew_lv_cc = [np.exp(-P_trip*(L/R)*((1/xnew_lv) - (1/T_trip)))]	# Clausius Clapeyron
	phase_data = pd.concat([phase_data, pd.DataFrame({'Temperature (K)': xnew_lv, 'Pressure (Pa)':ynew_lv, 'Phase':'Liquid-Vapor Coexistence', 'Dataset':'Extrapolated'})], ignore_index=True)

	# Solid-Vapor Phase Change Line
	xnew_sv = np.arange(min(saturated_temp)-169, T_trip, 1)
	ynew_sv = [saturated_pres_from_temp(x) for x in xnew_sv]
	phase_data = pd.concat([phase_data, pd.DataFrame({'Temperature (K)': xnew_sv, 'Pressure (Pa)':ynew_sv, 'Phase':'Solid-Vapor Coexistence', 'Dataset':'Extrapolated'})], ignore_index=True)

	# Solid-Liquid Phase Change Line (just a guess)
	xnew_sl = [T_trip, T_trip+.01]
	ynew_sl = [P_trip, P_trip + 1000000]
	phase_data = pd.concat([phase_data, pd.DataFrame({'Temperature (K)': xnew_sl, 'Pressure (Pa)':ynew_sl, 'Phase':'Solid-Liquid Coexistence', 'Dataset':'Extrapolated'})], ignore_index=True)
	phase_data = phase_data.sort_values('Temperature (K)').reset_index(drop=True)
	
	return saturated_pres_from_temp, saturated_temp_from_pres, phase_data



# --------------------------------------------------------------------------------
# Create a 2D function to return specific heat, Cp, using interpolation of NIST data
def create_cp_func(fluid_props):
	cp_data = pd.DataFrame(columns=['Temperature (K)'])
	pressures = []		
	for pres in list(fluid_props.keys()):
		f = fluid_props[pres][['Temperature (K)', 'Cp (J/g*K)', 'Phase']].rename(columns={'Cp (J/g*K)': pres})					# Pull out Temperature, Cp, and Phase data for each Pressure, and rename column to numercial pressure (in MPa)
		f = f.drop(f.index[f['Phase'] != 'vapor']).drop('Phase', axis=1) 														# Drop all rows that are not pertaining to vapor, then drop the 'Phase' column entirely so we're just left with Viscosity
		if not f.empty:
			cp_data = pd.merge(cp_data, f, how='outer', on='Temperature (K)')													# Merge on 'Temperature (K)'. Non-overlapping values will have NaNs.
			cp_data = cp_data.sort_values('Temperature (K)').reset_index(drop=True)													# Sort by temperature and reset index
			cp_data.interpolate(inplace=True)
			pressures.append(fluid_props[pres]['Pressure (MPa)'][0]*1E6)														# Interpolate across the non-overlapping values to fill in any NaNs that lie BETWEEN data points.
		else:
			pass
	cp_data = cp_data.fillna(method='bfill', axis=1)																			# We can't interpolate over NaN's, and it throws a NaN for all the pressures where there is no viscosity data (or when a phase change occurs), so...just backfill them for the time being
	cp_data = cp_data.dropna(axis=0, how='any', thresh=2)																		# Drop any rows that are completely NaN'ed out
	cp_data = cp_data.dropna(axis=1, how='all')																					# Drop the remaining all-Nan columns
	cp_func = interp2d(pressures[:len(cp_data.columns)-1], cp_data['Temperature (K)'].values, cp_data.iloc[:,1:].values)		# Cp, J/g*K
	return cp_func



# --------------------------------------------------------------------------------
# Create a 2D function to return specific heat, Cv, using interpolation of NIST data
def create_cv_func(fluid_props):
	cv_data = pd.DataFrame(columns=['Temperature (K)'])
	pressures = []		
	for pres in list(fluid_props.keys()):
		f = fluid_props[pres][['Temperature (K)', 'Cv (J/g*K)', 'Phase']].rename(columns={'Cp (J/g*K)': pres})					# Pull out Temperature, Cp, and Phase data for each Pressure, and rename column to numercial pressure (in MPa)
		f = f.drop(f.index[f['Phase'] != 'vapor']).drop('Phase', axis=1) 														# Drop all rows that are not pertaining to vapor, then drop the 'Phase' column entirely so we're just left with Viscosity
		if not f.empty:
			cv_data = pd.merge(cv_data, f, how='outer', on='Temperature (K)')													# Merge on 'Temperature (K)'. Non-overlapping values will have NaNs.
			cv_data = cv_data.sort_values('Temperature (K)').reset_index(drop=True)													# Sort by temperature and reset index
			cv_data.interpolate(inplace=True)
			pressures.append(fluid_props[pres]['Pressure (MPa)'][0]*1E6)														# Interpolate across the non-overlapping values to fill in any NaNs that lie BETWEEN data points.
		else:
			pass
	cv_data = cv_data.fillna(method='bfill', axis=1)																			# We can't interpolate over NaN's, and it throws a NaN for all the pressures where there is no viscosity data (or when a phase change occurs), so...just backfill them for the time being
	cv_data = cv_data.dropna(axis=0, how='any', thresh=2)																		# Drop any rows that are completely NaN'ed out
	cv_data = cv_data.dropna(axis=1, how='all')																					# Drop the remaining all-Nan columns
	cv_func = interp2d(pressures[:len(cv_data.columns)-1], cv_data['Temperature (K)'].values, cv_data.iloc[:,1:].values)		# Cv, J/g*K
	return cv_func



# --------------------------------------------------------------------------------
# Create a 2D function to return thermal conductivity, ktc, using interpolation of NIST data
def create_ktc_func(fluid_props):
	ktc_data = pd.DataFrame(columns=['Temperature (K)'])
	pressures = []		
	for pres in list(fluid_props.keys()):
		f = fluid_props[pres][['Temperature (K)', 'Therm. Cond. (W/m*K)', 'Phase']].rename(columns={'Therm. Cond. (W/m*K)': pres})					# Pull out Temperature, k, and Phase data for each Pressure, and rename column to numercial pressure (in MPa)
		f = f.drop(f.index[f['Phase'] != 'vapor']).drop('Phase', axis=1) 														# Drop all rows that are not pertaining to vapor, then drop the 'Phase' column entirely so we're just left with Viscosity
		if not f.empty:
			ktc_data = pd.merge(ktc_data, f, how='outer', on='Temperature (K)')													# Merge on 'Temperature (K)'. Non-overlapping values will have NaNs.
			ktc_data = ktc_data.sort_values('Temperature (K)').reset_index(drop=True)											# Sort by temperature and reset index
			ktc_data.interpolate(inplace=True)
			pressures.append(fluid_props[pres]['Pressure (MPa)'][0]*1E6)														# Interpolate across the non-overlapping values to fill in any NaNs that lie BETWEEN data points.
		else:
			pass
	ktc_data = ktc_data.fillna(method='bfill', axis=1)																			# We can't interpolate over NaN's, and it throws a NaN for all the pressures where there is no viscosity data (or when a phase change occurs), so...just backfill them for the time being
	ktc_data = ktc_data.dropna(axis=0, how='any', thresh=2)																		# Drop any rows that are completely NaN'ed out
	ktc_data = ktc_data.dropna(axis=1, how='all')																				# Drop the remaining all-Nan columns
	ktc_func = interp2d(pressures[:len(ktc_data.columns)-1], ktc_data['Temperature (K)'].values, ktc_data.iloc[:,1:].values)	# k, W/m-K
	return ktc_func




# --------------------------------------------------------------------------------
# Create a NEW 2D function to return enthalpy of a fluid at a given (P,T) based on NIST data and using linear interpolation
# def create_h_from_PT_gas_func(fluid_props):
# 	h_sp_data = pd.DataFrame(columns=['Temperature (K)'])
# 	pressures = []		
# 	for pres in list(fluid_props.keys()):
# 		f = fluid_props[pres][['Temperature (K)', 'Enthalpy (kJ/kg)', 'Phase']].rename(columns={'Enthalpy (kJ/kg)': pres})		# Pull out Temperature, Enthalpy, and Phase data for each Pressure, and rename column to numercial pressure (in MPa)
# 		f = f.drop(f.index[f['Phase'] != 'vapor']).drop('Phase', axis=1).dropna() 												# Drop all rows that are not pertaining to vapor, then drop the 'Phase' column entirely so we're just left with enthalpy
# 		if not f.empty:
# 			h_sp_data = pd.merge(h_sp_data, f, how='outer', on='Temperature (K)')												# Merge on 'Temperature (K)'. Non-overlapping values will have NaNs.
# 			h_sp_data = h_sp_data.sort_values('Temperature (K)').reset_index(drop=True)											# Sort by volume and reset index
# 			h_sp_data = h_sp_data.set_index('Temperature (K)').interpolate(method='index', limit_area='inside').reset_index()	# Fill in any missing values that fall between tabulated values
# 			pressures.append(fluid_props[pres]['Pressure (MPa)'][f.index[0]]*1E6)												# Interpolate across the non-overlapping values to fill in any NaNs that lie BETWEEN data points.
# 		else:
# 			pass

# 	start = 150
# 	stop = h_sp_data['Temperature (K)'].iloc[-1]
# 	step = h_sp_data['Temperature (K)'].diff().median()
# 	extrap_range = np.arange(start, stop, step)																					# Temperature range over which to extrapolate (pressure range is already okay, no need to extrapolate below 100 KPa)

# 	h_sp_data_new = pd.DataFrame(columns=['Temperature (K)'])
# 	for pres in h_sp_data.columns[1:]:
# 		h_sp_data_slice = h_sp_data[['Temperature (K)', pres]].dropna()															# Get the values over the range you wish to use as the basis for your interpolation (and drop NaNs)
# 		if not h_sp_data_slice.empty:
# 			f = interp1d(h_sp_data_slice['Temperature (K)'], h_sp_data_slice[pres].apply(np.log), fill_value='extrapolate')			# Make your interpolation function, first applying np.log() to all the values to make the extrapolation function from
# 			h_sp_data_fill = pd.DataFrame({'Temperature (K)': extrap_range, pres:np.exp(f(extrap_range))})							# Extrapolate over the range using said interp1d function
# 			h_sp_data_new = h_sp_data_new.merge(h_sp_data_fill, how='outer', on='Temperature (K)')
# 	h_sp_from_PT = interp2d(pressures[:len(h_sp_data_new.columns)-1], h_sp_data_new['Temperature (K)'].values, h_sp_data_new.iloc[:,1:].values)
# 	return h_sp_from_PT




# --------------------------------------------------------------------------------
# Create an EVEN NEWER 2D function to return enthalpy of a fluid at a given (P,T) based on NIST data and using linear interpolation
def create_h_from_PT_gas_func(fluid_props):
	h_sp_data = pd.DataFrame(columns=['Temperature (K)'])
	pressures = []		
	for pres in list(fluid_props.keys()):
		f = fluid_props[pres][['Temperature (K)', 'Enthalpy (kJ/kg)', 'Phase']].rename(columns={'Enthalpy (kJ/kg)': pres})		# Pull out Temperature, Internal Energy, and Phase data for each Pressure, and rename column to numercial pressure (in MPa)
		f = f.drop(f.index[f['Phase'] != 'vapor']).drop('Phase', axis=1).dropna() 												# Drop all rows that are not pertaining to vapor, then drop the 'Phase' column entirely so we're just left with internal energy
		if not f.empty:
			h_sp_data = pd.merge(h_sp_data, f, how='outer', on='Temperature (K)')												# Merge on 'Temperature (K)'. Non-overlapping values will have NaNs.
			h_sp_data = h_sp_data.sort_values('Temperature (K)').reset_index(drop=True)											# Sort by volume and reset index
			h_sp_data = h_sp_data.set_index('Temperature (K)').interpolate(method='index', limit_area='inside').reset_index()	# Fill in any missing values that fall between tabulated values
			pressures.append(fluid_props[pres]['Pressure (MPa)'][f.index[0]]*1E6)												# Interpolate across the non-overlapping values to fill in any NaNs that lie BETWEEN data points.
		else:
			pass

	start = 150
	stop = h_sp_data['Temperature (K)'].iloc[-1]
	step = h_sp_data['Temperature (K)'].diff().median()
	extrap_range = np.arange(start, stop, step)																					# Temperature range over which to extrapolate (pressure range is already okay, no need to extrapolate below 100 KPa)
	h_sp_data_new = pd.DataFrame(columns=['Temperature (K)'])
	h_sp_extrap = pd.DataFrame(columns=['Temperature (K)'])

	for pres in h_sp_data.columns[1:]:
		h_sp_data_slice = h_sp_data[['Temperature (K)', pres]].dropna()															# Get the values over the range you wish to use as the basis for your interpolation (and drop NaNs)
		if not h_sp_data_slice.empty:
			slope, intercept, r_value, p_value, std_err = stats.linregress(h_sp_data_slice['Temperature (K)'].values, h_sp_data_slice[pres].apply(np.log).values)
			def f(x):																											# Make your interpolation function, first applying np.log() to all the values to make the extrapolation function from
				return (slope*x + intercept)

			h_sp_data_fill = pd.DataFrame({'Temperature (K)': extrap_range, pres: np.exp(f(extrap_range))})									# Extrapolate over the range using said interp1d function
			h_sp_data_new = h_sp_data_new.merge(h_sp_data_fill, how='outer', on='Temperature (K)')
	h_sp_from_PT = interp2d(pressures[:len(h_sp_data_new.columns)-1], h_sp_data_new['Temperature (K)'].values, h_sp_data_new.iloc[:,1:].values)

	return h_sp_from_PT




# --------------------------------------------------------------------------------
# Create a NEW 2D function to return internal energy of a fluid at a given (P,T) based on NIST data and using linear interpolation
# Oh snap this might only be valid for single-phase substances
# In fact, that must be true because a phase change occurs at a fixed P-T for varying u
def create_u_from_PT_gas_func(fluid_props):
	u_sp_data = pd.DataFrame(columns=['Temperature (K)'])
	pressures = []		
	for pres in list(fluid_props.keys()):
		f = fluid_props[pres][['Temperature (K)', 'Internal Energy (kJ/kg)', 'Phase']].rename(columns={'Internal Energy (kJ/kg)': pres})		# Pull out Temperature, Internal Energy, and Phase data for each Pressure, and rename column to numercial pressure (in MPa)
		f = f.drop(f.index[f['Phase'] != 'vapor']).drop('Phase', axis=1).dropna() 												# Drop all rows that are not pertaining to vapor, then drop the 'Phase' column entirely so we're just left with internal energy
		if not f.empty:
			u_sp_data = pd.merge(u_sp_data, f, how='outer', on='Temperature (K)')												# Merge on 'Temperature (K)'. Non-overlapping values will have NaNs.
			u_sp_data = u_sp_data.sort_values('Temperature (K)').reset_index(drop=True)											# Sort by volume and reset index
			u_sp_data = u_sp_data.set_index('Temperature (K)').interpolate(method='index', limit_area='inside').reset_index()	# Fill in any missing values that fall between tabulated values
			pressures.append(fluid_props[pres]['Pressure (MPa)'][f.index[0]]*1E6)												# Interpolate across the non-overlapping values to fill in any NaNs that lie BETWEEN data points.
		else:
			pass

	start = 0
	stop = u_sp_data['Temperature (K)'].iloc[-1]
	step = u_sp_data['Temperature (K)'].diff().median()
	extrap_range = np.arange(start, stop, step)																					# Temperature range over which to extrapolate (pressure range is already okay, no need to extrapolate below 100 KPa)
	u_sp_data_new = pd.DataFrame(columns=['Temperature (K)'])
	u_sp_extrap = pd.DataFrame(columns=['Temperature (K)'])

	for pres in u_sp_data.columns[1:]:
		u_sp_data_slice = u_sp_data[['Temperature (K)', pres]].dropna()															# Get the values over the range you wish to use as the basis for your interpolation (and drop NaNs)
		if not u_sp_data_slice.empty:
			slope, intercept, r_value, p_value, std_err = stats.linregress(u_sp_data_slice['Temperature (K)'].values, u_sp_data_slice[pres].apply(np.log).values)
			def f(x):																											# Make your interpolation function, first applying np.log() to all the values to make the extrapolation function from
				return (slope*x + intercept)

			u_sp_data_fill = pd.DataFrame({'Temperature (K)': extrap_range, pres: np.exp(f(extrap_range))})									# Extrapolate over the range using said interp1d function
			u_sp_data_new = u_sp_data_new.merge(u_sp_data_fill, how='outer', on='Temperature (K)')
	u_sp_from_PT = interp2d(pressures[:len(u_sp_data_new.columns)-1], u_sp_data_new['Temperature (K)'].values, u_sp_data_new.iloc[:,1:].values)

	return u_sp_from_PT




# --------------------------------------------------------------------------------
# Create a 2D function to return internal energy of a fluid at a given (P,T) based on NIST data and using linear interpolation
# def create_u_from_PT_gas_func(fluid_props):
# 	u_sp_data = pd.DataFrame(columns=['Temperature (K)'])
# 	pressures = []		
# 	for pres in list(fluid_props.keys()):
# 		f = fluid_props[pres][['Temperature (K)', 'Internal Energy (kJ/kg)', 'Phase']].rename(columns={'Internal Energy (kJ/kg)': pres})		# Pull out Temperature, Enthalpy, and Phase data for each Pressure, and rename column to numercial pressure (in MPa)
# 		f = f.drop(f.index[f['Phase'] != 'vapor']).drop('Phase', axis=1) 														# Drop all rows that are not pertaining to vapor, then drop the 'Phase' column entirely so we're just left with Viscosity
# 		if not f.empty:
# 			u_sp_data = pd.merge(u_sp_data, f, how='outer', on='Temperature (K)')												# Merge on 'Temperature (K)'. Non-overlapping values will have NaNs.
# 			u_sp_data = u_sp_data.sort_values('Temperature (K)').reset_index(drop=True)											# Sort by temperature and reset index
# 			u_sp_data.interpolate(inplace=True)
# 			pressures.append(fluid_props[pres]['Pressure (MPa)'][0]*1E6)														# Interpolate across the non-overlapping values to fill in any NaNs that lie BETWEEN data points.
# 		else:
# 			pass
# 	u_sp_data = u_sp_data.fillna(method='bfill', axis=1)																		# We can't interpolate over NaN's, and it throws a NaN for all the pressures where there is no viscosity data (or when a phase change occurs), so...just backfill them for the time being
# 	u_sp_data = u_sp_data.dropna(axis=0, how='any', thresh=2)																	# Drop any rows that are completely NaN'ed out
# 	u_sp_data = u_sp_data.dropna(axis=1, how='all')																				# Drop the remaining all-Nan columns
# 	u_from_PT_gas_func = interp2d(pressures[:len(u_sp_data.columns)-1], u_sp_data['Temperature (K)'].values, u_sp_data.iloc[:,1:].values)
# 	return u_from_PT_gas_func




# --------------------------------------------------------------------------------
# Create a NEW 2D function to return density of a fluid at a given (P,T) based on NIST data and using linear interpolation
def create_r_from_PT_gas_func(fluid_props):
	r_data = pd.DataFrame(columns=['Temperature (K)'])
	pressures = []		
	for pres in list(fluid_props.keys()):
		f = fluid_props[pres][['Temperature (K)', 'Density (g/ml)', 'Phase']].rename(columns={'Density (g/ml)': pres})		# Pull out Temperature, Internal Energy, and Phase data for each Pressure, and rename column to numercial pressure (in MPa)
		f = f.drop(f.index[f['Phase'] != 'vapor']).drop('Phase', axis=1).dropna() 												# Drop all rows that are not pertaining to vapor, then drop the 'Phase' column entirely so we're just left with enthalpy
		if not f.empty:
			r_data = pd.merge(r_data, f, how='outer', on='Temperature (K)')												# Merge on 'Temperature (K)'. Non-overlapping values will have NaNs.
			r_data = r_data.sort_values('Temperature (K)').reset_index(drop=True)											# Sort by volume and reset index
			r_data = r_data.set_index('Temperature (K)').interpolate(method='index', limit_area='inside').reset_index()	# Fill in any missing values that fall between tabulated values
			pressures.append(fluid_props[pres]['Pressure (MPa)'][f.index[0]]*1E6)												# Interpolate across the non-overlapping values to fill in any NaNs that lie BETWEEN data points.
		else:
			pass

	start = 100
	stop = r_data['Temperature (K)'].iloc[-1]
	step = r_data['Temperature (K)'].diff().median()
	extrap_range = np.arange(start, stop, step)																					# Temperature range over which to extrapolate (pressure range is already okay, no need to extrapolate below 100 KPa)

	r_data_new = pd.DataFrame(columns=['Temperature (K)'])
	for pres in r_data.columns[1:]:
		r_data_slice = r_data[['Temperature (K)', pres]].dropna()															# Get the values over the range you wish to use as the basis for your interpolation (and drop NaNs)
		if not r_data_slice.empty:
			# f = interp1d(r_data_slice['Temperature (K)'], r_data_slice[pres].apply(np.log), kind='linear', fill_value='extrapolate')							# Make your interpolation function, first applying np.log() to all the values to make the extrapolation function from
			# slope, intercept, r_value, p_value, std_err = stats.linregress(r_data_slice['Temperature (K)'].values, r_data_slice[pres].apply(np.log).values)
			p = np.polyfit(r_data_slice['Temperature (K)'].apply(np.log).values, r_data_slice[pres].apply(np.log).values, 2)
			def f(x):
				x = np.log(x)																											# Make your interpolation function, first applying np.log() to all the values to make the extrapolation function from
				return np.exp(p[0]*x**2 + p[1]*x + p[2])

			r_data_fill = pd.DataFrame({'Temperature (K)': extrap_range, pres: f(extrap_range)})							# Extrapolate over the range using said interp1d function
			r_data_new = r_data_new.merge(r_data_fill, how='outer', on='Temperature (K)')
	r_from_PT = interp2d(pressures[:len(r_data_new.columns)-1], r_data_new['Temperature (K)'].values, r_data_new.iloc[:,1:].values)
	return r_from_PT




# --------------------------------------------------------------------------------
# Create a 2D function to return REAL density of a fluid at a given (P,T) based on NIST data and using linear interpolation
# def create_r_from_PT_gas_func(fluid_props):
# 	r_data = pd.DataFrame(columns=['Temperature (K)'])
# 	pressures = []		
# 	for pres in list(fluid_props.keys()):
# 		f = fluid_props[pres][['Temperature (K)', 'Density (g/ml)', 'Phase']].rename(columns={'Density (g/ml)': pres})			# Pull out Temperature, Enthalpy, and Phase data for each Pressure, and rename column to numercial pressure (in MPa)
# 		f = f.drop(f.index[f['Phase'] != 'vapor']).drop('Phase', axis=1) 														# Drop all rows that are not pertaining to vapor, then drop the 'Phase' column entirely so we're just left with Viscosity
# 		if not f.empty:
# 			r_data = pd.merge(r_data, f, how='outer', on='Temperature (K)')												# Merge on 'Temperature (K)'. Non-overlapping values will have NaNs.
# 			r_data = r_data.sort_values('Temperature (K)').reset_index(drop=True)											# Sort by temperature and reset index
# 			r_data.interpolate(inplace=True)
# 			pressures.append(fluid_props[pres]['Pressure (MPa)'][0]*1E6)														# Interpolate across the non-overlapping values to fill in any NaNs that lie BETWEEN data points.
# 		else:
# 			pass
# 	r_data = r_data.fillna(method='bfill', axis=1)																		# We can't interpolate over NaN's, and it throws a NaN for all the pressures where there is no viscosity data (or when a phase change occurs), so...just backfill them for the time being
# 	r_data = r_data.dropna(axis=0, how='any', thresh=2)																	# Drop any rows that are completely NaN'ed out
# 	r_data = r_data.dropna(axis=1, how='all')																				# Drop the remaining all-Nan columns
# 	r_from_PT_gas_func = interp2d(pressures[:len(r_data.columns)-1], r_data['Temperature (K)'].values, r_data.iloc[:,1:].values)
# 	return r_from_PT_gas_func







# --------------------------------------------------------------------------------
# Create a 2D function to return P,T of a fluid at a given (rho,u) based on NIST data and using linear interpolation
def create_PT_from_ru_gas_func(fluid_props):
	P_data = pd.DataFrame(columns=['Internal Energy (l+v, kJ/kg)'])
	T_data = pd.DataFrame(columns=['Internal Energy (l+v, kJ/kg)'])
	densities = []		
	for dens in list(fluid_props.keys()):
		f = fluid_props[dens][['Internal Energy (l+v, kJ/kg)', 'Pressure (MPa)', 'Quality (l+v)']].rename(columns={'Pressure (MPa)': dens})		# Pull out Temperature, Enthalpy, and Phase data for each Pressure, and rename column to numercial pressure (in MPa)
		g = fluid_props[dens][['Internal Energy (l+v, kJ/kg)', 'Temperature (K)', 'Quality (l+v)']].rename(columns={'Temperature (K)': dens})	# Pull out Temperature, Enthalpy, and Phase data for each Pressure, and rename column to numercial pressure (in MPa)

		f = f.drop(f.index[f['Quality (l+v)'] != 'undefined']).drop('Quality (l+v)', axis=1) 													# Drop all rows that are not pertaining to vapor, then drop the 'Quality' column entirely so we're just left with Viscosity
		g = g.drop(g.index[g['Quality (l+v)'] != 'undefined']).drop('Quality (l+v)', axis=1) 													# Drop all rows that are not pertaining to vapor, then drop the 'Quality' column entirely so we're just left with Viscosity
		if not f.empty:
			P_data = pd.merge(P_data, f, how='outer', on='Internal Energy (l+v, kJ/kg)')														# Merge on 'Temperature (K)'. Non-overlapping values will have NaNs.
			T_data = pd.merge(T_data, g, how='outer', on='Internal Energy (l+v, kJ/kg)')														# Merge on 'Temperature (K)'. Non-overlapping values will have NaNs.

			P_data = P_data.sort_values('Internal Energy (l+v, kJ/kg)').reset_index(drop=True)													# Sort by volume and reset index
			T_data = T_data.sort_values('Internal Energy (l+v, kJ/kg)').reset_index(drop=True)													# Sort by volume and reset index

			P_data = P_data.set_index('Internal Energy (l+v, kJ/kg)').interpolate(method='index', limit_area='inside').reset_index()			# Fill in any missing values that fall between tabulated values
			T_data = T_data.set_index('Internal Energy (l+v, kJ/kg)').interpolate(method='index', limit_area='inside').reset_index()

			densities.append(fluid_props[dens]['Density (l, kg/m3)'][f.index[0]])												# Interpolate across the non-overlapping values to fill in any NaNs that lie BETWEEN data points.
		else:
			pass

	start = 350
	stop = P_data['Internal Energy (l+v, kJ/kg)'].iloc[-1]
	step = P_data['Internal Energy (l+v, kJ/kg)'].diff().median()
	extrap_range = np.arange(start, stop, step)

	P_data_new = pd.DataFrame(columns=['Internal Energy (l+v, kJ/kg)'])
	T_data_new = pd.DataFrame(columns=['Internal Energy (l+v, kJ/kg)'])

	for dens in P_data.columns[1:]:
		P_data_slice = P_data[['Internal Energy (l+v, kJ/kg)', dens]].dropna()													# Get the x, y values over the range you wish to use as the basis for your interpolation (and drop NaNs)
		T_data_slice = T_data[['Internal Energy (l+v, kJ/kg)', dens]].dropna()													# Get the x, y values over the range you wish to use as the basis for your interpolation (and drop NaNs)
		
		f = interp1d(P_data_slice['Internal Energy (l+v, kJ/kg)'], P_data_slice[dens].apply(np.exp), fill_value='extrapolate')	# Make your interpolation function, first applying np.exp() to all the values to make the extrapolation function from
		g = interp1d(T_data_slice['Internal Energy (l+v, kJ/kg)'], T_data_slice[dens], fill_value='extrapolate')	# Make your interpolation function, first applying np.exp() to all the values to make the extrapolation function from

		P_data_fill = pd.DataFrame({'Internal Energy (l+v, kJ/kg)': extrap_range, dens:np.log(f(extrap_range))})
		T_data_fill = pd.DataFrame({'Internal Energy (l+v, kJ/kg)': extrap_range, dens: g(extrap_range)})
		
		P_data_new = P_data_new.merge(P_data_fill, how='outer', on='Internal Energy (l+v, kJ/kg)')
		T_data_new = T_data_new.merge(T_data_fill, how='outer', on='Internal Energy (l+v, kJ/kg)')

	P_from_ru_func = interp2d(densities[:len(P_data_new.columns)-1], P_data_new['Internal Energy (l+v, kJ/kg)'].values, P_data_new.iloc[:,1:].values)
	T_from_ru_func = interp2d(densities[:len(T_data_new.columns)-1], T_data_new['Internal Energy (l+v, kJ/kg)'].values, T_data_new.iloc[:,1:].values)

	return P_from_ru_func, T_from_ru_func



# --------------------------------------------------------------------------------
# Create a 2D function to return P,T of a fluid at a given (rho,h) based on NIST data and using linear interpolation
def create_PT_from_rh_gas_func(fluid_props):
	P_data = pd.DataFrame(columns=['Enthalpy (l+v, kJ/kg)'])
	T_data = pd.DataFrame(columns=['Enthalpy (l+v, kJ/kg)'])
	densities = []		
	for dens in list(fluid_props.keys()):
		f = fluid_props[dens][['Enthalpy (l+v, kJ/kg)', 'Pressure (MPa)', 'Quality (l+v)']].rename(columns={'Pressure (MPa)': dens})		# Pull out Temperature, Enthalpy, and Phase data for each Pressure, and rename column to numercial pressure (in MPa)
		g = fluid_props[dens][['Enthalpy (l+v, kJ/kg)', 'Temperature (K)', 'Quality (l+v)']].rename(columns={'Temperature (K)': dens})	# Pull out Temperature, Enthalpy, and Phase data for each Pressure, and rename column to numercial pressure (in MPa)

		f = f.drop(f.index[f['Quality (l+v)'] != 'undefined']).drop('Quality (l+v)', axis=1) 													# Drop all rows that are not pertaining to vapor, then drop the 'Quality' column entirely so we're just left with Viscosity
		g = g.drop(g.index[g['Quality (l+v)'] != 'undefined']).drop('Quality (l+v)', axis=1) 													# Drop all rows that are not pertaining to vapor, then drop the 'Quality' column entirely so we're just left with Viscosity
		if not f.empty:
			P_data = pd.merge(P_data, f, how='outer', on='Enthalpy (l+v, kJ/kg)')														# Merge on 'Temperature (K)'. Non-overlapping values will have NaNs.
			T_data = pd.merge(T_data, g, how='outer', on='Enthalpy (l+v, kJ/kg)')														# Merge on 'Temperature (K)'. Non-overlapping values will have NaNs.

			P_data = P_data.sort_values('Enthalpy (l+v, kJ/kg)').reset_index(drop=True)													# Sort by volume and reset index
			T_data = T_data.sort_values('Enthalpy (l+v, kJ/kg)').reset_index(drop=True)													# Sort by volume and reset index

			P_data = P_data.set_index('Enthalpy (l+v, kJ/kg)').interpolate(method='index', limit_area='inside').reset_index()			# Fill in any missing values that fall between tabulated values
			T_data = T_data.set_index('Enthalpy (l+v, kJ/kg)').interpolate(method='index', limit_area='inside').reset_index()

			densities.append(fluid_props[dens]['Density (l, kg/m3)'][f.index[0]])												# Interpolate across the non-overlapping values to fill in any NaNs that lie BETWEEN data points.
		else:
			pass

	start = 350
	stop = P_data['Enthalpy (l+v, kJ/kg)'].iloc[-1]
	step = P_data['Enthalpy (l+v, kJ/kg)'].diff().median()
	extrap_range = np.arange(start, stop, step)

	P_data_new = pd.DataFrame(columns=['Enthalpy (l+v, kJ/kg)'])
	T_data_new = pd.DataFrame(columns=['Enthalpy (l+v, kJ/kg)'])

	for dens in P_data.columns[1:]:
		P_data_slice = P_data[['Enthalpy (l+v, kJ/kg)', dens]].dropna()													# Get the x, y values over the range you wish to use as the basis for your interpolation (and drop NaNs)
		T_data_slice = T_data[['Enthalpy (l+v, kJ/kg)', dens]].dropna()													# Get the x, y values over the range you wish to use as the basis for your interpolation (and drop NaNs)
		
		f = interp1d(P_data_slice['Enthalpy (l+v, kJ/kg)'], P_data_slice[dens].apply(np.exp), fill_value='extrapolate')	# Make your interpolation function, first applying np.exp() to all the values to make the extrapolation function from
		g = interp1d(T_data_slice['Enthalpy (l+v, kJ/kg)'], T_data_slice[dens], fill_value='extrapolate')	# Make your interpolation function, first applying np.exp() to all the values to make the extrapolation function from

		P_data_fill = pd.DataFrame({'Enthalpy (l+v, kJ/kg)': extrap_range, dens:np.log(f(extrap_range))})
		T_data_fill = pd.DataFrame({'Enthalpy (l+v, kJ/kg)': extrap_range, dens: g(extrap_range)})
		
		P_data_new = P_data_new.merge(P_data_fill, how='outer', on='Enthalpy (l+v, kJ/kg)')
		T_data_new = T_data_new.merge(T_data_fill, how='outer', on='Enthalpy (l+v, kJ/kg)')

	P_from_rh_func = interp2d(densities[:len(P_data_new.columns)-1], P_data_new['Enthalpy (l+v, kJ/kg)'].values, P_data_new.iloc[:,1:].values)
	T_from_rh_func = interp2d(densities[:len(T_data_new.columns)-1], T_data_new['Enthalpy (l+v, kJ/kg)'].values, T_data_new.iloc[:,1:].values)

	return P_from_rh_func, T_from_rh_func




# --------------------------------------------------------------------------------
# Create a 2D function to return P of a fluid at a given (rho,T) based on NIST data and using linear interpolation
def create_P_from_rT_gas_func(fluid_props):
	P_data = pd.DataFrame(columns=['Temperature (K)'])
	densities = []		
	for dens in list(fluid_props.keys()):
		f = fluid_props[dens][['Temperature (K)', 'Pressure (MPa)', 'Quality (l+v)']].rename(columns={'Pressure (MPa)': dens})		# Pull out Temperature, Enthalpy, and Phase data for each Pressure, and rename column to numercial pressure (in MPa)
		f = f.drop(f.index[f['Quality (l+v)'] != 'undefined']).drop('Quality (l+v)', axis=1) 										# Drop all rows that are not pertaining to vapor, then drop the 'Quality' column entirely so we're just left with Viscosity
		if not f.empty:
			P_data = pd.merge(P_data, f, how='outer', on='Temperature (K)')												# Merge on 'Temperature (K)'. Non-overlapping values will have NaNs.
			P_data = P_data.sort_values('Temperature (K)').reset_index(drop=True)											# Sort by volume and reset index
			P_data.interpolate(inplace=True)																			# Fill in any missing values that fall between tabulated values

			densities.append(fluid_props[dens]['Density (l, kg/m3)'][f.index[0]])												# Interpolate across the non-overlapping values to fill in any NaNs that lie BETWEEN data points.
		else:
			pass
	P_data = P_data.fillna(method='bfill', axis=1)																		# We can't interpolate over NaN's, and it throws a NaN for all the densities where there is no viscosity data (or when a phase change occurs), so...just backfill them for the time being
	P_data = P_data.dropna(axis=0, how='any', thresh=2)																	# Drop any rows that are completely NaN'ed out
	P_data = P_data.dropna(axis=1, how='all')																			# Drop the remaining all-Nan columns
	P_from_rT_func = interp2d(densities[:len(P_data.columns)-1], P_data['Temperature (K)'].values, P_data.iloc[:,1:].values)
	return P_from_rT_func




# --------------------------------------------------------------------------------
# Create a 2D function to return P,T of a fluid at a given (rho,h) based on NIST data and using linear interpolation
# def create_PT_from_rh_gas_func(fluid_props):
# 	P_data = pd.DataFrame(columns=['Enthalpy (l+v, kJ/kg)'])
# 	T_data = pd.DataFrame(columns=['Enthalpy (l+v, kJ/kg)'])
# 	densities = []		
# 	for dens in list(fluid_props.keys()):
# 		f = fluid_props[dens][['Enthalpy (l+v, kJ/kg)', 'Pressure (MPa)', 'Quality (l+v)']].rename(columns={'Pressure (MPa)': dens})		# Pull out Temperature, Enthalpy, and Phase data for each Pressure, and rename column to numercial pressure (in MPa)
# 		g = fluid_props[dens][['Enthalpy (l+v, kJ/kg)', 'Temperature (K)', 'Quality (l+v)']].rename(columns={'Temperature (K)': dens})		# Pull out Temperature, Enthalpy, and Phase data for each Pressure, and rename column to numercial pressure (in MPa)

# 		f = f.drop(f.index[f['Quality (l+v)'] != 'undefined']).drop('Quality (l+v)', axis=1) 										# Drop all rows that are not pertaining to vapor, then drop the 'Quality' column entirely so we're just left with Viscosity
# 		g = g.drop(g.index[g['Quality (l+v)'] != 'undefined']).drop('Quality (l+v)', axis=1) 										# Drop all rows that are not pertaining to vapor, then drop the 'Quality' column entirely so we're just left with Viscosity
# 		if not f.empty:
# 			P_data = pd.merge(P_data, f, how='outer', on='Enthalpy (l+v, kJ/kg)')												# Merge on 'Temperature (K)'. Non-overlapping values will have NaNs.
# 			T_data = pd.merge(T_data, g, how='outer', on='Enthalpy (l+v, kJ/kg)')												# Merge on 'Temperature (K)'. Non-overlapping values will have NaNs.

# 			P_data = P_data.sort_values('Enthalpy (l+v, kJ/kg)').reset_index(drop=True)											# Sort by volume and reset index
# 			T_data = T_data.sort_values('Enthalpy (l+v, kJ/kg)').reset_index(drop=True)											# Sort by volume and reset index

# 			P_data.interpolate(inplace=True)																			# Fill in any missing values that fall between tabulated values
# 			T_data.interpolate(inplace=True)

# 			densities.append(fluid_props[dens]['Density (l, kg/m3)'][f.index[0]])												# Interpolate across the non-overlapping values to fill in any NaNs that lie BETWEEN data points.
# 		else:
# 			pass
# 	P_data = P_data.fillna(method='bfill', axis=1)																		# We can't interpolate over NaN's, and it throws a NaN for all the densities where there is no viscosity data (or when a phase change occurs), so...just backfill them for the time being
# 	T_data = T_data.fillna(method='bfill', axis=1)																		# We can't interpolate over NaN's, and it throws a NaN for all the densities where there is no viscosity data (or when a phase change occurs), so...just backfill them for the time being

# 	P_data = P_data.dropna(axis=0, how='any', thresh=2)																	# Drop any rows that are completely NaN'ed out
# 	T_data = T_data.dropna(axis=0, how='any', thresh=2)																	# Drop any rows that are completely NaN'ed out

# 	P_data = P_data.dropna(axis=1, how='all')																			# Drop the remaining all-Nan columns
# 	T_data = T_data.dropna(axis=1, how='all')																			# Drop the remaining all-Nan columns

# 	P_from_rh_func = interp2d(densities[:len(P_data.columns)-1], P_data['Enthalpy (l+v, kJ/kg)'].values, P_data.iloc[:,1:].values)
# 	T_from_rh_func = interp2d(densities[:len(T_data.columns)-1], T_data['Enthalpy (l+v, kJ/kg)'].values, T_data.iloc[:,1:].values)
# 	return P_from_rh_func, T_from_rh_func



def check_ru_phase(rho, u, phase_data, P_from_ru_func, T_from_ru_func, P_trip, T_trip):
	# Given a rho, what is the corresponding u that would constitute a phase change?
	idx_of_closest_under = abs(phase_data[phase_data['Density, v (kg/m^3)'] < rho]['Density, v (kg/m^3)'] - rho).idxmin()
	idx_of_closest_over = abs(phase_data[phase_data['Density, v (kg/m^3)'] > rho]['Density, v (kg/m^3)'] - rho).idxmin()

	u_x1 = phase_data['Internal Energy, v (kJ/kg)'].iloc[idx_of_closest_under]
	u_x2 = phase_data['Internal Energy, v (kJ/kg)'].iloc[idx_of_closest_over]
	u_y1 = phase_data['Density, v (kg/m^3)'].iloc[idx_of_closest_under]
	u_y2 = phase_data['Density, v (kg/m^3)'].iloc[idx_of_closest_over]

	# Linearly interpolate
	m = (u_x2-u_x1)/(u_y2-u_y1)
	u_crrspnd_to_rho = m*(rho - u_y1) + u_x1

	# And given a u, what is the corresponding rho that would constitute a phase change?
	idx_of_closest_under = abs(phase_data[phase_data['Internal Energy, v (kJ/kg)'] < u]['Internal Energy, v (kJ/kg)'] - u).idxmin()
	idx_of_closest_over = abs(phase_data[phase_data['Internal Energy, v (kJ/kg)'] > u]['Internal Energy, v (kJ/kg)'] - u).idxmin()

	r_x1 = phase_data['Internal Energy, v (kJ/kg)'].iloc[idx_of_closest_under]
	r_x2 = phase_data['Internal Energy, v (kJ/kg)'].iloc[idx_of_closest_over]
	r_y1 = phase_data['Density, v (kg/m^3)'].iloc[idx_of_closest_under]
	r_y2 = phase_data['Density, v (kg/m^3)'].iloc[idx_of_closest_over]

	m = (r_y2-r_y1)/(r_x2-r_x1)
	rho_crrspnd_to_u = m*(u - r_x1) + r_y1

	if rho < rho_crrspnd_to_u and u > u_crrspnd_to_rho:
		phase = 'vapor'
	elif rho > rho_crrspnd_to_u or u < u_crrspnd_to_rho:
		# Is it above or below the triple point?
		# Need a PT from ru function. Don't I already have one?
		# YES I DO!
		P_from_ru = P_from_ru_func(rho, u)[0]*1000000
		T_from_ru = T_from_ru_func(rho, u)[0]
		if P_from_ru < P_trip and T_from_ru < T_trip:
			phase = 'solid-vapor two-phase'
		elif P_from_ru > P_trip and T_from_ru > T_trip:
			phase = 'liquid-vapor two-phase'
		else:
			phase = 'what!?'

		# Temperature (used to estimate Enthalpy of Sublimation)
		# T_y1 = phase_data['Temperature (K)'].iloc[idx_of_closest_under]
		# T_y2 = phase_data['Temperature (K)'].iloc[idx_of_closest_over]
		
		# m = (T_y2-T_y1)/(r_x2-r_x1)
		# T_crrspnd_to_u = m*(u - r_x1) + T_y1
		
		# enth_of_sub = y = -0.2985*T_crrspnd_to_u + 644.19
		
		# Enthalpy
		# h_y1 = phase_data['Enthalpy (kJ/kg)'].iloc[idx_of_closest_under]
		# h_y2 = phase_data['Enthalpy (kJ/kg)'].iloc[idx_of_closest_over]

		# m = (h_y2-h_y1)/(r_x2-r_x1)
		# h_crrspnd_to_u = m*(u - r_x1) + h_y1

		# quality = h_crrspnd_to_u

	return phase