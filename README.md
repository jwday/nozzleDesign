# nozzleDesign
Simulates 1D gas flow through a converging-diverging (C-D) nozzle for CO2 and R134a, using a fixed-volume plenum chamber as the propellant source, with starting pressures at 100 and 87 psig, respectively. Simulation can be performed either using isentropic relations or using interpolated/extrapolated NIST data. Simulation outputs .csv file with simulation results. Plots are generated in separate files.

Many output options are available based on what properties the user wishes to view. Selection is done in-code by simply commenting out the undesired properties. The user may select from the following list:

- Plenum Total Pressure, Total Temperature, Total Density, Total Specific Volume, Specific Enthalpy, Specific Internal Energy, Vapor Quality
- Reynolds No., Nusselt No., Prandtl No., Viscosity Upstream of Valve
- Valve Wall Temperature
- Nozzle Inlet Total Pressure, Total Temperature, Static Temperature, Mach No.
- Nozzle Throat Flow Pressure, Temperature, Density, Reynold's No., Mach No., Velocity
- Mass Flow Rate
- Thrust contributions from mass flow rate (momentum) and pressure differential
- Instantaneous Thrust, Ideal Thrust Coefficient, Viscous Loss Percentage, Effective Thrust Coefficient, Effective Thrust, Time Average Thrust
- Net Impulse, Net Effective Impulse, Specific Impulse
- Area Ratio at Shock Location
- Saturated Vapor Pressure at Current Temperature
- Saturated Vapor Temperature at Current Pressure
- Saturated Vapor Pressure at Exit
- Saturated Vapor Temperature at Exit


## Running this code
Ideally, a user should be able to run nozzle_master.py directly from your favorite editor with Python extension (such as VSCode). Currently support has not been added for command line execution. Currently only CO2 and R134a are supported. The following libraries are required:

- matplotlib
- scipy
- pandas
- numpy
- seaborn

This was developed using Anaconda 3 (Python 3.8.5). If you are using Windows Subsystem for Linux, you must have a working X server in order to view the plots.

## Operation
This code integrates over time to estimate the various gas properties (pressure, temperature, thrust, etc.), then the script will generate a number of plots:

1. Time-rate of change of Pressure (P) vs. temperature (T), and density (rho) vs. internal energy (u) plotted against real gas data to identify the point where phase change (condensation) begins within the plenum due to the rapidly decreasing pressure.
2. Throat Reynold's Number (based on pipe flow), used to identify the point of transition from turbulent to laminar flow (~< 4E3)
3. Upstream and Downstream flow temperatures as it enters and exits the valve. These results indicate that there is very little heat transfered from the valve to the flow.
4. Plenum pressure, thrust (both instantaneous and time-averaged), and net impulse vs. time for a given gas
5. Time-logarithmic comparison of in-space vs. on-ground properties as specified by whichever properties are uncommented in the 'data' dictionary object (line 1009).
6. Flow properties over the length of the specified nozzle geometry taken at three points in time

### Ex. Performance vs. Time
Nozzle performance is illustrated by showing the time-rate of change of plenum pressure, corresponding thrust (both instantaneous and time-average), and the net impulse generated over the duration of the discharge.

<img src="https://github.com/jwday/nozzleDesign/blob/master/Sim_Thrust_and_Impulse_CO2.png" alt="Sim CO2 P and Thr" height="400">

### Ex. Thrust Coefficient vs. Time
One of the many parameters which can be viewed is the Thrust Coefficient, including viscous losses as modeled based on the work by [Spisz et. al from NASA Lewis (Glenn) Research Center](https://ntrs.nasa.gov/api/citations/19650027295/downloads/19650027295.pdf).

<img src="https://github.com/jwday/nozzleDesign/blob/master/Thrust_Coeff_and_Visc_Losses_CO2.png" alt="Sim CO2 Cf vs Time" height="400">

#### Ex. Density vs. Internal Energy (rho-u) Plot
The change in density and internal energy was used to determine the phase and estimate the plenum gas quality (X) as it transitioned into the two-phase region. Of course it was later found that despite a phase change potentially taking place, the quality remained very high (>99%) for the duration of the discharge, indicating that the isentropic assumptions which were originally used to model the gas would still remain valid and continue to provide a good estimate.

<img src="https://github.com/jwday/nozzleDesign/blob/master/rho-u_CO2_Phase-Region-Plot.png" alt="Sim CO2 Phase Change" height="400">
