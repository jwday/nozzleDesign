# Thrust vs. Burn Time for 0.25 N-s impulse requirement 
import numpy as np
import matplotlib.pyplot as plt

t = np.linspace(0, 60, 6001)

req_impulse = 0.25  # N-s
max_time = 56  # Seconds
min_const_thrust = req_impulse/max_time  # N

num_thrust_pulses = 4
thrust_pulse_impulse = req_impulse/num_thrust_pulses  # This is equivalent to the area under the curve
pulse_startstop_times = np.linspace(0, max_time, 2*num_thrust_pulses)
pulse_startstop_times = np.append(pulse_startstop_times, [60])

# Assume pulses decrease linearly in thrust over time (triangular pulse)
thrust_pulse_duration = max_time/(2*num_thrust_pulses-1)  # Assuming plenum refil is equal in time to pulse duration
max_pulse_thrust = 2*thrust_pulse_impulse/thrust_pulse_duration  # The "height" of the pulse triangle

pulse_thrust = []
n = 0
for i in t:
    if n < num_thrust_pulses:
        if i >= pulse_startstop_times[2*n] and i < pulse_startstop_times[2*n+1]:
            pulse_thrust.append( i*(-max_pulse_thrust/thrust_pulse_duration) + max_pulse_thrust*(2*n+1))
        elif i >= pulse_startstop_times[2*n+1] and i < pulse_startstop_times[2*n+2]:
            pulse_thrust.append(0)
        else:
            n += 1
            pulse_thrust.append( i*(-max_pulse_thrust/thrust_pulse_duration) + max_pulse_thrust*(2*n+1))  # To catch the top of the new thrust curve at the beginning of each pulse


plt.plot(t, pulse_thrust)
plt.show()