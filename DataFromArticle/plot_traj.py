from ruckig import InputParameter, Result, Ruckig, Trajectory  # pip install ruckig
import seaborn as sns

sns.set_style("darkgrid")

otg = Ruckig(1)
inp = InputParameter(1)
trajectory = Trajectory(1)

inp.current_position = [0.0]
inp.current_velocity = [0.7]
inp.current_acceleration = [0.0]
		
inp.target_velocity = [0.1]
inp.target_position = [1.0]
inp.target_acceleration = [0.0]
		
inp.max_velocity = [5.0]
inp.max_acceleration = [0.5]
inp.max_jerk = [0.5]

inp.min_velocity = [-1.0e-16]  # This would be 0, but there are some rounding errors that cause issues, it seems.
inp.min_acceleration = [-0.5]

#inp.minimum_duration = elst
if not otg.validate_input(inp, check_current_state_within_limits=True, check_target_state_within_limits=True):
	assert False

# Calculate the trajectory in an offline manner
try:
	result = otg.calculate(inp, trajectory)
except RuntimeError:
	assert False
if result == Result.ErrorInvalidInput:
	assert False

import numpy as np
import matplotlib.pyplot as plt
ts = np.linspace(0, trajectory.duration, 1000)
xs = []
accs = []
vs = []
for t in ts:
	new_position, new_velocity, new_acceleration = trajectory.at_time(t)
	xs.append(new_position)
	vs.append(new_velocity)
	accs.append(new_acceleration)

palette = sns.color_palette("rocket_r")

jerk = [(accs[i+1][0] - accs[i][0])/(ts[i+1] - ts[i]) for i in range(len(accs) - 1)]
jerk.append(jerk[-1])

plt.plot(ts, xs, label="Position", color=palette[0])
plt.plot(ts, vs, label="Velocity", color=palette[2], linestyle="--")
plt.plot(ts, accs, label="Acceleration", color=palette[3], linestyle=":")
plt.plot(ts, jerk, label="Jerk", color=palette[5], linestyle="-.")
#plt.plot(ts, xs, label="Position", color=palette[5])
#plt.plot(ts, vs, label="Velocity", color=palette[5], linestyle="--")
#plt.plot(ts, accs, label="Acceleration", color=palette[5], linestyle=":")
#plt.plot(ts, jerk, label="Jerk", color=palette[5], linestyle="-.")
plt.title("Example of Trajectory from Ruckig")
plt.xlabel("Time [s]")
plt.ylabel("Position, Velocity\n Acceleration & Jerk")
plt.legend()

plt.show()
