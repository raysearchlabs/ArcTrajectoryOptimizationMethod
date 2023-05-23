from ATOM import atom
from random import random, seed

seed(0)
n = 10
vel_res = 256
irr_times = [random() for _ in range(n)]
elsts = [random() for _ in range(n-1)]
angle_distances = [1.0] * (n - 1)
maximum_window_size = 0.5

parameters = {"v_max": 5.0, "a_max": 0.5, "a_min": -0.5, "j_max": 0.5}
delivery_time, vels = atom(irr_times, elsts, angle_distances, maximum_window_size, parameters, vel_res)

print("Delivery time [s]:", delivery_time )
print("Velocities [angle / s]:", vels)
