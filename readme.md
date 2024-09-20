The official repository for the ATOM (Arc Trajectory Optimization Method) algorithm. Given a Proton Arc Therapy plan for cancer treatment, ATOM finds the fastest delivery of the radiotherapy plan to the patient, while respecting the mechanical machine constraints.

The radiotherapy plan consists of a sequence of so called energy layers. Each such layer has a minimum delivery duration (irradiation time) and an angular range in which it has to be delivered (gantry window). The machine will have to change energy between the layers, and this change is usually not instantanious.

This means that the inputs to the problem are:
- A list of irradiation times [s].
- The size of the gantry windows [deg].
- A list of distances between the middle angles of the gantry windows [deg].
- A list of energy change times [s].
Note that the last two lists will be one element shorter than the first, since the energy change occurs between the energy layers.

Furthermore, the movement constraints of the machine needs to be given:
- Maximum velocity. It is assumed that the minimum velocity is 0.
- Maximum acceleration.
- Minimum acceleration.
- Maximum jerk. The minimum is assumed to be the maximum but with negative sign.

```
from ATOM import atom

irr_times = [0.1, 0.2, 0.15]
elsts = [0.9, 1.5]
angle_distances = [2.0, 2.0]  # The angles are [0, 2, 4]
maximum_window_size = 1.0  # Which means that the ranges are (-0.5, 0.5), (1.5, 2.5), (3.5, 4.5). However, the gantry will never pass outside [0, 4]

parameters = {"v_max": 5.0, "a_max": 0.5, "a_min": -0.5, "j_max": 0.5}
delivery_time, discrete_vels = atom(irr_times, elsts, angle_distances, maximum_window_size, parameters)
vels = [discrete_vels[i] * parameters["v_max"] / (256-1) for i in range(len(discrete_vels))]

print("Delivery time [s]:", delivery_time )
print("Velocities [angle / s]:", vels)
```

NOTE: It requires the pip library ruckig, and has been tested using version 0.9.2.
