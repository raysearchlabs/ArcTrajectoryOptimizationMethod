"""
Script to calculate the delivery time of a dynamic proton arc plan from RayStation.

Note that these plans are currently (June 2023) not clinically supported, and a special research build is required.

This script might take 1-5 min, depending on the plan. Patience is a virtue. Or so I've heard.
"""


#########################################################################
# PARAMETERS!
# Change them to accurately reflect the machine of interest.
#########################################################################

# NOTE: Make sure that the maximum window size is smaller than the angular spacing between 
MAX_WINDOW_SIZE = 1.99  # TODO: In the future this value will be read automatically from RayStation. 

# Maximum velocity, acceleration and jerk of the gantry.
VEL_MAX = 5.0
ACC_MAX = 0.5
JERK_MAX = 0.25

# Energy switching times
DOWN_SWITCH_TIME = 0.5
UP_SWITCH_TIME = 5.0

# Irradiation times
TIME_PER_SPOT_SWITCH = 0.008
SPOT_DELIVERY_S_PER_MU = 0.005

# This is the velocity resolution used in ATOM.
#Higher value gives more accurate result, but increases calculation time.
VEL_RES = 256

#########################################################################
# END OF PARAMETERS!
# Thanks for your attention (:
#########################################################################

from ATOM import atom

from connect import *


def calculate_spot_scanning_time(beamMU, segment):
	# Constant delivery rate
	totalTime = 0
	for w in segment.Spots.Weights:
		spotMu = beamMU * w
		spot_time = SPOT_DELIVERY_S_PER_MU * spotMu
		totalTime += spot_time

	return totalTime


def calculat_spot_switching_time(segment):
	# Constant spot switching time
	numberOfSpots = len(segment.Spots.Weights)
	numberOfSpotSwitches = numberOfSpots - 1 if numberOfSpots > 0 else 0
	return TIME_PER_SPOT_SWITCH * numberOfSpotSwitches


def estimate_switch_time(fromNrj, toNrj):
	# Constant up-switching time & constant down-switching
	if abs(fromNrj - toNrj) < 1.0e-8:
		return 0.0
	elif fromNrj < toNrj:
		return UP_SWITCH_TIME
	else:
		return DOWN_SWITCH_TIME


def calc_is_clockwise(fromAngle, toAngle):
	assert fromAngle >= 0
	assert toAngle >= 0
	assert fromAngle < 360
	assert toAngle < 360
	assert fromAngle != toAngle

	shiftedToAngle = toAngle - fromAngle

	while shiftedToAngle < 0:
		shiftedToAngle += 360

	return shiftedToAngle < 180


plan = get_current("Plan")

beams = list(get_current("BeamSet").Beams)[:]
parameters = {"v_max": VEL_MAX, "a_max": ACC_MAX, "a_min": -ACC_MAX, "j_max": JERK_MAX}

for b in range(len(beams)):
	total_time = 0.0
	n = len(beams[b].Segments)

	nrjs = [beams[b].Segments[i].NominalEnergy for i in range(n)]

	elsts = [estimate_switch_time(nrjs[i], nrjs[i+1]) for i in range(len(nrjs)-1)]
	irr_times = [calculate_spot_scanning_time(beams[b].BeamMU, beams[b].Segments[i]) + calculat_spot_switching_time(beams[b].Segments[i]) for i in range(n)]

	angles = [beams[b].Segments[i].IonArcSegmentProperties.DeltaGantryAngle for i in range(n)]

	angle_distances = [min(abs(angles[i] - angles[i+1]), abs(360-abs(angles[i] - angles[i+1]))) for i in range(len(angles) - 1)]
	assert min(angle_distances) > 0  # TODO: Handle case where an angle has many energy layers.
	assert max(angle_distances) < 180  # We don't support this case.
	assert len(angle_distances) > 2
	is_clockwise = [calc_is_clockwise(angles[i], angles[i+1]) for i in range(n-1)]
	direction_change = [is_clockwise[i] and not is_clockwise[i+1] for i in range(n-2)] + [True]

	current_idx = 0
	while any(direction_change):
		next_dir_change = direction_change.index(True)
		direction_change[next_dir_change] = False

		irr_times_local = irr_times[current_idx:next_dir_change]
		elsts_local = elsts[current_idx:next_dir_change - 1]
		angle_distances_local = angle_distances[current_idx:next_dir_change - 1]

		delivery_time, discrete_vels = atom(irr_times_local, elsts_local, angle_distances_local, MAX_WINDOW_SIZE, parameters, VEL_RES)
		vels = [discrete_vels[i] * parameters["v_max"] / (VEL_RES - 1) for i in range(len(discrete_vels))]

		current_idx = next_dir_change

		total_time += delivery_time
	print("TOTAL TIME: ", round(total_time, 3), "[s] for beam", b+1)
