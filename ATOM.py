from math import sqrt
from ruckig import InputParameter, Result, Ruckig, Trajectory  # pip install ruckig


class State():
	def __init__(self, v_idx, ang_idx, disc_vels, irr_times, is_final=False):
		self.v_idx = v_idx
		self.ang_idx = ang_idx
		self.v = disc_vels[v_idx]
		self.vel_square = self.v * self.v
		if is_final:
			assert self.v == 0
			self.half_window = 0
		else:
			self.half_window = self.v * irr_times[self.ang_idx] * 0.5

		self.key = str(self.v_idx) + "_" + str(ang_idx)

	def __repr__(self):
		return self.key

	def __str__(self):
		return self.key

	def __hash__(self):
		return hash(self.key)


def linspace(a, b, n):
	assert a < b
	assert n > 1
	delta = (b - a) / (n-1)

	return [a + i*delta for i in range(n)]


def calc_time_between_segments(irr_times, switch_times, angle_distances, max_window, idx, v0, v1, parameters):
	assert idx > 0
	idx0 = idx - 1
	idx1 = idx

	window_size_0 = v0 * irr_times[idx0]
	window_size_1 = v1 * irr_times[idx1]

	if window_size_0 > max_window or window_size_1 > max_window:
		return float('inf')
	remaining_angle = angle_distances[idx0] - (window_size_0+window_size_1)/2
	assert remaining_angle > 0

	elst = switch_times[idx0]

	otg = Ruckig(1)
	inp = InputParameter(1)
	trajectory = Trajectory(1)

	inp.current_position = [0.0]
	inp.current_velocity = [v0]
	inp.current_acceleration = [0.0]
	
	inp.target_velocity = [v1]
	inp.target_position = [remaining_angle]
	inp.target_acceleration = [0.0]
	
	inp.max_velocity = [parameters["v_max"]]
	inp.max_acceleration = [parameters["a_max"]]
	inp.max_jerk = [parameters["j_max"]]

	inp.min_velocity = [-1.0e-16]  # This would be 0, but there are some rounding errors that cause issues, it seems.
	inp.min_acceleration = [parameters["a_min"]]

	inp.minimum_duration = elst
	if not otg.validate_input(inp, check_current_state_within_limits=True, check_target_state_within_limits=True):
		return float('inf')

	# Calculate the trajectory in an offline manner
	try:
		result = otg.calculate(inp, trajectory)
	except RuntimeError:
		return float('inf')
	if result == Result.ErrorInvalidInput:
		return float('inf')

	return trajectory.duration + irr_times[idx0]


def get_neighs(all_states, current, has_been_visited):
	# TODO: This is a disgrace. Make it faster.
	return list(filter(lambda el: not has_been_visited[current.ang_idx+1][all_states[current.ang_idx+1].index(el)], all_states[current.ang_idx+1]))


def find_current(open_set, f_score):
	# TODO: This might be too slow. Use a more clever data structure. A fibonacci heap, probably. 
	best_state = None
	best_f_score = float('inf')
	for state in open_set:
		if f_score[state.ang_idx][state.v_idx] < best_f_score:
			best_f_score = f_score[state.ang_idx][state.v_idx]
			best_state = state
		elif f_score[state.ang_idx][state.v_idx] == best_f_score:
			if state.v > best_state.v:
				best_f_score = f_score[state.ang_idx][state.v_idx]
				best_state = state

	assert best_f_score < float('inf')
	return best_state


def atom(irr_times, elsts, angle_distances, maximum_window_size, parameters, vel_res=256):
	n = len(irr_times)
	assert n > 1
	assert n - 1 == len(elsts)
	assert n - 1 == len(angle_distances)
	assert maximum_window_size < min(angle_distances)
	assert min(irr_times) >= 0
	assert min(elsts) >= 0
	assert vel_res > 1

	disc_vels = linspace(0, parameters["v_max"], vel_res)
	max_v = parameters["v_max"]
	max_acc = parameters["a_max"]
	min_acc = parameters["a_min"]

	assert max_v > 0
	assert min_acc < 0
	assert max_acc > 0
	assert parameters["j_max"] > 0

	heuristic_from_ang_idx = [sum(irr_times[i:]) - irr_times[-1] + sum(elsts[i:]) for i in range(n)]

	came_from = {}

	g_score = []
	f_score = []

	all_states = []
	initial_state = State(0, 0, disc_vels, irr_times)
	final_state = State(0, n - 1, disc_vels, irr_times, is_final=True)
	all_states.append([initial_state])
	has_been_visited = []
	g_score.append([0])
	f_score.append([heuristic_from_ang_idx[initial_state.ang_idx]])
	has_been_visited.append([False])
	for ang_idx in range(1, n - 1):
		states_for_angle = []
		has_been_visited_for_angle = []
		g_score_for_angle = []
		f_score_for_angle = []
		for v_idx in range(len(disc_vels)):
			state = State(v_idx, ang_idx, disc_vels, irr_times)
			states_for_angle.append(state)
			has_been_visited_for_angle.append(False)
			g_score_for_angle.append(float('inf'))
			f_score_for_angle.append(float('inf'))
		g_score.append(g_score_for_angle)
		f_score.append(f_score_for_angle)
		has_been_visited.append(has_been_visited_for_angle)
		all_states.append(states_for_angle)

	has_been_visited.append([False])
	all_states.append([final_state])
	g_score.append([float('inf')])
	f_score.append([float('inf')])

	open_set = set([initial_state])

	found_it = False
	while len(open_set) > 0:
		current = find_current(open_set, f_score)
		assert not has_been_visited[current.ang_idx][current.v_idx]
		has_been_visited[current.ang_idx][current.v_idx] = True
		if current == final_state:
			assert heuristic_from_ang_idx[current.ang_idx] == 0
			found_it = True
			break

		open_set.remove(current)
		current_g_score = g_score[current.ang_idx][current.v_idx]

		neighs = get_neighs(all_states, current, has_been_visited)
		for neigh in neighs:
			assert neigh.ang_idx == current.ang_idx + 1
		if len(neighs) > 0:
			angle_dist_between_center_points = angle_distances[current.ang_idx]

			local_max_vel = sqrt(2 * max_acc * (angle_dist_between_center_points - current.half_window) + current.vel_square)
			local_min_vel = 0 if 2 * min_acc * (angle_dist_between_center_points - current.half_window) + current.vel_square < 0 else sqrt(2 * min_acc * (angle_dist_between_center_points - current.half_window) + current.vel_square)

			g_scores_neighs = g_score[neighs[0].ang_idx]
			f_scores_neighs = f_score[neighs[0].ang_idx]
			elst = elsts[current.ang_idx]
			irr_time_from = irr_times[current.ang_idx]
			for neigh_state in neighs:
				if neigh_state.v > local_max_vel or neigh_state.v < local_min_vel:
					continue
				d = calc_time_between_segments(irr_times, elsts, angle_distances, maximum_window_size, neigh_state.ang_idx, current.v, neigh_state.v, parameters)
				if d != float('inf'):
					tentative_g_score = current_g_score + d
					if tentative_g_score <= g_scores_neighs[neigh_state.v_idx]:
						came_from[neigh_state] = current
						g_scores_neighs[neigh_state.v_idx] = tentative_g_score
						f_scores_neighs[neigh_state.v_idx] = tentative_g_score + heuristic_from_ang_idx[neigh_state.ang_idx]

						if has_been_visited[neigh_state.ang_idx][neigh_state.v_idx] == False and neigh_state not in open_set:
							open_set.add(neigh_state)

	assert found_it
	time_val = g_score[current.ang_idx][current.v_idx] + irr_times[-1]
	traj = [current.v_idx]
	while current in came_from:
		traj.append(max_v * came_from[current].v_idx / (vel_res - 1))
		current = came_from[current]
	traj = list(reversed(traj))
	return time_val, traj
