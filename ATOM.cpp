#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <ruckig/ruckig.hpp>

#include <chrono>
#include <cmath>
#include <utility>
#include <vector>
#include <cassert>
#include <set>
#include <limits>
#include <random>

using namespace std::chrono;
using namespace ruckig;

namespace {
#define VEL_DISC 256

class State {
  public:
	size_t vIdx;
	size_t angIdx;
	double v;
	double velSquare;

	size_t key;
	double halfWindow;

	State() :
		vIdx(0u),
		angIdx(0u),
		v(0.0),
		velSquare(0.0),
		key(0u),
		halfWindow(0.0)
	{};

	State(size_t vIdx_, size_t angIdx_, double v_, double irrTime) :
		vIdx(vIdx_),
		angIdx(angIdx_),
		v(v_),
		velSquare(v_ * v_),
		key(vIdx + angIdx * VEL_DISC),
		halfWindow(v_ * irrTime * 0.5 )
	{}

	bool operator==(const State& other) const {
		return key == other.key;
	}

	size_t operator() () const {
		return key;
	}

	bool operator< (const State& right) const {
		return key < right.key;
	}
};

double calcCost(
    const std::vector<double>& irrTimes,
    double elst,
    const std::vector<double>& angleDistances,
    double maximumWindowSize,
    size_t neighAngIdx,
    double currentV,
    double neighV,
    double vMax,
    double aMax,
    double jMax) {
	assert(neighAngIdx >= 1);
	const size_t currentAngleIdx = neighAngIdx - 1;

	const auto windowSize0 = currentV * irrTimes.at(currentAngleIdx);
	const auto windowSize1 = neighV * irrTimes.at(neighAngIdx);

	if (windowSize0 > maximumWindowSize || windowSize1 > maximumWindowSize) {
		return std::numeric_limits<double>::max();
	}

	const double remainingAngle = angleDistances[currentAngleIdx] -  0.5 * (windowSize0 + windowSize1);
	assert(remainingAngle > 0);

	InputParameter<1> input;
	input.current_position = { 0.0 };
	input.current_velocity = { currentV };
	input.current_acceleration = { 0.0 };

	input.target_position = { remainingAngle };
	input.target_velocity = { neighV };
	input.target_acceleration = { 0.0 };

	input.max_velocity = { vMax };
	input.max_acceleration = { aMax };
	input.max_jerk = { jMax };

	input.min_velocity = { -vMax};
	input.minimum_duration = elst;

	// We don't need to pass the control rate (cycle time) when using only offline features
	Ruckig<1> otg;
	Trajectory<1> trajectory;

	// Calculate the trajectory in an offline manner (outside of the control loop)
	const Result result = otg.calculate(input, trajectory);
	if (result != Result::Working ) {
		return std::numeric_limits<double>::max();
	}
	assert(result != Result::Finished);
	assert(trajectory.get_duration() >= elst - 0.00001);
	return trajectory.get_duration() + irrTimes.at(currentAngleIdx);
}

State findCurrent(std::unordered_map<size_t, State>& stateSet, std::set<std::pair<double, size_t>>& stateSetF) {
	size_t keyOfBest = stateSetF.begin()->second;
	const State bestStateCopy = stateSet.find(keyOfBest)->second;

	stateSet.erase(keyOfBest);
	stateSetF.erase(stateSetF.begin());
	return bestStateCopy;
}

std::vector<double> getDiscVels(double vMax) {
	std::vector<double> discVels(VEL_DISC);
	for ( size_t i = 0; i < VEL_DISC; ++i) {
		discVels.at(i) = i * vMax / (VEL_DISC - 1);
	}
	return discVels;
}
}

double atom(
    size_t n,
    std::vector<double>& irrTimes,
    std::vector<double>& elsts,
    std::vector<double>& angleDistances,
    double maximumWindowSize,
    double vMax,
    double aMax,
    double jMax) {
	const auto discVels = getDiscVels(vMax);

	double cumsumIrr = 0.0;
	double cumsumElst = 0.0;

	assert(!irrTimes.empty());
	assert(n > 0);
	const auto irrSum = std::accumulate( irrTimes.cbegin(), --irrTimes.cend(), 0.0 );  // NOTE: We're skipping the last element.
	const auto elstSum = std::accumulate(elsts.begin(), elsts.end(), 0.0);

	std::vector<double> heuristicFromAngIdx(n);
	for (size_t i = 0; i < n; ++i) {
		heuristicFromAngIdx.at(i) = elstSum + irrSum - cumsumIrr - cumsumElst;
		if (i + 1 != n) {
			cumsumElst += elsts.at(i);
			cumsumIrr += irrTimes.at(i);
		}
	}

	auto initialState = State(0, 0,  0.0, irrTimes.front());
	auto finalState = State(0, n - 1,  0.0, irrTimes.at(n - 1));

	std::vector<std::vector<double>> gScores = { {0.0} };
	std::vector<std::vector<double>> fScores = { {heuristicFromAngIdx.at(initialState.angIdx)} };
	std::vector<std::vector<State>> allStates = { {initialState} };
	std::vector<std::vector<bool>> hasBeenVisited = { {false} };

	for (size_t angIdx = 1; angIdx < n - 1; angIdx++) {
		std::vector<double> gScoreForAngle(VEL_DISC, std::numeric_limits<double>::max());
		std::vector<double> fScoreForAngle(VEL_DISC, std::numeric_limits<double>::max());
		std::vector<State> statesForAngle;
		statesForAngle.reserve(VEL_DISC);
		std::vector<bool> hasBeenVisitedForAngle(VEL_DISC, false);
		for (size_t vIdx = 0; vIdx < VEL_DISC; vIdx++) {
			auto state = State(vIdx, angIdx, discVels.at(vIdx), irrTimes.at(angIdx));
			statesForAngle.push_back(state);
		}
		gScores.push_back(gScoreForAngle);
		fScores.push_back(fScoreForAngle);
		hasBeenVisited.push_back(hasBeenVisitedForAngle);
		allStates.push_back(statesForAngle);
	}

	hasBeenVisited.push_back({ false });
	allStates.push_back({ finalState });

	gScores.push_back({ std::numeric_limits<double>::max() });
	fScores.push_back({ std::numeric_limits<double>::max() });

	std::unordered_map<size_t, State> stateSet = {{0, initialState}};
	std::set<std::pair<double, size_t>> stateSetF = {{ fScores.at(0).at(0), 0 }};

	std::unordered_map<size_t, State> bestPrev;
	std::set<std::pair<double, size_t>> openSetSorted = { {double{0}, size_t{0}} };

	bool foundIt = false;
	while ( !stateSet.empty()) {
		const auto currentState = findCurrent(stateSet, stateSetF);
		hasBeenVisited.at(currentState.angIdx).at(currentState.vIdx) = true;

		if (currentState.key == finalState.key) {
			foundIt = true;
			break;
		}
		assert(currentState.angIdx < n);
		assert(currentState.angIdx < gScores.size());
		assert(currentState.vIdx < gScores.at(currentState.angIdx).size());
		auto currentGScore = gScores.at(currentState.angIdx).at(currentState.vIdx);

		const auto angleDistBetweenCenterPoints = angleDistances.at(currentState.angIdx);
		const auto localMinVel =  -2.0 * aMax * (angleDistBetweenCenterPoints - currentState.halfWindow) + currentState.velSquare < 0.0 ?
		                          0.0  :
		                          sqrt( -2.0  * aMax * (angleDistBetweenCenterPoints - currentState.halfWindow) + currentState.velSquare);
		const auto localMaxVel = sqrt(2.0 * aMax * (angleDistBetweenCenterPoints - currentState.halfWindow) + currentState.velSquare);

		const auto elst = elsts.at(currentState.angIdx);
		const size_t a = currentState.angIdx + 1;
		auto& gScoresNeighs = gScores.at(a);
		auto& fScoresNeighs = fScores.at(a);

		for ( size_t i = 0; i < allStates.at(a).size(); ++i) {
			if ( hasBeenVisited.at(a).at(i)) {
				continue;
			}
			const auto& neighState = allStates.at(a).at(i);
			if (neighState.v < localMinVel) {
				continue;
			}
			if ( neighState.v > localMaxVel) {
				break;
			}

			const auto dist = calcCost(irrTimes,
			                           elst,
			                           angleDistances,
			                           maximumWindowSize,
			                           neighState.angIdx,
			                           currentState.v,
			                           neighState.v,
			                           vMax,
			                           aMax,
			                           jMax);

			if (dist < std::numeric_limits<double>::max()) {
				const auto tentativeGScore = currentGScore + dist;
				assert(neighState.vIdx < gScoresNeighs.size());

				if (tentativeGScore <= gScoresNeighs.at(neighState.vIdx) ) {
					assert(neighState.angIdx < heuristicFromAngIdx.size());

					const auto oldF = fScoresNeighs.at(neighState.vIdx);

					gScoresNeighs.at(neighState.vIdx) = tentativeGScore;
					fScoresNeighs.at(neighState.vIdx) = tentativeGScore + heuristicFromAngIdx.at(neighState.angIdx);

					assert(neighState.angIdx < hasBeenVisited.size());
					assert(neighState.vIdx < hasBeenVisited.at(neighState.angIdx).size());
					assert(!hasBeenVisited.at(neighState.angIdx).at(neighState.vIdx));
					if (stateSet.count(neighState.key) == 0) {
						stateSet.insert({ neighState.key, neighState });
					} else {
						stateSetF.erase({ oldF, neighState.key });
					}
					stateSetF.emplace(fScoresNeighs.at(neighState.vIdx), neighState.key);
					bestPrev.insert({ allStates.at(a).at(i).key, currentState });
				}
			}
		}
	}
	auto currentStateKey = finalState.key;
	while (bestPrev.count(currentStateKey) == 1) {
		const auto& currentState = bestPrev.at(currentStateKey);
		currentStateKey = currentState.key;
	}

	assert(foundIt);
	auto timeVal = gScores.at(finalState.angIdx).at(finalState.vIdx) + irrTimes.back();
	return timeVal;
}

int main() {
	srand(1);

	milliseconds ms = duration_cast<milliseconds>(
	                      system_clock::now().time_since_epoch());
	std::default_random_engine generator(0ULL);

	auto s = 0.0;
	for (size_t iter = 0; iter < 10; iter++) {
		size_t n_ = 360;
		auto ELST_down_ =  0.5;
		auto ELST_up_ =  10 * ELST_down_;
		std::vector<double> angleDistances_(n_ - 1, 2.0);
		auto maxWindow_ =  2.0;

		std::uniform_double_distribution<double> distributionElst(0.0, 1.0);
		std::uniform_double_distribution<double> distributionIrrs(0.0, 1.26);

		std::vector<double> irrTimes_(n_);
		for (size_t i = 0; i < n_; i++) {
			irrTimes_.at(i) = distributionIrrs(generator);
		}

		std::vector<double> switchTimes_(n_ - 1);
		for (size_t i = 0; i < n_ - 1; i++) {
			if (distributionElst(generator) <  0.1 ) {
				switchTimes_.at(i) = ELST_up_;
			} else {
				switchTimes_.at(i) = ELST_down_;
			}
		}

		auto t2 = atom(n_, irrTimes_, switchTimes_, angleDistances_, maxWindow_,  5.0, 0.5, 0.5);
		std::cout << t2 << "," <<
		          std::accumulate(irrTimes_.begin(), irrTimes_.end(), 0.0) + std::accumulate(switchTimes_.begin(), switchTimes_.end(), 0.0)
		          << std::endl;
		for (auto el : irrTimes_) {
			std::cout << el << ", ";
		}
		std::cout << "\n";
		for (auto el : switchTimes_) {
			std::cout << el << ", ";
		}
		std::cout << "\n";
		s += t2;
	}

	milliseconds ms2 = duration_cast<milliseconds>(
	                       system_clock::now().time_since_epoch());

	std::cout << "Time taken to run a for loop = " << ms2.count() - ms.count() << " milliseconds.\n";
}
