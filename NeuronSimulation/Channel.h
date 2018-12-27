#pragma once
#include "Utilities.h"

namespace channel {
	enum Type {
		NAP, KP, TYPES_COUNT
	};
	enum Gating {
		VOLTAGE_GATED, LIGAND_GATED
	};
	enum State {
		OPEN, CLOSED, INACTIVE, STATES_COUNT, NONE = 0
	};
}

class Channel
{
	friend class Simulation;
	friend class Neuron;

	channel::Type type;
	channel::Gating gating;
	channel::State state;

	float timeLeft;
	float U;
	float radius;
	float width;

	float xIn;
	float yIn;
	float zIn;

	float xOut;
	float yOut;
	float zOut;

public:
	Channel(float coordsIn[3], float coordsOut[3], channel::Type type, channel::Gating _gating);
};

