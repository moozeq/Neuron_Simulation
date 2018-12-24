#pragma once

namespace channel {
	enum Ion {
		NAP, KP, CLM
	};
	enum Gating {
		VOLTAGE_GATED, LIGAND_GATED
	};
	enum State {
		OPEN, CLOSED, INACTIVE, STATES_COUNT
	};
}

class Channel
{
	
public:
	Channel();
};

