#pragma once
#include "Utilities.h"

/*
To add new particle:
	par::ParticleType	new type
	newParticle()		case
	setupStructures		loop with newParticle
	setupBuffers	
				texture (path in Config.h)
				texture (gluint in Simulation.h)
				VAO (var in Simulation.h)
				bind VAO to buffer
	render				bind texture and draw
*/

namespace particle {
	// must be in same order as in channels namespace
	enum Type {
		NONE = -1, NAP, KP, CLM, ORGANIC_ANION, NEUROTRANSMITTER, TYPES_COUNT
	};
}

struct Particle
{
	float x;
	float y;
	float z;
	float vx;
	float vy;
	float vz;

	float charge;
	float mass;
	float r;
};