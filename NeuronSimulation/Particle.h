#pragma once
#include "Utilities.h"

/*
To add new particle:
	par::ParticleType	new type
	newParticle()		case
	loadConfig			update particlesBufferSize
	setupStructures		loop with newParticle
	setupBuffers	
				texture (path in Config.h)
				texture (gluint in Simulation.h)
				VAO (var in Simulation.h)
				bind VAO to buffer
	render				bind texture and draw
*/

namespace particle {
	enum Type {
		NAP, KP, CLM, MASSIVEION, TYPES_COUNT
	};
}

namespace collision {
	enum Axis {
		X, Y, Z
	};
	enum Direction {
		NEG_X = 1, POS_X = 2, NEG_Y = 4, POS_Y = 8, NEG_Z = 16, POS_Z = 32
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
};