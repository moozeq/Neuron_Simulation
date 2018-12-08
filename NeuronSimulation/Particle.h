#pragma once
#include "Utilities.h"

/*
To add new particle:
	par::ParType	new type
	newParticle()	case
	loadConfig		update particlesBufferSize
	setupStructures	loop with newParticle
	setupBuffers	
			texture (path in Config.h)
			texture (gluint in Simulation.h)
			VAO (var in Simulation.h)
			bind VAO to buffer
	render			bind texture and draw
*/

namespace par {
	enum ParType {
		NAP, KP, CLM, MASSIVEION
	};
}

struct Particle
{
	double x;
	double y;
	double z;
	double vx;
	double vy;
	double vz;

	double charge;
	double mass;

	size_t index;
	
};