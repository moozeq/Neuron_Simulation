#pragma once
#include "Utilities.h"

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