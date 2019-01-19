#pragma once
#include <string>
#include "Barrier.h"

struct Config
{
	double metricFactor;
	double timeFactor;

	size_t buffersNum;

	int width;
	int height;
	std::string logPath;

	// neuron
	double nucleusRadius;
	double somaRadius;
	double axonRadius;
	double axonLength;
	double axonHillockAreaFactor;
	double dendriteRadius;
	double dendriteLength;

	double particlesFlow[channel::TYPES_COUNT];

	// ions
	size_t particlesCount[particle::TYPES_COUNT];
	std::string particlesTextures[particle::TYPES_COUNT];

	// channels
	double NapIonsChannelsDensity[barrier::DENSITY_TYPES_COUNT];
	double KpIonsChannelsDensity[barrier::DENSITY_TYPES_COUNT];
	std::string channelsTextures[channel::TYPES_COUNT];
};