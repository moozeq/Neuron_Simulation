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

	double NapInflow;

	// ions
	size_t NapIonsNum;
	size_t KpIonsNum;
	size_t ClmIonsNum;
	size_t otherParticlesNum;
	size_t maxNeurotransmittersNum;

	std::string NapIonTexturePath;
	std::string KpIonTexturePath;
	std::string ClmIonTexturePath;
	std::string otherParticlesTexturePath;
	std::string neurotransmittersTexturePath;

	// channels
	double NapIonsChannelsDensity[barrier::DENSITY_TYPES_COUNT];
	double KpIonsChannelsDensity[barrier::DENSITY_TYPES_COUNT];

	std::string NapChannelTexturePath;
	std::string KpChannelTexturePath;
};