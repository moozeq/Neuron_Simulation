#pragma once
#include <string>
#include "Barrier.h"

struct Config
{
	double metricFactor;
	double timeFactor;

	size_t buffersNum;
	double NapInflow;

	int width;
	int height;
	std::string logPath;

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