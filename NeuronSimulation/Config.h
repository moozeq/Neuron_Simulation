#pragma once
#include <string>
#include "Channel.h"

struct Config
{
	double metricFactor;
	double timeFactor;

	size_t buffersNum;

	int width;
	int height;
	std::string logPath;

	// ions
	size_t NapIonsNum;
	size_t KpIonsNum;
	size_t ClmIonsNum;
	size_t otherParticlesNum;

	std::string NapIonTexturePath;
	std::string KpIonTexturePath;
	std::string ClmIonTexturePath;
	std::string otherParticlesTexturePath;

	// channels
	size_t NapIonsChannelsNum;
	size_t KpIonsChannelsNum;

	std::string NapChannelTexturePath;
	std::string KpChannelTexturePath;
};