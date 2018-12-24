#pragma once
#include <string>
#include "Channel.h"

struct Config
{
	double metricFactorSq;
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

	float ionRadius;
	std::string NapIonTexturePath;
	std::string KpIonTexturePath;
	std::string ClmIonTexturePath;
	std::string otherParticlesTexturePath;

	// channels
	size_t NapIonsChannelsNum;
	size_t KpIonsChannelsNum;

	float channelRadius;
	std::string NapChannelTexturePath[channel::STATES_COUNT];
	std::string KpChannelTexturePath[channel::STATES_COUNT];
};