#pragma once
#include <string>

struct Config
{
	double metricFactor;
	double timeFactor;

	int width;
	int height;
	size_t NapIonsNum;
	size_t KpIonsNum;
	size_t ClmIonsNum;
	size_t otherParticlesNum;
	std::string logPath;

	float ionRadius;
	std::string NapIonTexturePath;
	std::string KpIonTexturePath;
	std::string ClmIonTexturePath;
	std::string otherParticlesTexturePath;

	float channelRadius;
	std::string NapChannelTexturePath;
	std::string KpChannelTexturePath;
};