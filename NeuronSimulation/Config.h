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
	std::string logPath;

	double ionRadius;
	std::string NapIonTexturePath;
	std::string KpIonTexturePath;
	std::string ClmIonTexturePath;

	double channelRadius;
	std::string NapChannelTexturePath;
	std::string KpChannelTexturePath;
};