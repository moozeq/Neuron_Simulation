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
	std::string logPath;

	double ionRadius;
	std::string NapIonTexturePath;
	std::string KpIonTexturePath;

	double channelRadius;
	std::string NapChannelTexturePath;
	std::string KpChannelTexturePath;
};