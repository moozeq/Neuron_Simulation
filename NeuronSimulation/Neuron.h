#pragma once
#include "Utilities.h"
#include "Barrier.h"
#include "Channel.h"

class Neuron
{
	friend class Simulation;

	ShaderProgram* barriersRenderProgram;

	std::vector<Channel> channels;

	GLuint channelsPosBuf;

	double lipidBilayerWidth;
	double metricFactor;
	double timeFactor;

	void addBarrier(float x, float y, float z, float radius, float width, float NapChannelsDensity, float KpChannelsDensity);

	void setupPrograms();
	void setupStructures();

public:
	Neuron(double _metricFactor, double _timeFactor);
	~Neuron();
	void render();
	std::vector<float> getChannels();
};

