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

	void setupPrograms();
	void setupStructures();

public:
	Neuron();
	~Neuron();
	void render();
	std::vector<float> getChannels();
};

