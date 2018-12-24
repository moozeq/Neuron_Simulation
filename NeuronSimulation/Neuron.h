#pragma once
#include "Utilities.h"
#include "Barrier.h"

class Neuron
{
	ShaderProgram* barriersRenderProgram;
	void setupPrograms();
	void setupStructures();
	void setupBuffers();

public:
	Neuron();
	void render();
};

