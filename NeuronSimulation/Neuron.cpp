#include "Neuron.h"

void Neuron::setupPrograms() {
	std::vector<const GLchar*> barriersPaths({
			   "Barriers.vert",
			   "Barriers.frag"
		});
	barriersRenderProgram = new ShaderProgram(shader::VF, barriersPaths);
}

void Neuron::setupStructures()
{
}

void Neuron::setupBuffers()
{
}

Neuron::Neuron()
{
	setupPrograms();
}


void Neuron::render()
{

}