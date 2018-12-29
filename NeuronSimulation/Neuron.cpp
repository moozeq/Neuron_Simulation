#include "Neuron.h"

void Neuron::addBarrier(float x, float y, float z, float radius, float length, float NapChannelsDensity, float KpChannelsDensity)
{
	float coords[3] = { x, y, z };
	insideLayer.push_back(Barrier(coords, radius, length));
	outsideLayer.push_back(Barrier(coords, radius + lipidBilayerWidth, length));
}

void Neuron::setupPrograms() {
	std::vector<const GLchar*> barriersPaths({
			   "Barriers.vert",
			   "Barriers.frag"
		});
	barriersRenderProgram = new ShaderProgram(shader::VF, barriersPaths);
}

void Neuron::setupStructures()
{
	addBarrier(0.0f, 0.0f, 0.0f, 0.125f, 0.25f, 10.0f, 0.0f);
	addBarrier(0.25f, 0.0f, 0.0f, 0.25f, 0.25f, 10.0f, 0.0f);


	float start[3] = { 0.0f, 0.0f, 0.0f };
	float stop[3] = { 0.0f - lipidBilayerWidth, 0.0f, 0.0f };
	float start2[3] = { -0.5f, 0.0f, 0.0f };
	float stop2[3] = { -0.6f, 0.0f, 0.0f };
	float start3[3] = { -0.1f, 0.0f, 0.0f };
	float stop3[3] = { -0.0f, 0.0f, 0.0f };
	channels.push_back(Channel(start, stop, channel::NAP, channel::VOLTAGE_GATED));
	//channels.push_back(Channel(start2, stop2, channel::NAP, channel::VOLTAGE_GATED));
	//channels.push_back(Channel(start3, stop3, channel::NAP, channel::VOLTAGE_GATED));
}

Neuron::Neuron(double _metricFactor, double _timeFactor) :
	metricFactor(_metricFactor), timeFactor(_timeFactor)
{
	lipidBilayerWidth = phy::lipidBilayerWidth / metricFactor;
	setupPrograms();
	setupStructures();
}

Neuron::~Neuron()
{
	delete barriersRenderProgram;
}


void Neuron::render(shader::Uniforms uniforms)
{
	barriersRenderProgram->use();
	barriersRenderProgram->setUniforms(uniforms);

	for (Barrier& barrier : outsideLayer)
		barrier.render();

	for (Barrier& barrier : insideLayer)
		barrier.render();
}

std::vector<float> Neuron::getChannels()
{
	std::vector<float> channelsAttribs;
	channelsAttribs.reserve(channels.size() * 4);
	for (Channel& channel : channels) {
		channelsAttribs.push_back((channel.xIn + channel.xOut) / 2);
		channelsAttribs.push_back((channel.yIn + channel.yOut) / 2);
		channelsAttribs.push_back((channel.zIn + channel.zOut) / 2);
		channelsAttribs.push_back(0.0f);
	}
	return channelsAttribs;
}
