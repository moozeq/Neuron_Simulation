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

Neuron::Neuron()
{
	setupPrograms();
	float start[3] = { -1.0f, 0.0f, 0.0f };
	float stop[3] = { -1.1f, 0.0f, 0.0f };
	float start2[3] = { -0.5f, 0.0f, 0.0f };
	float stop2[3] = { -0.6f, 0.0f, 0.0f };
	float start3[3] = { -0.1f, 0.0f, 0.0f };
	float stop3[3] = { -0.0f, 0.0f, 0.0f };
	channels.push_back(Channel(start, stop, channel::NAP, channel::VOLTAGE_GATED));
	//channels.push_back(Channel(start2, stop2, channel::NAP, channel::VOLTAGE_GATED));
	//channels.push_back(Channel(start3, stop3, channel::NAP, channel::VOLTAGE_GATED));
}

Neuron::~Neuron()
{

}


void Neuron::render()
{

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
