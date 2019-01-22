#define STB_IMAGE_IMPLEMENTATION
#include "STBImage.h"

#include "Simulation.h"
#include "Config.h"
#include <Windows.h>
extern "C" {
	_declspec(dllexport) DWORD NvOptimusEnablement = 0x00000001;
}

int main(void)
{
	Config config;

	config.width = 1200;
	config.height = 800;

	// distance 1.0 in simulation is 1.0 um in real world
	config.metricFactor = 1e-6;
	// time 1.0 in simulation is 100.0 ns in real world
	config.timeFactor = 1e-7;

	// nucleus
	config.nucleusRadius = 0.125;
	// soma
	config.somaRadius = 0.4;
	// axon
	config.axonRadius = 0.25;
	config.axonLength = 6.0;
	config.axonHillockAreaFactor = 0.0625;
	// dendrite
	config.dendriteRadius = 0.0625;
	config.dendriteLength = 0.5;

	config.particlesFlow[channel::NAP] = 1e6;
	config.particlesFlow[channel::KP] = 4e6;

	config.particlesCount[particle::NAP] = 10000;
	config.particlesCount[particle::KP] = 10000;
	config.particlesCount[particle::CLM] = 0;
	config.particlesCount[particle::ORGANIC_ANION] = 0;
	config.particlesCount[particle::NEUROTRANSMITTER] = 1000;

	config.NapIonsChannelsDensity[barrier::SOMA] = 0.0f;
	config.KpIonsChannelsDensity[barrier::SOMA] = 0.0f;

	config.NapIonsChannelsDensity[barrier::AXON] = 180.0f;
	config.KpIonsChannelsDensity[barrier::AXON] = 180.0f;

	config.NapIonsChannelsDensity[barrier::DENDRITE] = 0.0f;
	config.KpIonsChannelsDensity[barrier::DENDRITE] = 0.0f;

	config.NapIonsChannelsDensity[barrier::AXON_HILLOCK] = 16000.0f;
	config.KpIonsChannelsDensity[barrier::AXON_HILLOCK] = 16000.0f;

	config.NapIonsChannelsDensity[barrier::SYNAPSE] = 16000.0f;
	config.KpIonsChannelsDensity[barrier::SYNAPSE] = 0.0f;

	config.logPath = "simulation.log";

	config.particlesTextures[particle::NAP] = "NapIon.png";
	config.particlesTextures[particle::KP] = "KpIon.png";
	config.particlesTextures[particle::CLM] = "ClmIon.png";
	config.particlesTextures[particle::ORGANIC_ANION] = "otherParticles.png";
	config.particlesTextures[particle::NEUROTRANSMITTER] = "neurotransmitters.png";

	config.channelsTextures[channel::NAP] = "NapChannel.png";
	config.channelsTextures[channel::KP] = "KpChannel.png";

	try {
		Simulation simulation(config);
		simulation.start();
	}
	catch (std::exception& e) {
		std::ofstream logfile;
		logfile.open(config.logPath, std::ofstream::app | std::ofstream::binary);
		if (logfile.good())
			log(logfile, e.what());
		return 1;
	}

	return 0;
}