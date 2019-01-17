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
	// distance 1.0 in simulation is 1.0 um in real world
	config.metricFactor = 1e-6;
	// time 1.0 in simulation is 100.0 ns in real world
	config.timeFactor = 1e-7;

	config.buffersNum = 2;

	// nucleus
	config.nucleusRadius = 0.125;
	// soma
	config.somaRadius = 0.4;
	// axon
	config.axonRadius = 0.25;
	config.axonLength = 12.0;
	config.axonHillockAreaFactor = 0.03125;
	// dendrite
	config.dendriteRadius = 0.0625;
	config.dendriteLength = 0.5;

	config.NapInflow = 1.0;

	config.width = 1200;
	config.height = 800;
	config.logPath = "simulation.log";

	config.maxNeurotransmittersNum = 4000;
	config.NapIonsNum = 16000;
	config.KpIonsNum = 0;
	config.ClmIonsNum = 0;
	config.otherParticlesNum = 0;

	config.NapIonsChannelsDensity[barrier::SOMA] = 000.0f;
	config.KpIonsChannelsDensity[barrier::SOMA] = 00.0f;

	config.NapIonsChannelsDensity[barrier::AXON] = 1000.0f;
	config.KpIonsChannelsDensity[barrier::AXON] = 100.0f;

	config.NapIonsChannelsDensity[barrier::DENDRITE] = 000.0f;
	config.KpIonsChannelsDensity[barrier::DENDRITE] = 00.0f;

	config.NapIonsChannelsDensity[barrier::AXON_HILLOCK] = 64000.0f;
	config.KpIonsChannelsDensity[barrier::AXON_HILLOCK] = 0.0f;

	config.NapIonsChannelsDensity[barrier::SYNAPSE] = 15000.0f;
	config.KpIonsChannelsDensity[barrier::SYNAPSE] = 0.0f;

	config.NapIonTexturePath = "NapIon.png";
	config.KpIonTexturePath = "KpIon.png";
	config.ClmIonTexturePath = "ClmIon.png";
	config.otherParticlesTexturePath = "otherParticles.png";
	config.neurotransmittersTexturePath = "neurotransmitters.png";

	config.NapChannelTexturePath = "NapChannel.png";
	config.KpChannelTexturePath = "KpChannel.png";

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