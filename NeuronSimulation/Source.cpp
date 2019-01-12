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
	// time 1.0 in simulation is 1.0 ns in real world
	config.timeFactor = 1e-7;

	config.buffersNum = 2;

	config.width = 1200;
	config.height = 800;
	config.logPath = "simulation.log";

	config.maxNeurotransmittersNum = 5000;
	config.NapIonsNum = 440;
	config.KpIonsNum = 4000;
	config.ClmIonsNum = 160;
	config.otherParticlesNum = 400;

	config.NapIonsChannelsDensity[barrier::SOMA] = 1000.0f;
	config.KpIonsChannelsDensity[barrier::SOMA] = 10.0f;

	config.NapIonsChannelsDensity[barrier::AXON] = 1000.0f;
	config.KpIonsChannelsDensity[barrier::AXON] = 10.0f;

	config.NapIonsChannelsDensity[barrier::DENDRITE] = 1000.0f;
	config.KpIonsChannelsDensity[barrier::DENDRITE] = 10.0f;

	config.NapIonsChannelsDensity[barrier::AXON_HILLOCK] = 3000.0f;
	config.KpIonsChannelsDensity[barrier::AXON_HILLOCK] = 0.0f;

	config.NapIonsChannelsDensity[barrier::SYNAPSE] = 5000.0f;
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