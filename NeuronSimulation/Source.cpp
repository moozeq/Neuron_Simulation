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
	// time 1.0 in simulation is 1.0 ps in real world
	config.timeFactor = 1e-12;

	config.buffersNum = 2;

	config.width = 1200;
	config.height = 800;
	config.logPath = "simulation.log";

	config.NapIonsNum = 100;
	config.KpIonsNum = 100;
	config.ClmIonsNum = 200;
	config.otherParticlesNum = 0;

	config.NapIonsChannelsNum = 1;
	config.KpIonsChannelsNum = 0;

	config.NapIonTexturePath = "NapIon.png";
	config.KpIonTexturePath = "KpIon.png";
	config.ClmIonTexturePath = "ClmIon.png";
	config.otherParticlesTexturePath = "otherParticles.png";

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