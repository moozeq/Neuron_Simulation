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
	config.timeFactor = 1e-9;

	config.buffersNum = 2;

	config.width = 1200;
	config.height = 800;
	config.logPath = "simulation.log";

	config.NapIonsNum = 3000;
	config.KpIonsNum = 3000;
	config.ClmIonsNum = 3000;
	config.otherParticlesNum = 0;

	config.NapIonsChannelsDensity = 1000.0f;
	config.KpIonsChannelsDensity = 100.0f;

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