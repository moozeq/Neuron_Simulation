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
	config.metricFactor = 1e-9;
	config.timeFactor = 1e-3;

	config.width = 1200;
	config.height = 800;
	config.NapIonsNum = 1000;
	config.KpIonsNum = 1000;
	config.ClmIonsNum = 0;
	config.otherParticlesNum = 1;
	config.logPath = "simulation.log";

	config.ionRadius = 0.0078125;
	config.NapIonTexturePath = "NapIon.png";
	config.KpIonTexturePath = "KpIon.png";
	config.ClmIonTexturePath = "ClmIon.png";
	config.otherParticlesTexturePath = "otherParticles.png";

	config.channelRadius = 0.02;
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