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
	config.metricFactorSq = 1e-9;
	config.timeFactor = 1e-9;

	config.buffersNum = 2;

	config.width = 1200;
	config.height = 800;
	config.logPath = "simulation.log";

	config.NapIonsNum = 200;
	config.KpIonsNum = 0;
	config.ClmIonsNum = 0;
	config.otherParticlesNum = 0;

	config.NapIonsChannelsNum = 1;
	config.KpIonsChannelsNum = 0;

	config.ionRadius = 0.0078125;
	config.NapIonTexturePath = "NapIon.png";
	config.KpIonTexturePath = "KpIon.png";
	config.ClmIonTexturePath = "ClmIon.png";
	config.otherParticlesTexturePath = "otherParticles.png";

	config.channelRadius = 0.015625;
	config.NapChannelTexturePath[channel::OPEN] = "NapChannelOpen.png";
	config.NapChannelTexturePath[channel::CLOSED] = "NapChannelClosed.png";
	config.NapChannelTexturePath[channel::INACTIVE] = "NapChannelInactive.png";
	config.KpChannelTexturePath[channel::OPEN] = "KpChannelOpen.png";
	config.KpChannelTexturePath[channel::CLOSED] = "KpChannelClosed.png";
	config.KpChannelTexturePath[channel::INACTIVE] = "KpChannelInactive.png";

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