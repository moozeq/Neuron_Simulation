#include "Channel.h"

Channel::Channel(float coordsIn[3], float coordsOut[3], channel::Type type, channel::Gating _gating) :
	xIn(coordsIn[0]), yIn(coordsIn[1]), zIn(coordsIn[2]), xOut(coordsOut[0]), yOut(coordsOut[1]), zOut(coordsOut[2]), gating(_gating), width(phy::lipidBilayerWidth), U(0.0f)
{
	switch (type)
	{
	case channel::NAP:
		radius = phy::NapR;
		break;
	case channel::KP:
		radius = phy::KpR;
		break;
	default:
		radius = 0.0f;
		break;
	}
}
