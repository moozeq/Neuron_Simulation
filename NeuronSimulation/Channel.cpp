#include "Channel.h"

// need for vector.resize()
Channel::Channel()
{
}

Channel::Channel(float coordsIn[3], float coordsOut[3], channel::Type _type, channel::Gating _gating) :
	xIn(coordsIn[0]), yIn(coordsIn[1]), zIn(coordsIn[2]),
	xOut(coordsOut[0]), yOut(coordsOut[1]), zOut(coordsOut[2]),
	type(_type), gating(_gating), state(channel::CLOSED),
	width(phy::lipidBilayerWidth), U(0.0f), timeLeft(0.0f)
{
	switch (_type)
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
