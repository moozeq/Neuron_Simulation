#include "Neuron.h"

float Neuron::addBarrier(float x, float y, float z, float radius, float length, barrier::Type barrierType)
{
	float midPoint[3] = { x, y, z };
	float area;

	switch (barrierType) {
	case barrier::SOMA:
		area = 4 * phy::pi * radius * radius;
		barriers.push_back(new Soma(midPoint, radius, lipidBilayerWidth));
		break;
	case barrier::AXON:
		area = 2 * phy::pi * radius * length;
		barriers.push_back(new Axon(midPoint, radius, length, lipidBilayerWidth));
		break;
	case barrier::DENDRITE:
		area = 2 * phy::pi * radius * length;
		break;
	}

	// returns area in metric factor squared
	return area;
}

void Neuron::setupPrograms() {
	std::vector<const GLchar*> barriersPaths({
			   "Barriers.vert",
			   "Barriers.frag"
		});
	barriersRenderProgram = new ShaderProgram(shader::VF, barriersPaths);
}

void Neuron::setupChannels(unsigned barrierIndex)
{
	if (barrierIndex >= barriers.size())
		return;

	Barrier* barrier = barriers[barrierIndex];

	// add Nap channels
	for (unsigned i = barrier->NapChannelsIndexFrom; i < barrier->NapChannelsIndexTo; ++i) {
		float insideCoords[3];
		float outsideCoords[3];

		glm::vec3 inOutVec;

		barrier->getRandPointOnInnerLayer(insideCoords, inOutVec);
		glm::vec3 lipidBilayerLengthVec = glm::vec3(lipidBilayerWidth) * inOutVec;

		outsideCoords[0] = insideCoords[0] + lipidBilayerLengthVec[0];
		outsideCoords[1] = insideCoords[1] + lipidBilayerLengthVec[1];
		outsideCoords[2] = insideCoords[2] + lipidBilayerLengthVec[2];

		channels[i] = Channel(insideCoords, outsideCoords, channel::NAP, channel::VOLTAGE_GATED);
	}

	// add Kp channels
	for (unsigned i = barrier->KpChannelsIndexFrom; i < barrier->KpChannelsIndexTo; ++i) {
		float insideCoords[3];
		float outsideCoords[3];

		glm::vec3 inOutVec;

		barrier->getRandPointOnInnerLayer(insideCoords, inOutVec);
		glm::vec3 lipidBilayerLengthVec = glm::vec3(lipidBilayerWidth) * inOutVec;

		outsideCoords[0] = insideCoords[0] + lipidBilayerLengthVec[0];
		outsideCoords[1] = insideCoords[1] + lipidBilayerLengthVec[1];
		outsideCoords[2] = insideCoords[2] + lipidBilayerLengthVec[2];

		channels[i] = Channel(insideCoords, outsideCoords, channel::KP, channel::VOLTAGE_GATED);
	}
}

void Neuron::setupStructures()
{
	float somaRadius = 0.5f;
	float somaArea;
	unsigned somaNapChannelsCount = 0;
	unsigned somaKpChannelsCount = 0;

	// axon real radius = 0.125um, length = 5um, area = 2pi * 0.125um * 0.5um = ~4um2
	float axonRadius = 0.125f;
	float axonLength = 5.0f;
	float axonArea;
	unsigned axonNapChannelsCount = 0;
	unsigned axonKpChannelsCount = 0;

	somaArea = addBarrier(-axonLength / 2.0f - somaRadius, 0.0f, 0.0f, somaRadius, 0.0f, barrier::SOMA);
	axonArea = addBarrier(0.0f, 0.0f, 0.0f, axonRadius, axonLength, barrier::AXON);
	
	somaNapChannelsCount = somaArea * NapChannelsDensity;
	somaKpChannelsCount = somaArea * KpChannelsDensity;

	axonNapChannelsCount = axonArea * NapChannelsDensity;
	axonKpChannelsCount = axonArea * KpChannelsDensity;

	NapChannelsCount = somaNapChannelsCount + axonNapChannelsCount;
	KpChannelsCount = somaKpChannelsCount + axonKpChannelsCount;

	barriers[0]->NapChannelsIndexFrom = 0;
	barriers[0]->NapChannelsIndexTo = somaNapChannelsCount;

	barriers[0]->KpChannelsIndexFrom = NapChannelsCount;
	barriers[0]->KpChannelsIndexTo = NapChannelsCount + somaKpChannelsCount;


	barriers[1]->NapChannelsIndexFrom = somaNapChannelsCount;
	barriers[1]->NapChannelsIndexTo = somaNapChannelsCount + axonNapChannelsCount;

	barriers[1]->KpChannelsIndexFrom = NapChannelsCount + somaKpChannelsCount;
	barriers[1]->KpChannelsIndexTo = NapChannelsCount + somaKpChannelsCount + axonKpChannelsCount;
	
	channels.resize(NapChannelsCount + KpChannelsCount);

	for (unsigned i = 0; i < barriers.size(); ++i)
		setupChannels(i);
}

Neuron::Neuron(double _metricFactor, double _timeFactor, double _NapChannelsDensity, double _KpChannelsDensity) :
	metricFactor(_metricFactor), timeFactor(_timeFactor),
	NapChannelsDensity(_NapChannelsDensity), KpChannelsDensity(_KpChannelsDensity),
	NapChannelsCount(0), KpChannelsCount(0)
{
	lipidBilayerWidth = phy::lipidBilayerWidth / metricFactor;
	setupPrograms();
	setupStructures();
}

Neuron::~Neuron()
{
	delete barriersRenderProgram;
	for (Barrier* barrier : barriers)
		delete barrier;
}

void Neuron::render(shader::Uniforms uniforms) const
{
	barriersRenderProgram->use();
	barriersRenderProgram->setUniforms(uniforms);

	for (const Barrier* barrier : barriers)
		barrier->render();
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

bool Neuron::checkCollision(Particle& nextParticleState, Particle& oldParticleState, const particle::Type type) const
{
	for (const Barrier* barrier : barriers) {
		float newCoords[3] = { nextParticleState.x, nextParticleState.y, nextParticleState.z };
		float oldCoords[3] = { oldParticleState.x, oldParticleState.y, oldParticleState.z };

		// barriers do not overlap, detect if collide with current barrier
		collision::Type collisionType = barrier->checkCollision(newCoords, oldCoords);

		if (collisionType) {
			float collisionPoint[3];
			bool collisionStrike = barrier->getCollisionPoint(collisionPoint, newCoords, oldCoords, collisionType);
			if (!collisionStrike)
				return false;

			// checking only channels with proper type
			unsigned channelsIndexFrom, channelsIndexTo;
			if (type == particle::NAP) {
				channelsIndexFrom = barrier->NapChannelsIndexFrom;
				channelsIndexTo = barrier->NapChannelsIndexTo;
			}
			else if (type == particle::KP) {
				channelsIndexFrom = barrier->KpChannelsIndexFrom;
				channelsIndexTo = barrier->KpChannelsIndexTo;
			}
			else
				channelsIndexFrom = channelsIndexTo = 0;

			// check if collide with channel from current barrier
			for (unsigned i = channelsIndexFrom; i < channelsIndexTo; ++i) {
				const Channel& currentChannel = channels[i];

				// check if channel is open and if ion type is appropriate to channel type
				if (currentChannel.state != channel::OPEN)
					continue;

				float dx, dy, dz, d;

				// coordsIn because its >>> inside layer <<<
				if (collisionType == collision::INSIDE) {
					dx = currentChannel.xIn - collisionPoint[0];
					dy = currentChannel.yIn - collisionPoint[1];
					dz = currentChannel.zIn - collisionPoint[2];
				}
				// coordsIn because its >>> outside layer <<<
				else if (collisionType == collision::OUTSIDE) {
					dx = currentChannel.xOut - collisionPoint[0];
					dy = currentChannel.yOut - collisionPoint[1];
					dz = currentChannel.zOut - collisionPoint[2];
				}
				d = metricFactor * sqrt(dx * dx + dy * dy + dz * dz);

				// collision point is within channel radius and channel is open
				if (d < currentChannel.radius) {
					glm::vec3 channelVec;

					// particle passed through channel >>> inside layer <<<
					if (collisionType == collision::INSIDE) {
						nextParticleState.x = currentChannel.xOut;
						nextParticleState.y = currentChannel.yOut;
						nextParticleState.z = currentChannel.zOut;
						// change direction of particle's velocity vector  >>> inside layer <<<
						channelVec = glm::vec3(currentChannel.xOut - currentChannel.xIn, currentChannel.yOut - currentChannel.yIn, currentChannel.zOut - currentChannel.zIn);
					}
					// particle passed through channel >>> outside layer <<<
					else if (collisionType == collision::OUTSIDE) {
						nextParticleState.x = currentChannel.xIn;
						nextParticleState.y = currentChannel.yIn;
						nextParticleState.z = currentChannel.zIn;
						// change direction of particle's velocity vector  >>> outside layer <<<
						channelVec = glm::vec3(currentChannel.xIn - currentChannel.xOut, currentChannel.yIn - currentChannel.yOut, currentChannel.zIn - currentChannel.zOut);
					}

					channelVec = glm::normalize(channelVec);
					float v = sqrt(nextParticleState.vx * nextParticleState.vx + nextParticleState.vy * nextParticleState.vy + nextParticleState.vz * nextParticleState.vz);
					oldParticleState.vx = nextParticleState.vx = v * channelVec[0];
					oldParticleState.vy = nextParticleState.vy = v * channelVec[1];
					oldParticleState.vz = nextParticleState.vz = v * channelVec[2];

					// collide with channel
					return true;
				}
			}

			// TODO collide with barrier change coords and velocity
			glm::vec3 n;
			barrier->getCollisionNormalVec(collisionPoint, n, collisionType);
			glm::vec3 newVelocity = glm::reflect(glm::vec3(nextParticleState.vx, nextParticleState.vy, nextParticleState.vz), n);
			nextParticleState.x = collisionPoint[0];
			nextParticleState.y = collisionPoint[1];
			nextParticleState.z = collisionPoint[2];

			// need to update velocities in both particle states
			oldParticleState.vx = nextParticleState.vx = newVelocity[0];
			oldParticleState.vy = nextParticleState.vy = newVelocity[1];
			oldParticleState.vz = nextParticleState.vz = newVelocity[2];
			return true;
		}
	}

	return false;
}
