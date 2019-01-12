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
		barriers.push_back(new Dendrite(midPoint, radius, length, lipidBilayerWidth));
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

		while (!(barrier->getRandPointOnInnerLayer(insideCoords, inOutVec)));
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

		while (!(barrier->getRandPointOnInnerLayer(insideCoords, inOutVec)));
		glm::vec3 lipidBilayerLengthVec = glm::vec3(lipidBilayerWidth) * inOutVec;

		outsideCoords[0] = insideCoords[0] + lipidBilayerLengthVec[0];
		outsideCoords[1] = insideCoords[1] + lipidBilayerLengthVec[1];
		outsideCoords[2] = insideCoords[2] + lipidBilayerLengthVec[2];

		channels[i] = Channel(insideCoords, outsideCoords, channel::KP, channel::VOLTAGE_GATED);
	}
}

void Neuron::setupSoma()
{
	float somaRadius = 0.5f;
	float somaArea;
	unsigned somaNapChannelsCount = 0;
	unsigned somaKpChannelsCount = 0;

	somaArea = addBarrier(0.0f, 0.0f, 0.0f, somaRadius, 0.0f, barrier::SOMA);

	NapChannelsCount[barrier::SOMA] = somaArea * NapChannelsDensity[barrier::SOMA];
	KpChannelsCount[barrier::SOMA] = somaArea * KpChannelsDensity[barrier::SOMA];
}

void Neuron::setupAxon()
{
	// axon real radius = 0.125um, length = 5um, area = 2pi * 0.125um * 0.5um = ~4um2
	float axonRadius = 0.25f;
	float axonLength = 5.0f;
	float axonArea;
	double axonHillockAreaFactor = 0.1;

	float somaRadius = barriers[0]->radius;
	float somaInnerRadius = barriers[0]->innerRadius;
	float somaAxonGap = getGap(somaInnerRadius, axonRadius - lipidBilayerWidth / 2.0f);

	axonArea = addBarrier(somaRadius + axonLength / 2.0f - somaAxonGap, 0.0f, 0.0f, axonRadius, axonLength, barrier::AXON);
	unsigned axonNapChannelsCount = axonArea * ((1.0 - axonHillockAreaFactor) * NapChannelsDensity[barrier::AXON] + axonHillockAreaFactor * NapChannelsDensity[barrier::AXON_HILLOCK]);
	unsigned axonKpChannelsCount = axonArea * ((1.0 - axonHillockAreaFactor) * KpChannelsDensity[barrier::AXON] + axonHillockAreaFactor * KpChannelsDensity[barrier::AXON_HILLOCK]);
	NapChannelsCount[barrier::AXON] = axonNapChannelsCount;
	KpChannelsCount[barrier::AXON] = axonKpChannelsCount;

	float synapseArea = phy::pi * axonRadius * axonRadius;
	unsigned axonSynapseNapChannelsCount = synapseArea * NapChannelsDensity[barrier::SYNAPSE];
	unsigned axonSynapseKpChannelsCount = synapseArea * KpChannelsDensity[barrier::SYNAPSE];
	NapChannelsCount[barrier::AXON] += axonSynapseNapChannelsCount;
	KpChannelsCount[barrier::AXON] += axonSynapseKpChannelsCount;

	float probability = (float)axonSynapseNapChannelsCount / (float)(axonSynapseNapChannelsCount + axonNapChannelsCount);
	Axon* axon = static_cast<Axon*>(barriers[1]);
	axon->setSynapseProbability(probability);
}

void Neuron::setupDendrites()
{
	// dendrite real radius = 0.0625um, length = 1um
	float dendriteRadius = 0.125f;
	float dendriteLength = 1.0f;
	float dendriteArea;

	float somaRadius = barriers[0]->radius;
	float somaInnerRadius = barriers[0]->innerRadius;
	float dendriteSomaGap = getGap(somaInnerRadius, dendriteRadius - lipidBilayerWidth / 2.0f);

	dendriteArea = addBarrier(-somaRadius - dendriteLength / 2.0f + dendriteSomaGap, 0.0f, 0.0f, dendriteRadius, dendriteLength, barrier::DENDRITE);
	unsigned dendriteNapChannelsCount = dendriteArea * NapChannelsDensity[barrier::DENDRITE];
	unsigned dendriteKpChannelsCount = dendriteArea * KpChannelsDensity[barrier::DENDRITE];
	NapChannelsCount[barrier::DENDRITE] = dendriteNapChannelsCount;
	KpChannelsCount[barrier::DENDRITE] = dendriteKpChannelsCount;
	
	float synapseArea = phy::pi * dendriteRadius * dendriteRadius;
	unsigned dendriteSynapseNapChannelsCount = synapseArea * NapChannelsDensity[barrier::SYNAPSE];
	unsigned dendriteSynapseKpChannelsCount = synapseArea * KpChannelsDensity[barrier::SYNAPSE];
	NapChannelsCount[barrier::DENDRITE] += dendriteSynapseNapChannelsCount;
	KpChannelsCount[barrier::DENDRITE] += dendriteSynapseKpChannelsCount;

	float probability = (float)dendriteSynapseNapChannelsCount / (float)(dendriteSynapseNapChannelsCount + dendriteNapChannelsCount);
	Dendrite* dendrite = static_cast<Dendrite*>(barriers[2]);
	dendrite->setSynapseProbability(probability);
}

void Neuron::setupChannelsIndexes()
{
	// soma
	barriers[0]->NapChannelsIndexFrom = 0;
	barriers[0]->NapChannelsIndexTo = NapChannelsCount[barrier::SOMA];

	barriers[0]->KpChannelsIndexFrom = allNapChannelsCount;
	barriers[0]->KpChannelsIndexTo = allNapChannelsCount + KpChannelsCount[barrier::SOMA];

	// axon
	barriers[1]->NapChannelsIndexFrom = NapChannelsCount[barrier::SOMA];
	barriers[1]->NapChannelsIndexTo = NapChannelsCount[barrier::SOMA] + NapChannelsCount[barrier::AXON];

	barriers[1]->KpChannelsIndexFrom = allNapChannelsCount + KpChannelsCount[barrier::SOMA];
	barriers[1]->KpChannelsIndexTo = allNapChannelsCount + KpChannelsCount[barrier::SOMA] + KpChannelsCount[barrier::AXON];

	// dendrite
	barriers[2]->NapChannelsIndexFrom = NapChannelsCount[barrier::SOMA] + NapChannelsCount[barrier::AXON];
	barriers[2]->NapChannelsIndexTo = NapChannelsCount[barrier::SOMA] + NapChannelsCount[barrier::AXON] + NapChannelsCount[barrier::DENDRITE];

	barriers[2]->KpChannelsIndexFrom = allNapChannelsCount + KpChannelsCount[barrier::SOMA] + KpChannelsCount[barrier::AXON];
	barriers[2]->KpChannelsIndexTo = allNapChannelsCount + KpChannelsCount[barrier::SOMA] + KpChannelsCount[barrier::AXON] + KpChannelsCount[barrier::DENDRITE];

}

void Neuron::setupConnections()
{
	// add connection points
	Soma* soma = static_cast<Soma*>(barriers[0]);
	float somaInnerRadius = barriers[0]->innerRadius;
	float axonInnerRadius = barriers[1]->innerRadius;
	float dendriteInnerRadius = barriers[2]->innerRadius;

	float somaAxonConnection[3] = { soma->x0 + somaInnerRadius, 0.0f, 0.0f };
	float dendriteSomaConnection[3] = { soma->x0 - somaInnerRadius, 0.0f, 0.0f };

	soma->addConnection(somaAxonConnection, axonInnerRadius, barrier::SOMA_AXON);
	soma->addConnection(dendriteSomaConnection, dendriteInnerRadius, barrier::DENDRITE_SOMA);
}

void Neuron::setupStructures()
{
	setupSoma();
	setupAxon();
	setupDendrites();

	for (int i = 0; i < barrier::TYPES_COUNT; ++i) {
		allNapChannelsCount += NapChannelsCount[i];
		allKpChannelsCount += KpChannelsCount[i];
	}

	setupChannelsIndexes();
	setupConnections();

	// distribute channels on barriers
	channels.resize(allNapChannelsCount + allKpChannelsCount);
	for (unsigned i = 0; i < barriers.size(); ++i)
		setupChannels(i);
}

Neuron::Neuron(double _metricFactor, double _timeFactor, double _NapChannelsDensity[barrier::DENSITY_TYPES_COUNT], double _KpChannelsDensity[barrier::DENSITY_TYPES_COUNT]) :
	metricFactor(_metricFactor), timeFactor(_timeFactor),
	allNapChannelsCount(0), allKpChannelsCount(0)
{
	lipidBilayerWidth = phy::lipidBilayerWidth / metricFactor;
	for (int i = 0; i < barrier::DENSITY_TYPES_COUNT; ++i) {
		NapChannelsDensity[i] = _NapChannelsDensity[i];
		KpChannelsDensity[i] = _KpChannelsDensity[i];
	}
	for (int i = 0; i < barrier::TYPES_COUNT; ++i) {
		NapChannelsCount[i] = 0;
		KpChannelsCount[i] = 0;
	}
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
			barrier->getCollisionPoint(collisionPoint, newCoords, oldCoords, collisionType);

			// checking only channels with proper type
			unsigned channelsIndexFrom, channelsIndexTo;

			// Nap and neurotransmitters can only get to cell from outside
			if ((type == particle::NAP || type == particle::NEUROTRANSMITTER) && (collisionType == collision::OUTSIDE || collisionType == collision::DISC_OUTSIDE)) {
				channelsIndexFrom = barrier->NapChannelsIndexFrom;
				channelsIndexTo = barrier->NapChannelsIndexTo;
			}
			// Kp can travel through lipid bilayer in both directions
			else if (type == particle::KP) {
				channelsIndexFrom = barrier->KpChannelsIndexFrom;
				channelsIndexTo = barrier->KpChannelsIndexTo;
			}
			// other cases - simple reflection, no channels involved
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
				if (collisionType == collision::INSIDE || collisionType == collision::DISC_INSIDE) {
					dx = currentChannel.xIn - collisionPoint[0];
					dy = currentChannel.yIn - collisionPoint[1];
					dz = currentChannel.zIn - collisionPoint[2];
				}
				// coordsIn because its >>> outside layer <<<
				else if (collisionType == collision::OUTSIDE || collisionType == collision::DISC_OUTSIDE) {
					dx = currentChannel.xOut - collisionPoint[0];
					dy = currentChannel.yOut - collisionPoint[1];
					dz = currentChannel.zOut - collisionPoint[2];
				}
				else
					dx = dy = dz = 0.0f;
				d = metricFactor * sqrt(dx * dx + dy * dy + dz * dz);

				// collision point is within channel radius and channel is open
				if (d < currentChannel.radius) {
					//glm::vec3 channelVecTemp;

					// particle passed through channel >>> inside layer <<<
					if (collisionType == collision::INSIDE || collisionType == collision::DISC_INSIDE) {
						nextParticleState.x = currentChannel.xOut;
						nextParticleState.y = currentChannel.yOut;
						nextParticleState.z = currentChannel.zOut;
						// change direction of particle's velocity vector  >>> inside layer <<<
						//channelVecTemp = glm::vec3(currentChannel.xOut - currentChannel.xIn, currentChannel.yOut - currentChannel.yIn, currentChannel.zOut - currentChannel.zIn);
					}
					// particle passed through channel >>> outside layer <<<
					else if (collisionType == collision::OUTSIDE || collisionType == collision::DISC_OUTSIDE) {
						nextParticleState.x = currentChannel.xIn;
						nextParticleState.y = currentChannel.yIn;
						nextParticleState.z = currentChannel.zIn;
						// change direction of particle's velocity vector  >>> outside layer <<<
						//channelVecTemp = glm::vec3(currentChannel.xIn - currentChannel.xOut, currentChannel.yIn - currentChannel.yOut, currentChannel.zIn - currentChannel.zOut);
					}

					//glm::vec3 channelVec = glm::normalize(channelVecTemp);

					//// TODO detect if channelVec is inf/nan
					//if (glm::isinf(channelVec[0]) || glm::isinf(channelVec[1]) || glm::isinf(channelVec[2])) {
					//	std::cout << "\n current channel i = " + std::to_string(i) + " xin = " + std::to_string(currentChannel.xIn) + " yin = " + std::to_string(currentChannel.yIn) + " zin = " + std::to_string(currentChannel.zIn);
					//	std::cout << "\n current channel i = " + std::to_string(i) + " yout = " + std::to_string(currentChannel.yOut) + " zout = " + std::to_string(currentChannel.zOut);
					//	std::cout << "\n temp vec x = " + std::to_string(channelVecTemp[0]) + " y = " + std::to_string(channelVecTemp[1]) + " z = " + std::to_string(channelVecTemp[2]);
					//}

					/*float v = sqrt(nextParticleState.vx * nextParticleState.vx + nextParticleState.vy * nextParticleState.vy + nextParticleState.vz * nextParticleState.vz);
					oldParticleState.vx = nextParticleState.vx = v * channelVec[0];
					oldParticleState.vy = nextParticleState.vy = v * channelVec[1];
					oldParticleState.vz = nextParticleState.vz = v * channelVec[2];*/

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

float* Neuron::getSynapsePosition()
{
	float* position;
	Dendrite* dendrite = static_cast<Dendrite*>(barriers[2]);
	position = dendrite->getSynapsePoint();
	return position;
}