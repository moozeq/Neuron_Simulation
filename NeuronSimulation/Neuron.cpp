#include "Neuron.h"

double Neuron::addBarrier(double x, double y, double z, double radius, double length, barrier::Type barrierType)
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
	case barrier::NUCLEUS:
		area = 4 * phy::pi * radius * radius;
		barriers.push_back(new Nucleus(midPoint, radius, lipidBilayerWidth));
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

		if (barrierIndex == barrier::DENDRITE_LOC) {
			Dendrite* dendrite = static_cast<Dendrite*>(barrier);
			float* synapse = dendrite->getSynapsePoint();
			if (outsideCoords[0] == synapse[0]) {
				channels[i] = Channel(insideCoords, outsideCoords, channel::NAP, channel::LIGAND_GATED);
				continue;
			}
		}
		else if (barrierIndex == barrier::DENDRITE_LOC) {

		}

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

		if (barrierIndex == barrier::DENDRITE_LOC) {
			Dendrite* dendrite = static_cast<Dendrite*>(barrier);
			float* synapse = dendrite->getSynapsePoint();
			if (outsideCoords[0] == synapse[0]) {
				--i;
				continue;
			}
		}

		channels[i] = Channel(insideCoords, outsideCoords, channel::KP, channel::CONST_OPEN);
	}
}

void Neuron::setupNucleus()
{
	areas[barrier::NUCLEUS_LOC] = addBarrier(0.0f, 0.0f, 0.0f, config.nucleusRadius, 0.0f, barrier::NUCLEUS);
}

void Neuron::setupSoma()
{
	areas[barrier::SOMA_LOC] = addBarrier(0.0f, 0.0f, 0.0f, config.somaRadius, 0.0f, barrier::SOMA);

	NapChannelsCount[barrier::SOMA] += areas[barrier::SOMA_LOC] * NapChannelsDensity[barrier::SOMA];
	KpChannelsCount[barrier::SOMA] += areas[barrier::SOMA_LOC] * KpChannelsDensity[barrier::SOMA];
}

void Neuron::setupAxon()
{
	float somaRadius = barriers[barrier::SOMA_LOC]->radius;
	float somaInnerRadius = barriers[barrier::SOMA_LOC]->innerRadius;
	float somaAxonGap = getGap(somaInnerRadius, config.axonRadius - lipidBilayerWidth / 2.0f);

	areas[barrier::AXON_LOC] = addBarrier(somaRadius + config.axonLength / 2.0f - somaAxonGap, 0.0f, 0.0f, config.axonRadius, config.axonLength, barrier::AXON);
	
	unsigned axonNapChannelsCount = areas[barrier::AXON_LOC] * ((1.0 - config.axonHillockAreaFactor) * NapChannelsDensity[barrier::AXON]);
	unsigned axonKpChannelsCount = areas[barrier::AXON_LOC] * ((1.0 - config.axonHillockAreaFactor) * KpChannelsDensity[barrier::AXON]);
	NapChannelsCount[barrier::AXON] += axonNapChannelsCount;
	KpChannelsCount[barrier::AXON] += axonKpChannelsCount;

	unsigned axonHillockNapChannelsCount = config.axonHillockAreaFactor * NapChannelsDensity[barrier::AXON_HILLOCK];
	unsigned axonHillockKpChannelsCount = config.axonHillockAreaFactor * KpChannelsDensity[barrier::AXON_HILLOCK];
	NapChannelsCount[barrier::AXON] += axonHillockNapChannelsCount;
	KpChannelsCount[barrier::AXON] += axonHillockKpChannelsCount;

	/*float synapseArea = phy::pi * config.axonRadius * config.axonRadius;
	unsigned axonSynapseNapChannelsCount = synapseArea * NapChannelsDensity[barrier::SYNAPSE];
	unsigned axonSynapseKpChannelsCount = synapseArea * KpChannelsDensity[barrier::SYNAPSE];
	NapChannelsCount[barrier::AXON] += axonSynapseNapChannelsCount;
	KpChannelsCount[barrier::AXON] += axonSynapseKpChannelsCount;*/

	float synapseProbability = 0;
	float axonHillockProbability = (float)axonHillockNapChannelsCount / (float)(axonNapChannelsCount + axonHillockNapChannelsCount);
	Axon* axon = static_cast<Axon*>(barriers[barrier::AXON_LOC]);
	axon->setSynapseProbability(synapseProbability);
	axon->setAxonHillockProbability(axonHillockProbability);
	axon->setAxonHillockAreaFactor(config.axonHillockAreaFactor);
}

void Neuron::setupDendrites()
{
	float somaRadius = barriers[barrier::SOMA_LOC]->radius;
	float somaInnerRadius = barriers[barrier::SOMA_LOC]->innerRadius;
	float dendriteSomaGap = getGap(somaInnerRadius, config.dendriteRadius - lipidBilayerWidth / 2.0f);

	areas[barrier::DENDRITE_LOC] = addBarrier(-somaRadius - config.dendriteLength / 2.0f + dendriteSomaGap, 0.0f, 0.0f, config.dendriteRadius, config.dendriteLength, barrier::DENDRITE);
	unsigned dendriteNapChannelsCount = areas[barrier::DENDRITE_LOC] * NapChannelsDensity[barrier::DENDRITE];
	unsigned dendriteKpChannelsCount = areas[barrier::DENDRITE_LOC] * KpChannelsDensity[barrier::DENDRITE];
	NapChannelsCount[barrier::DENDRITE] = dendriteNapChannelsCount;
	KpChannelsCount[barrier::DENDRITE] = dendriteKpChannelsCount;
	
	float synapseArea = phy::pi * config.dendriteRadius * config.dendriteRadius;
	unsigned dendriteSynapseNapChannelsCount = synapseArea * NapChannelsDensity[barrier::SYNAPSE];
	unsigned dendriteSynapseKpChannelsCount = synapseArea * KpChannelsDensity[barrier::SYNAPSE];
	NapChannelsCount[barrier::DENDRITE] += dendriteSynapseNapChannelsCount;
	KpChannelsCount[barrier::DENDRITE] += dendriteSynapseKpChannelsCount;

	float probability = (float)dendriteSynapseNapChannelsCount / (float)(dendriteSynapseNapChannelsCount + dendriteNapChannelsCount);
	Dendrite* dendrite = static_cast<Dendrite*>(barriers[barrier::DENDRITE_LOC]);
	dendrite->setSynapseProbability(probability);
}

void Neuron::setupChannelsIndexes()
{
	// nuclues
	barriers[barrier::NUCLEUS_LOC]->NapChannelsIndexFrom = 0;
	barriers[barrier::NUCLEUS_LOC]->NapChannelsIndexTo = 0;

	barriers[barrier::NUCLEUS_LOC]->KpChannelsIndexFrom = 0;
	barriers[barrier::NUCLEUS_LOC]->KpChannelsIndexTo = 0;

	// soma
	barriers[barrier::SOMA_LOC]->NapChannelsIndexFrom = 0;
	barriers[barrier::SOMA_LOC]->NapChannelsIndexTo = NapChannelsCount[barrier::SOMA];

	barriers[barrier::SOMA_LOC]->KpChannelsIndexFrom = allNapChannelsCount;
	barriers[barrier::SOMA_LOC]->KpChannelsIndexTo = allNapChannelsCount + KpChannelsCount[barrier::SOMA];

	// axon
	barriers[barrier::AXON_LOC]->NapChannelsIndexFrom = NapChannelsCount[barrier::SOMA];
	barriers[barrier::AXON_LOC]->NapChannelsIndexTo = NapChannelsCount[barrier::SOMA] + NapChannelsCount[barrier::AXON];

	barriers[barrier::AXON_LOC]->KpChannelsIndexFrom = allNapChannelsCount + KpChannelsCount[barrier::SOMA];
	barriers[barrier::AXON_LOC]->KpChannelsIndexTo = allNapChannelsCount + KpChannelsCount[barrier::SOMA] + KpChannelsCount[barrier::AXON];

	// dendrite
	barriers[barrier::DENDRITE_LOC]->NapChannelsIndexFrom = NapChannelsCount[barrier::SOMA] + NapChannelsCount[barrier::AXON];
	barriers[barrier::DENDRITE_LOC]->NapChannelsIndexTo = NapChannelsCount[barrier::SOMA] + NapChannelsCount[barrier::AXON] + NapChannelsCount[barrier::DENDRITE];

	barriers[barrier::DENDRITE_LOC]->KpChannelsIndexFrom = allNapChannelsCount + KpChannelsCount[barrier::SOMA] + KpChannelsCount[barrier::AXON];
	barriers[barrier::DENDRITE_LOC]->KpChannelsIndexTo = allNapChannelsCount + KpChannelsCount[barrier::SOMA] + KpChannelsCount[barrier::AXON] + KpChannelsCount[barrier::DENDRITE];
}

void Neuron::setupConnections()
{
	// add connection points
	Soma* soma = static_cast<Soma*>(barriers[barrier::SOMA_LOC]);
	float somaInnerRadius = barriers[barrier::SOMA_LOC]->innerRadius;
	float axonInnerRadius = barriers[barrier::AXON_LOC]->innerRadius;
	float dendriteInnerRadius = barriers[barrier::DENDRITE_LOC]->innerRadius;

	float somaAxonConnection[3] = { soma->x0 + somaInnerRadius, 0.0f, 0.0f };
	float dendriteSomaConnection[3] = { soma->x0 - somaInnerRadius, 0.0f, 0.0f };

	soma->addConnection(somaAxonConnection, axonInnerRadius, barrier::SOMA_AXON);
	soma->addConnection(dendriteSomaConnection, dendriteInnerRadius, barrier::DENDRITE_SOMA);
}

void Neuron::setupStructures()
{
	setupNucleus();
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

Neuron::Neuron(double _metricFactor, double _timeFactor, const Config& _config) :
	metricFactor(_metricFactor), timeFactor(_timeFactor),
	allNapChannelsCount(0), allKpChannelsCount(0)
{
	config = _config;

	lipidBilayerWidth = phy::lipidBilayerWidth / metricFactor;
	for (int i = 0; i < barrier::DENSITY_TYPES_COUNT; ++i) {
		NapChannelsDensity[i] = config.NapIonsChannelsDensity[i];
		KpChannelsDensity[i] = config.KpIonsChannelsDensity[i];
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

bool Neuron::checkCollision(Particle& nextParticleState, Particle& oldParticleState, const particle::Type type)
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
				Channel& currentChannel = channels[i];

				//// check if channel is open and if ion type is appropriate to channel type
				//if (currentChannel.state != channel::OPEN)
				//	continue;

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
					// particle passed through channel >>> outside layer <<< or when ligand gated stay in channel position
					else if (collisionType == collision::OUTSIDE || collisionType == collision::DISC_OUTSIDE) {

						// when neurotransmitter hits ligand gated channel
						if (type == particle::NEUROTRANSMITTER && currentChannel.gating == channel::LIGAND_GATED && currentChannel.state == channel::CLOSED) {
							nextParticleState.x = currentChannel.xOut;
							nextParticleState.y = currentChannel.yOut;
							nextParticleState.z = currentChannel.zOut;

							oldParticleState.vx = nextParticleState.vx = 0.0f;
							oldParticleState.vy = nextParticleState.vy = 0.0f;
							oldParticleState.vz = nextParticleState.vz = 0.0f;

							currentChannel.state = channel::OPEN;
						}
						else if (type == particle::NEUROTRANSMITTER && currentChannel.gating == channel::LIGAND_GATED && currentChannel.state == channel::OPEN) {
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
						else {
							nextParticleState.x = currentChannel.xIn;
							nextParticleState.y = currentChannel.yIn;
							nextParticleState.z = currentChannel.zIn;
						}
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
	Dendrite* dendrite = static_cast<Dendrite*>(barriers[barrier::DENDRITE_LOC]);
	position = dendrite->getSynapsePoint();
	return position;
}