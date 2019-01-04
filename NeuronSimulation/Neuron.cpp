#include "Neuron.h"

float Neuron::addBarrier(float x, float y, float z, float radius, float length)
{
	float coords[3] = { x, y, z };
	float area = 2 * phy::pi * radius * length;

	insideLayer.push_back(Barrier(coords, radius, length));
	outsideLayer.push_back(Barrier(coords, radius + lipidBilayerWidth, length));

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

void Neuron::setupStructures()
{
	float barrierRadius = 0.125f;
	float barrierLength = 2.0f;

	// barrier real radius = 0.125um, length = 5um, area = 2pi * 0.125um * 0.5um = ~4um2
	float area = addBarrier(0.0f, 0.0f, 0.0f, barrierRadius, barrierLength);
	unsigned NapChannelsCount = area * NapChannelsDensity;
	unsigned KpChannelsCount = area * KpChannelsDensity;

	outsideLayer[0].NapChannelsIndexFrom = insideLayer[0].NapChannelsIndexFrom = 0;
	outsideLayer[0].NapChannelsIndexTo = insideLayer[0].NapChannelsIndexTo = NapChannelsCount;

	outsideLayer[0].KpChannelsIndexFrom = insideLayer[0].KpChannelsIndexFrom = NapChannelsCount;
	outsideLayer[0].KpChannelsIndexTo = insideLayer[0].KpChannelsIndexTo = NapChannelsCount + KpChannelsCount;

	Barrier& barrier = insideLayer[0];

	// add Nap channels
	for (unsigned i = 0; i < NapChannelsCount; ++i) {
		float insideCoords[3];
		float outsideCoords[3];

		// draw x and angle
		float x = getRandDouble(barrier.startCoords[0], barrier.stopCoords[0]);
		float angle = getRandDouble(0.0, 2 * phy::pi);

		insideCoords[0] = x;
		insideCoords[1] = barrierRadius * sin(angle);
		insideCoords[2] = barrierRadius * cos(angle);

		glm::vec3 radiusVec = glm::normalize(glm::vec3(insideCoords[0] - x, insideCoords[1] - barrier.y0, insideCoords[2] - barrier.z0));
		glm::vec3 lipidBilayerLengthVec = glm::vec3(lipidBilayerWidth) * radiusVec;

		outsideCoords[0] = insideCoords[0] + lipidBilayerLengthVec[0];
		outsideCoords[1] = insideCoords[1] + lipidBilayerLengthVec[1];
		outsideCoords[2] = insideCoords[2] + lipidBilayerLengthVec[2];

		channels.push_back(Channel(insideCoords, outsideCoords, channel::NAP, channel::VOLTAGE_GATED));
	}
	
	// add Kp channels
	for (unsigned i = NapChannelsCount; i < NapChannelsCount + KpChannelsCount; ++i) {
		float insideCoords[3];
		float outsideCoords[3];

		// draw x and angle
		float x = getRandDouble(barrier.startCoords[0], barrier.stopCoords[0]);
		float angle = getRandDouble(0.0, 2 * phy::pi);

		insideCoords[0] = x;
		insideCoords[1] = barrierRadius * sin(angle);
		insideCoords[2] = barrierRadius * cos(angle);

		glm::vec3 radiusVec = glm::normalize(glm::vec3(insideCoords[0] - x, insideCoords[1] - barrier.y0, insideCoords[2] - barrier.z0));
		glm::vec3 lipidBilayerLengthVec = glm::vec3(lipidBilayerWidth) * radiusVec;

		outsideCoords[0] = insideCoords[0] + lipidBilayerLengthVec[0];
		outsideCoords[1] = insideCoords[1] + lipidBilayerLengthVec[1];
		outsideCoords[2] = insideCoords[2] + lipidBilayerLengthVec[2];

		channels.push_back(Channel(insideCoords, outsideCoords, channel::KP, channel::VOLTAGE_GATED));
	}
}

Neuron::Neuron(double _metricFactor, double _timeFactor, double _NapChannelsDensity, double _KpChannelsDensity) :
	metricFactor(_metricFactor), timeFactor(_timeFactor), NapChannelsDensity(_NapChannelsDensity), KpChannelsDensity(_KpChannelsDensity)
{
	lipidBilayerWidth = phy::lipidBilayerWidth / metricFactor;
	setupPrograms();
	setupStructures();
}

Neuron::~Neuron()
{
	delete barriersRenderProgram;
}

void Neuron::render(shader::Uniforms uniforms)
{
	barriersRenderProgram->use();
	barriersRenderProgram->setUniforms(uniforms);

	for (Barrier& barrier : outsideLayer)
		barrier.render();

	for (Barrier& barrier : insideLayer)
		barrier.render();
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

bool Neuron::checkCollision(Particle& nextParticleState, const Particle& oldParticleState, const particle::Type type)
{
	for (Barrier& barrier : outsideLayer) {
		float newCoords[3] = { nextParticleState.x, nextParticleState.y, nextParticleState.z };
		float oldCoords[3] = { oldParticleState.x, oldParticleState.y, oldParticleState.z };

		// barriers do not overlap, detect if collide with current barrier
		if (barrier.checkCollision(newCoords, oldCoords)) {
			float collisionPoint[3];
			barrier.getCollisionPoint(collisionPoint, newCoords, oldCoords);

			// checking only channels with proper type
			unsigned channelsIndexFrom, channelsIndexTo;
			if (type == particle::NAP) {
				channelsIndexFrom = barrier.NapChannelsIndexFrom;
				channelsIndexTo = barrier.NapChannelsIndexTo;
			}
			else if (type == particle::KP) {
				channelsIndexFrom = barrier.KpChannelsIndexFrom;
				channelsIndexTo = barrier.KpChannelsIndexTo;
			}
			else
				channelsIndexFrom = channelsIndexTo = 0;

			// check if collide with channel from current barrier
			for (unsigned i = channelsIndexFrom; i < channelsIndexTo; ++i) {
				const Channel& currentChannel = channels[i];

				// check if channel is open and if ion type is appropriate to channel type
				if (currentChannel.state != channel::OPEN)
					continue;

				// coordsOut because its >>> outside layer <<<
				float dx, dy, dz, d;
				dx = currentChannel.xOut - collisionPoint[0];
				dy = currentChannel.yOut - collisionPoint[1];
				dz = currentChannel.zOut - collisionPoint[2];
				d = metricFactor * sqrt(dx * dx + dy * dy + dz * dz);

				// collision point is within channel radius and channel is open
				if (d < currentChannel.radius) {

					// particle passed through channel >>> outside layer <<<
					nextParticleState.x = currentChannel.xIn;
					nextParticleState.y = currentChannel.yIn;
					nextParticleState.z = currentChannel.zIn;

					// change direction of particle's velocity vector  >>> outside layer <<<
					glm::vec3 channelVec = glm::vec3(currentChannel.xIn - currentChannel.xOut, currentChannel.yIn - currentChannel.yOut, currentChannel.zIn - currentChannel.zOut);
					channelVec = glm::normalize(channelVec);
					float v = sqrt(oldParticleState.vx * oldParticleState.vx + oldParticleState.vy * oldParticleState.vy + oldParticleState.vz * oldParticleState.vz);
					nextParticleState.vx = v * channelVec[0];
					nextParticleState.vy = v * channelVec[1];
					nextParticleState.vz = v * channelVec[2];

					// collide with channel
					return true;
				}
			}

			// TODO collide with barrier change coords and velocity
			float center[3] = { barrier.x0, barrier.y0, barrier.z0 };
			float h = getPointOnLineDistanceFromCenter(collisionPoint, center, barrier.radius);
			glm::vec3 n;
			if (nextParticleState.x <= barrier.x0)
				n = glm::normalize(glm::vec3(barrier.x0 - h - collisionPoint[0], collisionPoint[1], collisionPoint[2]));
			else
				n = glm::normalize(glm::vec3(barrier.x0 + h - collisionPoint[0], collisionPoint[1], collisionPoint[2]));
			glm::vec3 newVelocity = glm::reflect(glm::vec3(nextParticleState.vx, nextParticleState.vy, nextParticleState.vz), n);
			nextParticleState.x = oldCoords[0];
			nextParticleState.y = oldCoords[1];
			nextParticleState.z = oldCoords[2];
			
			nextParticleState.vx = newVelocity[0];
			nextParticleState.vy = newVelocity[1];
			nextParticleState.vz = newVelocity[2];
			return true;
		}
	}

	for (Barrier& barrier : insideLayer) {
		float newCoords[3] = { nextParticleState.x, nextParticleState.y, nextParticleState.z };
		float oldCoords[3] = { oldParticleState.x, oldParticleState.y, oldParticleState.z };

		// barriers do not overlap, detect if collide with current barrier
		if (barrier.checkCollision(newCoords, oldCoords)) {
			float collisionPoint[3];
			barrier.getCollisionPoint(collisionPoint, newCoords, oldCoords);

			// checking only channels with proper type
			unsigned channelsIndexFrom, channelsIndexTo;
			if (type == particle::NAP) {
				channelsIndexFrom = barrier.NapChannelsIndexFrom;
				channelsIndexTo = barrier.NapChannelsIndexTo;
			}
			else if (type == particle::KP) {
				channelsIndexFrom = barrier.KpChannelsIndexFrom;
				channelsIndexTo = barrier.KpChannelsIndexTo;
			}
			else
				channelsIndexFrom = channelsIndexTo = 0;

			// check if collide with channel from current barrier
			for (unsigned i = channelsIndexFrom; i < channelsIndexTo; ++i) {
				const Channel& currentChannel = channels[i];

				// check if channel is open and if ion type is appropriate to channel type
				if (currentChannel.state != channel::OPEN)
					continue;

				// coordsOut because its >>> inside layer <<<
				float dx, dy, dz, d;
				dx = currentChannel.xIn - collisionPoint[0];
				dy = currentChannel.yIn - collisionPoint[1];
				dz = currentChannel.zIn - collisionPoint[2];
				d = metricFactor * sqrt(dx * dx + dy * dy + dz * dz);

				// collision point is within channel radius and channel is open
				if (d < currentChannel.radius) {

					// particle passed through channel >>> inside layer <<<
					nextParticleState.x = currentChannel.xOut;
					nextParticleState.y = currentChannel.yOut;
					nextParticleState.z = currentChannel.zOut;

					// change direction of particle's velocity vector  >>> inside layer <<<
					glm::vec3 channelVec = glm::vec3(currentChannel.xOut - currentChannel.xIn, currentChannel.yOut - currentChannel.yIn, currentChannel.zOut - currentChannel.zIn);
					channelVec = glm::normalize(channelVec);
					float v = sqrt(oldParticleState.vx * oldParticleState.vx + oldParticleState.vy * oldParticleState.vy + oldParticleState.vz * oldParticleState.vz);
					nextParticleState.vx = v * channelVec[0];
					nextParticleState.vy = v * channelVec[1];
					nextParticleState.vz = v * channelVec[2];

					// collide with channel
					return true;
				}
			}

			// TODO collide with barrier change coords and velocity
			float center[3] = { barrier.x0, barrier.y0, barrier.z0 };
			float h = getPointOnLineDistanceFromCenter(collisionPoint, center, barrier.radius);
			glm::vec3 n;
			if (nextParticleState.x <= barrier.x0)
				n = glm::normalize(glm::vec3(barrier.x0 - h - collisionPoint[0], collisionPoint[1], collisionPoint[2]));
			else
				n = glm::normalize(glm::vec3(barrier.x0 + h - collisionPoint[0], collisionPoint[1], collisionPoint[2]));
			glm::vec3 newVelocity = glm::reflect(glm::vec3(nextParticleState.vx, nextParticleState.vy, nextParticleState.vz), n);
			nextParticleState.x = oldCoords[0];
			nextParticleState.y = oldCoords[1];
			nextParticleState.z = oldCoords[2];
			
			nextParticleState.vx = newVelocity[0];
			nextParticleState.vy = newVelocity[1];
			nextParticleState.vz = newVelocity[2];
			return true;
		}
	}

	return false;
}
