#include "Neuron.h"

void Neuron::addBarrier(float x, float y, float z, float radius, float length, float NapChannelsDensity, float KpChannelsDensity)
{
	float coords[3] = { x, y, z };
	insideLayer.push_back(Barrier(coords, radius, length));
	outsideLayer.push_back(Barrier(coords, radius + lipidBilayerWidth, length));
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
	addBarrier(0.0f, 0.0f, 0.0f, 0.125f, 0.25f, 10.0f, 0.0f);
	addBarrier(0.25f, 0.0f, 0.0f, 0.25f, 0.25f, 10.0f, 0.0f);


	float start[3] = { 0.0f, 0.0f, 0.125f };
	float stop[3] = { 0.0f, 0.0f, 0.125f + lipidBilayerWidth };
	float start2[3] = { -0.5f, 0.0f, 0.0f };
	float stop2[3] = { -0.6f, 0.0f, 0.0f };
	float start3[3] = { -0.1f, 0.0f, 0.0f };
	float stop3[3] = { -0.0f, 0.0f, 0.0f };
	channels.push_back(Channel(start, stop, channel::NAP, channel::VOLTAGE_GATED));
	//channels.push_back(Channel(start2, stop2, channel::NAP, channel::VOLTAGE_GATED));
	//channels.push_back(Channel(start3, stop3, channel::NAP, channel::VOLTAGE_GATED));
}

Neuron::Neuron(double _metricFactor, double _timeFactor) :
	metricFactor(_metricFactor), timeFactor(_timeFactor)
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

			// check if collide with channel from current barrier
			for (unsigned i = barrier.channelsIndexFrom; i < barrier.channelsIndexTo; ++i) {
				const Channel& currentChannel = channels[i];

				// check if channel is open and if ion type is appropriate to channel type
				if (currentChannel.state != channel::OPEN || currentChannel.type != type)
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

			// check if collide with channel from current barrier
			for (unsigned i = barrier.channelsIndexFrom; i < barrier.channelsIndexTo; ++i) {
				const Channel& currentChannel = channels[i];

				// check if channel is open and if ion type is appropriate to channel type
				if (currentChannel.state != channel::OPEN || currentChannel.type != type)
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
			return true;
		}
	}

	return false;
}
