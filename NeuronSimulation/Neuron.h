#pragma once
#include "Utilities.h"
#include "Barrier.h"

class Neuron
{
	friend class Simulation;

	ShaderProgram* barriersRenderProgram;

	std::vector<Channel> channels;
	std::vector<Barrier> insideLayer;
	std::vector<Barrier> outsideLayer;
	GLuint channelsPosBuf;

	double lipidBilayerWidth;
	double metricFactor;
	double timeFactor;

	double NapChannelsDensity;
	double KpChannelsDensity;

	float addBarrier(float x, float y, float z, float radius, float length);

	void setupPrograms();
	void setupStructures();

public:
	Neuron(double _metricFactor, double _timeFactor, double _NapChannelsDensity, double _KpChannelsDensity);
	~Neuron();
	void render(shader::Uniforms uniforms);
	std::vector<float> getChannels();
	bool checkCollision(Particle& nextParticleState, const Particle& oldParticleState, const particle::Type type);
};