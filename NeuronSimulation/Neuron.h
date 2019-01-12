#pragma once
#include "Utilities.h"
#include "Barrier.h"
#include "Axon.h"
#include "Soma.h"
#include "Dendrite.h"

class Neuron
{
	friend class Simulation;

	ShaderProgram* barriersRenderProgram;

	std::vector<Channel> channels;
	std::vector<Barrier*> barriers;
	GLuint channelsPosBuf;

	double lipidBilayerWidth;
	double metricFactor;
	double timeFactor;

	double NapChannelsDensity[barrier::DENSITY_TYPES_COUNT];
	double KpChannelsDensity[barrier::DENSITY_TYPES_COUNT];

	unsigned NapChannelsCount[barrier::TYPES_COUNT];
	unsigned KpChannelsCount[barrier::TYPES_COUNT];

	unsigned allNapChannelsCount;
	unsigned allKpChannelsCount;

	float addBarrier(float x, float y, float z, float radius, float length, barrier::Type barrierType);

	void setupPrograms();
	void setupChannels(unsigned barrierIndex);
	void setupSoma();
	void setupAxon();
	void setupDendrites();
	void setupChannelsIndexes();
	void setupConnections();
	void setupStructures();

public:
	Neuron(double _metricFactor, double _timeFactor, double _NapChannelsDensity[barrier::DENSITY_TYPES_COUNT], double _KpChannelsDensity[barrier::DENSITY_TYPES_COUNT]);
	~Neuron();
	void render(shader::Uniforms uniforms) const;
	std::vector<float> getChannels();
	bool checkCollision(Particle& nextParticleState, Particle& oldParticleState, const particle::Type type) const;
	float* getSynapsePosition();
};