#pragma once
#include "Utilities.h"
#include "Config.h"
#include "Axon.h"
#include "Soma.h"
#include "Dendrite.h"
#include "Nucleus.h"

class Neuron
{
	friend class Simulation;

	Config config;
	ShaderProgram* barriersRenderProgram;

	std::vector<Channel> channels;
	std::vector<Barrier*> barriers;
	GLuint channelsPosBuf;

	double lipidBilayerWidth;
	double metricFactor;
	double timeFactor;

	double areas[barrier::PARTS_COUNT];

	double NapChannelsDensity[barrier::DENSITY_TYPES_COUNT];
	double KpChannelsDensity[barrier::DENSITY_TYPES_COUNT];

	unsigned NapChannelsCount[barrier::TYPES_COUNT];
	unsigned KpChannelsCount[barrier::TYPES_COUNT];

	unsigned allNapChannelsCount;
	unsigned allKpChannelsCount;

	double addBarrier(double x, double y, double z, double radius, double length, barrier::Type barrierType);

	void setupPrograms();
	void setupChannels(unsigned barrierIndex);
	void setupNucleus();
	void setupSoma();
	void setupAxon();
	void setupDendrites();
	void setupChannelsIndexes();
	void setupConnections();
	void setupStructures();

public:
	Neuron(double _metricFactor, double _timeFactor, const Config& config);
	~Neuron();
	void render(shader::Uniforms uniforms) const;
	std::vector<float> getChannels();
	bool checkCollision(Particle& nextParticleState, Particle& oldParticleState, const particle::Type type);
	float* getSynapsePosition();
};