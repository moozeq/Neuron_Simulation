#pragma once
class Barrier
{
	friend class Simulation;
	friend class Neuron;

	unsigned channelsIndexFrom;
	unsigned channelsIndexTo;

public:
	Barrier();
	bool checkCollision(float* coordinates);
};

