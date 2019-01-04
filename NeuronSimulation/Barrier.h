#pragma once
#include "Utilities.h"
#include "Channel.h"

class Barrier
{
	friend class Simulation;
	friend class Neuron;

	unsigned channelsIndexFrom;
	unsigned channelsIndexTo;

	unsigned slices;
	float radius;
	float length;

	float x0;
	float y0;
	float z0;

	float startCoords[3];
	float stopCoords[3];
	
	GLuint circleVAO;
	GLuint layerVAO;

	GLuint circlesTexture;
	GLuint bilayerTexture;

	std::vector<FCoord8> circlesVertices;
	std::vector<FCoord8> bilayerVertices;
	std::vector<UCoord3> circlesIndices;
	std::vector<UCoord3> bilayerIndices;

	GLuint generateCircles();
	GLuint generateBilayer();
	void addCircle();

public:
	Barrier(float coords[3], float _radius, float _length);
	bool checkCollision(float newCoords[3], float oldCoords[3]);
	bool getCollisionPoint(float* point, float newCoords[3], float oldCoords[3]);
	void render();
};

