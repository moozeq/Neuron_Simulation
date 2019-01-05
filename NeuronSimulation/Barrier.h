#pragma once
#include "Utilities.h"
#include "Channel.h"

namespace collision {
	enum Type {
		NONE, INSIDE, OUTSIDE
	};
}

class Barrier
{
	friend class Simulation;
	friend class Neuron;

	unsigned NapChannelsIndexFrom;
	unsigned NapChannelsIndexTo;
	unsigned KpChannelsIndexFrom;
	unsigned KpChannelsIndexTo;

	unsigned slices;
	float radius;
	float length;
	float lipidBilayerWidth;

	float x0;
	float y0;
	float z0;

	float startCoords[3];
	float stopCoords[3];
	
	GLuint innerLayerVAO;
	GLuint outerLayerVAO;
	GLuint bilayerTexture;

	std::vector<FCoord8> innerLayerVertices;
	std::vector<UCoord3> innerLayerIndices;

	std::vector<FCoord8> outerLayerVertices;
	std::vector<UCoord3> outerLayerIndices;

	GLuint generateBilayer(std::vector<FCoord8>& vertices, std::vector<UCoord3>& indices);
	void addCircle(std::vector<FCoord8>& vertices);

public:
	Barrier(float coords[3], float _radius, float _length, float _lipidBilayerWidth);
	collision::Type checkCollision(float newCoords[3], float oldCoords[3]);
	bool getCollisionPoint(float* point, float newCoords[3], float oldCoords[3], collision::Type collisionType);
	glm::vec3 getCollisionNormalVec(float collisionPoint[3], collision::Type collisionType);
	void render();
};

