#pragma once
#include "Utilities.h"

namespace collision {
	enum Type {
		NONE, INSIDE, OUTSIDE
	};
}

namespace barrier {
	enum Type {
		SOMA, AXON, DENDRITE
	};
}

class Barrier
{
	friend class Simulation;
	friend class Neuron;

protected:
	Barrier();

	unsigned NapChannelsIndexFrom;
	unsigned NapChannelsIndexTo;
	unsigned KpChannelsIndexFrom;
	unsigned KpChannelsIndexTo;

	float x0;
	float y0;
	float z0;

	GLuint innerLayerVAO;
	GLuint outerLayerVAO;
	GLuint bilayerTexture;

	std::vector<FCoord8> innerLayerVertices;
	std::vector<UCoord3> innerLayerIndices;

	std::vector<FCoord8> outerLayerVertices;
	std::vector<UCoord3> outerLayerIndices;

public:
	virtual collision::Type checkCollision(const float newCoords[3], const float oldCoords[3]) const = 0;
	virtual bool getCollisionPoint(float* point, float newCoords[3], float oldCoords[3], collision::Type collisionType) const = 0;
	virtual bool getCollisionNormalVec(float collisionPoint[3], glm::vec3& n, collision::Type collisionType) const = 0;
	virtual bool getRandPointOnInnerLayer(float* point, glm::vec3& inOutVec) const = 0;
	void render() const;
};

