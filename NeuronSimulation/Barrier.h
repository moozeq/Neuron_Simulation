#pragma once
#include "Utilities.h"

namespace collision {
	enum Type {
		NONE, INSIDE, OUTSIDE, DISC
	};
}

namespace barrier {
	enum Type {
		SOMA, AXON, DENDRITE, AXON_HILLOCK, TYPES_COUNT = 3, DENSITY_TYPES_COUNT
	};
	enum Connection {
		SOMA_AXON, DENDRITE_SOMA
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

	float precision;
	float lipidBilayerWidth;
	float innerRadius;
	float outerRadius;
	float radius;
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
	
	virtual GLuint generateLayer(std::vector<FCoord8>& vertices, std::vector<UCoord3>& indices) = 0;

public:
	virtual collision::Type checkCollision(const float newCoords[3], const float oldCoords[3]) const = 0;
	virtual int getCollisionPoint(float* point, float newCoords[3], float oldCoords[3], collision::Type collisionType) const = 0;
	virtual bool getCollisionNormalVec(float collisionPoint[3], glm::vec3& n, collision::Type collisionType) const = 0;
	virtual bool getRandPointOnInnerLayer(float* point, glm::vec3& inOutVec) const = 0;
	void render() const;
};

