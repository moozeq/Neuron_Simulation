#pragma once
#include "Barrier.h"
#include "Utilities.h"

class Axon :
	public Barrier
{
	unsigned slices;
	float length;
	float discPoint[3];
	float synapseProbability;

	float startCoords[3];
	float stopCoords[3];

	GLuint generateLayer(std::vector<FCoord8>& vertices, std::vector<UCoord3>& indices);
	void addCircle(std::vector<FCoord8>& vertices);
	void addDisc(std::vector<UCoord3>& indices, unsigned midPointIndex);

public:
	Axon(float coords[3], float _radius, float _length, float _lipidBilayerWidth);
	collision::Type checkCollision(const float newCoords[3], const float oldCoords[3]) const;
	int getCollisionPoint(float* point, float newCoords[3], float oldCoords[3], collision::Type collisionType) const;
	bool getCollisionNormalVec(float collisionPoint[3], glm::vec3& n, collision::Type collisionType) const;
	bool getRandPointOnInnerLayer(float* point, glm::vec3& inOutVec) const;

	void setSynapseProbability(float probability);
};