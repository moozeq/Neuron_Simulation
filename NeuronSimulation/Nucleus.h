#pragma once
#include "Barrier.h"
#include "Utilities.h"

class Nucleus :
	public Barrier
{
	GLuint generateLayer(std::vector<FCoord8>& vertices, std::vector<UCoord3>& indices);

public:
	Nucleus(float coords[3], float _radius, float _lipidBilayerWidth);
	collision::Type checkCollision(const float newCoords[3], const float oldCoords[3]) const;
	int getCollisionPoint(float* point, float newCoords[3], float oldCoords[3], collision::Type collisionType) const;
	bool getCollisionNormalVec(float collisionPoint[3], glm::vec3& n, collision::Type collisionType) const;
	bool getRandPointOnInnerLayer(float* point, glm::vec3& inOutVec) const;
};