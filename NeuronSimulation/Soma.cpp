#include "Soma.h"

Soma::Soma(float coords[3], float _radius, float _lipidBilayerWidth) :
	radius(_radius), lipidBilayerWidth(_lipidBilayerWidth)
{
	x0 = coords[0];
	y0 = coords[1];
	z0 = coords[2];

	radius -= lipidBilayerWidth / 2;
	//innerLayerVAO = generateBilayer(innerLayerVertices, innerLayerIndices);

	radius += lipidBilayerWidth;
	//outerLayerVAO = generateBilayer(outerLayerVertices, outerLayerIndices);

	radius -= lipidBilayerWidth / 2;
}

collision::Type Soma::checkCollision(const float newCoords[3], const float oldCoords[3]) const
{
	//// not in barrier length
	//if (!(newCoords[0] <= stopCoords[0] && newCoords[0] >= startCoords[0] &&
	//	oldCoords[0] <= stopCoords[0] && oldCoords[0] >= startCoords[0]))
	//	return collision::NONE;

	//float oldD = getPointLineDistance(oldCoords, startCoords, stopCoords);
	//float newD = getPointLineDistance(newCoords, startCoords, stopCoords);
	//float halfLipidBilayerWidth = lipidBilayerWidth / 2;

	//if (oldD <= radius - halfLipidBilayerWidth && newD > radius - halfLipidBilayerWidth)
	//	return collision::INSIDE;
	//if (newD < radius + halfLipidBilayerWidth && oldD >= radius + halfLipidBilayerWidth)
	//	return collision::OUTSIDE;

	return collision::NONE;
}

bool Soma::getCollisionPoint(float* point, float newCoords[3], float oldCoords[3], collision::Type collisionType) const
{
	/*float precision = 0.0001f;
	float barrierRadius;
	if (collisionType == collision::INSIDE)
		barrierRadius = radius - lipidBilayerWidth / 2;
	else if (collisionType == collision::OUTSIDE)
		barrierRadius = radius + lipidBilayerWidth / 2;
	else
		return false;

	glm::vec3 currPoint = glm::vec3(oldCoords[0], oldCoords[1], oldCoords[2]);
	glm::vec3 dVec = glm::vec3(newCoords[0], newCoords[1], newCoords[2]) - currPoint;
	dVec *= precision;
	int iterations = 1.0f / precision;

	for (int i = 0; i < iterations; ++i, currPoint += dVec) {
		float d = getPointLineDistance(&currPoint[0], startCoords, stopCoords);
		if (fabs(d - barrierRadius) < precision) {
			point[0] = currPoint[0];
			point[1] = currPoint[1];
			point[2] = currPoint[2];
			return true;
		}
	}*/

	return false;
}

bool Soma::getCollisionNormalVec(float collisionPoint[3], glm::vec3& n, collision::Type collisionType) const
{
	/*float center[3] = { x0, y0, z0 };
	float h = getPointOnLineDistanceFromCenter(collisionPoint, center, radius);
	if (collisionPoint[0] <= x0)
		n = glm::normalize(glm::vec3(x0 - h - collisionPoint[0], collisionPoint[1], collisionPoint[2]));
	else
		n = glm::normalize(glm::vec3(x0 + h - collisionPoint[0], collisionPoint[1], collisionPoint[2]));*/
	return true;
}

bool Soma::getRandPointOnInnerLayer(float* point, glm::vec3& inOutVec) const
{
	//// draw x and angle
	//float x = getRandDouble(startCoords[0], stopCoords[0]);
	//float angle = getRandDouble(0.0, 2 * phy::pi);
	//point[0] = x;
	//point[1] = (radius - lipidBilayerWidth / 2) * sin(angle);
	//point[2] = (radius - lipidBilayerWidth / 2) * cos(angle);

	//inOutVec = glm::normalize(glm::vec3(point[0] - x, point[1] - y0, point[2] - z0));
	return true;
}
