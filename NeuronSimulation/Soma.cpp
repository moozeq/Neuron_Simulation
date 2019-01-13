#include "Soma.h"

GLuint Soma::generateLayer(std::vector<FCoord8>& vertices, std::vector<UCoord3>& indices)
{
	// create 12 vertices of a icosahedron
	float scale = radius / 1.5f;
	float t = ((1 + sqrt(5.0f)) / 2.0f) * scale;
	float half = 1.0f * scale;

	// TODO add texture and proper normals
	vertices.push_back(FCoord8({ -half + x0, t + y0, 0 + z0, 0, 0, 0, 0, 0 }));
	vertices.push_back(FCoord8({ half + x0, t + y0, 0 + z0, 0, 0, 0, 0, 0 }));
	vertices.push_back(FCoord8({ -half + x0, -t + y0, 0 + z0, 0, 0, 0, 0, 0 }));
	vertices.push_back(FCoord8({ half + x0, -t + y0, 0 + z0, 0, 0, 0, 0, 0 }));

	vertices.push_back(FCoord8({ 0 + x0, -half + y0, t + z0, 0, 0, 0, 0, 0 }));
	vertices.push_back(FCoord8({ 0 + x0, half + y0, t + z0, 0, 0, 0, 0, 0 }));
	vertices.push_back(FCoord8({ 0 + x0, -half + y0, -t + z0, 0, 0, 0, 0, 0 }));
	vertices.push_back(FCoord8({ 0 + x0, half + y0, -t + z0, 0, 0, 0, 0, 0 }));

	vertices.push_back(FCoord8({ t + x0, 0 + y0, -half + z0, 0, 0, 0, 0, 0 }));
	vertices.push_back(FCoord8({ t + x0, 0 + y0, half + z0, 0, 0, 0, 0, 0 }));
	vertices.push_back(FCoord8({ -t + x0, 0 + y0, -half + z0, 0, 0, 0, 0, 0 }));
	vertices.push_back(FCoord8({ -t + x0, 0 + y0, half + z0, 0, 0, 0, 0, 0 }));

	// 5 faces around point 0
	indices.push_back(UCoord3({ 0, 11, 5 }));
	indices.push_back(UCoord3({ 0, 5, 1 }));
	indices.push_back(UCoord3({ 0, 1, 7 }));
	indices.push_back(UCoord3({ 0, 7, 10 }));
	indices.push_back(UCoord3({ 0, 10, 11 }));

	// 5 adjacent faces
	indices.push_back(UCoord3({ 1, 5, 9 }));
	indices.push_back(UCoord3({ 5, 11, 4 }));
	indices.push_back(UCoord3({ 11, 10, 2 }));
	indices.push_back(UCoord3({ 10, 7, 6 }));
	indices.push_back(UCoord3({ 7, 1, 8 }));

	// 5 faces around point 3
	indices.push_back(UCoord3({ 3, 9, 4 }));
	indices.push_back(UCoord3({ 3, 4, 2 }));
	indices.push_back(UCoord3({ 3, 2, 6 }));
	indices.push_back(UCoord3({ 3, 6, 8 }));
	indices.push_back(UCoord3({ 3, 8, 9 }));

	// 5 adjacent faces
	indices.push_back(UCoord3({ 4, 9, 5 }));
	indices.push_back(UCoord3({ 2, 4, 11 }));
	indices.push_back(UCoord3({ 6, 2, 10 }));
	indices.push_back(UCoord3({ 8, 6, 7 }));
	indices.push_back(UCoord3({ 9, 8, 1 }));

	return generateVAO(&vertices, &indices);
}

Soma::Soma(float coords[3], float _radius, float _lipidBilayerWidth)
{
	precision = 0.0001f;
	radius = _radius;
	lipidBilayerWidth = _lipidBilayerWidth;
	innerRadius = _radius - _lipidBilayerWidth / 2.0f;
	outerRadius = _radius + _lipidBilayerWidth / 2.0f;
	x0 = coords[0];
	y0 = coords[1];
	z0 = coords[2];

	radius -= _lipidBilayerWidth / 2.0f;
	innerLayerVAO = generateLayer(innerLayerVertices, innerLayerIndices);

	radius += _lipidBilayerWidth;
	outerLayerVAO = generateLayer(outerLayerVertices, outerLayerIndices);

	radius = _radius;
}

collision::Type Soma::checkCollision(const float newCoords[3], const float oldCoords[3]) const
{
	// TODO not only X coord but point
	for (const float pointX : connectPointsXRight) {
		if (newCoords[0] > pointX && oldCoords[0] > pointX)
			return collision::NONE;
	}
	for (const float pointX : connectPointsXLeft) {
		if (newCoords[0] < pointX&& oldCoords[0] < pointX)
			return collision::NONE;
	}

	float newD = glm::distance(glm::vec3(newCoords[0], newCoords[1], newCoords[2]), glm::vec3(x0, y0, z0));
	float oldD = glm::distance(glm::vec3(oldCoords[0], oldCoords[1], oldCoords[2]), glm::vec3(x0, y0, z0));

	if (oldD < innerRadius && newD > innerRadius)
		return collision::INSIDE;
	if (newD < outerRadius && oldD > outerRadius)
		return collision::OUTSIDE;

	// TODO better inside bilayer collision
	if (oldD < radius && newD > radius)
		return collision::INSIDE;

	return collision::NONE;
}

int Soma::getCollisionPoint(float* point, float newCoords[3], float oldCoords[3], collision::Type collisionType) const
{
	// TODO better finding collision point
	float barrierRadius;
	if (collisionType == collision::INSIDE)
		barrierRadius = innerRadius;
	else if (collisionType == collision::OUTSIDE)
		barrierRadius = outerRadius;
	else
		return -2;

	glm::vec3 currPoint = glm::vec3(oldCoords[0], oldCoords[1], oldCoords[2]);
	glm::vec3 dVec = glm::vec3(newCoords[0], newCoords[1], newCoords[2]) - currPoint;
	glm::vec3 midPoint = glm::vec3(x0, y0, z0);
	dVec *= precision;
	int iterations = 1 / precision;

	for (int i = 0; i < iterations; ++i, currPoint += dVec) {
		float d = glm::distance(currPoint, midPoint);
		if (fabs(d - barrierRadius) < precision) {
			point[0] = currPoint[0];
			point[1] = currPoint[1];
			point[2] = currPoint[2];
			return i;
		}
	}

	point[0] = (newCoords[0] - oldCoords[0]) / 2.0f;
	point[1] = (newCoords[1] - oldCoords[1]) / 2.0f;
	point[2] = (newCoords[2] - oldCoords[2]) / 2.0f;
	return -1;
}

bool Soma::getCollisionNormalVec(float collisionPoint[3], glm::vec3& n, collision::Type collisionType) const
{
	glm::vec3 midPoint = glm::vec3(x0, y0, z0);
	n = glm::normalize(glm::vec3(collisionPoint[0], collisionPoint[1], collisionPoint[2]) - midPoint);
	
	if (collisionType == collision::OUTSIDE)
		n *= -1.0f;
	return true;
}

bool Soma::getRandPointOnInnerLayer(float* point, glm::vec3& inOutVec) const
{
	float u = getRandDouble(0.0, 1.0);
	float v = getRandDouble(0.0, 1.0);
	float theta = 2 * phy::pi * u;
	float phi = acos(2 * v - 1);
	point[0] = x0 + (innerRadius * sin(phi) * cos(theta));
	point[1] = y0 + (innerRadius * sin(phi) * sin(theta));
	point[2] = z0 + (innerRadius * cos(phi));
	inOutVec = glm::normalize(glm::vec3(point[0] - x0, point[1] - y0, point[2] - z0));

	for (const float pointX : connectPointsXRight) {
		if (point[0] > pointX)
			return false;
	}
	for (const float pointX : connectPointsXLeft) {
		if (point[0] < pointX)
			return false;
	}
	return true;
}

void Soma::addConnection(float connectionPoint[3], float connectionRadius, barrier::Connection connectionType)
{
	float h = getGap(radius, connectionRadius);

	if (connectionType == barrier::SOMA_AXON)
		connectPointsXRight.push_back(connectionPoint[0] - h);
	else if (connectionType == barrier::DENDRITE_SOMA)
		connectPointsXLeft.push_back(connectionPoint[0] + h);

	connectedBarrierRadius.push_back(connectionRadius);
}
