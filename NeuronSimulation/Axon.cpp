#include "Axon.h"

GLuint Axon::generateBilayer(std::vector<FCoord8>& vertices, std::vector<UCoord3>& indices)
{
	GLuint offset = slices + 2; // 2 midpoints
	for (GLuint i = 1; i < slices; ++i) {
		indices.push_back(UCoord3({ i, i + 1, i + offset }));
		indices.push_back(UCoord3({ i + 1, i + offset, i + offset + 1 }));
	}
	indices.push_back(UCoord3({ slices, slices + 1, slices + offset + 1 }));
	indices.push_back(UCoord3({ slices, slices + offset, slices + offset + 1 }));

	return generateVAO(&vertices, &indices);
}

void Axon::addCircle(std::vector<FCoord8>& vertices)
{
	GLuint offset = vertices.size();
	GLfloat x = -length / 2.0f;
	GLfloat normX = -1.0f;
	GLfloat corner = 0.0f;
	if (offset > 0) {
		x = -x;
		normX = -normX;
		corner = 1.0f;
	}
	GLfloat theta = (2.0f * phy::pi) / (float)slices;
	GLfloat alfa = -phy::pi / 2.0f;
	GLfloat y, z, sa, ca;
	GLfloat multi = 1.0f / (2.0f * radius);
	vertices.push_back(FCoord8({ x + x0, y0, z0, 0.5f, 0.5f, normX, 0.0f, 0.0f }));

	for (GLuint i = 1; i < slices + 2; ++i, alfa += theta) {
		sa = sin(alfa);
		ca = cos(alfa);
		z = radius * ca;
		y = radius * sa;
		vertices.push_back(FCoord8({ x + x0, y + y0, z + z0, radius * alfa, corner, normX, sa, ca }));
	}
}

Axon::Axon(float coords[3], float _radius, float _length, float _lipidBilayerWidth) :
	radius(_radius), length(_length), slices(30), lipidBilayerWidth(_lipidBilayerWidth)
{
	x0 = coords[0];
	y0 = coords[1];
	z0 = coords[2];

	startCoords[0] = coords[0] - _length / 2;
	startCoords[1] = coords[1];
	startCoords[2] = coords[2];

	stopCoords[0] = coords[0] + _length / 2;
	stopCoords[1] = coords[1];
	stopCoords[2] = coords[2];

	radius -= lipidBilayerWidth / 2;

	addCircle(innerLayerVertices);
	addCircle(innerLayerVertices);
	innerLayerVAO = generateBilayer(innerLayerVertices, innerLayerIndices);

	radius += lipidBilayerWidth;

	addCircle(outerLayerVertices);
	addCircle(outerLayerVertices);
	outerLayerVAO = generateBilayer(outerLayerVertices, outerLayerIndices);

	radius -= lipidBilayerWidth / 2;
}

collision::Type Axon::checkCollision(const float newCoords[3], const float oldCoords[3]) const
{
	// not in barrier length
	if (!(newCoords[0] <= stopCoords[0] && newCoords[0] >= startCoords[0] &&
		oldCoords[0] <= stopCoords[0] && oldCoords[0] >= startCoords[0]))
		return collision::NONE;

	float oldD = getPointLineDistance(oldCoords, startCoords, stopCoords);
	float newD = getPointLineDistance(newCoords, startCoords, stopCoords);
	float halfLipidBilayerWidth = lipidBilayerWidth / 2;

	if (oldD <= radius - halfLipidBilayerWidth && newD > radius - halfLipidBilayerWidth)
		return collision::INSIDE;
	if (newD < radius + halfLipidBilayerWidth && oldD >= radius + halfLipidBilayerWidth)
		return collision::OUTSIDE;

	return collision::NONE;
}

bool Axon::getCollisionPoint(float* point, float newCoords[3], float oldCoords[3], collision::Type collisionType) const
{
	float precision = 0.0001f;
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
	}

	return false;
}

bool Axon::getCollisionNormalVec(float collisionPoint[3], glm::vec3& n, collision::Type collisionType) const
{
	float center[3] = { x0, y0, z0 };
	float h = getPointOnLineDistanceFromCenter(collisionPoint, center, radius);
	if (collisionPoint[0] <= x0)
		n = glm::normalize(glm::vec3(x0 - h - collisionPoint[0], collisionPoint[1], collisionPoint[2]));
	else
		n = glm::normalize(glm::vec3(x0 + h - collisionPoint[0], collisionPoint[1], collisionPoint[2]));
	return true;
}

bool Axon::getRandPointOnInnerLayer(float* point, glm::vec3& inOutVec) const
{
	// draw x and angle
	float x = getRandDouble(startCoords[0], stopCoords[0]);
	float angle = getRandDouble(0.0, 2 * phy::pi);
	point[0] = x;
	point[1] = (radius - lipidBilayerWidth / 2) * sin(angle);
	point[2] = (radius - lipidBilayerWidth / 2) * cos(angle);

	inOutVec = glm::normalize(glm::vec3(point[0] - x, point[1] - y0, point[2] - z0));
	return true;
}