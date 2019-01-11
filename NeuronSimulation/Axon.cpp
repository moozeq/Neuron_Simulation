#include "Axon.h"

GLuint Axon::generateLayer(std::vector<FCoord8>& vertices, std::vector<UCoord3>& indices)
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

void Axon::addDisc(std::vector<UCoord3>& indices, unsigned midPointIndex)
{
	for (unsigned i = 1; i < slices + 1; ++i)
		indices.push_back(UCoord3({ midPointIndex, i + midPointIndex, i + 1 + midPointIndex }));
}

Axon::Axon(float coords[3], float _radius, float _length, float _lipidBilayerWidth) :
	length(_length), slices(30)
{
	precision = 0.01f;
	radius = _radius;
	lipidBilayerWidth = _lipidBilayerWidth;
	innerRadius = _radius - _lipidBilayerWidth / 2.0f;
	outerRadius = _radius + _lipidBilayerWidth / 2.0f;

	x0 = coords[0];
	y0 = coords[1];
	z0 = coords[2];

	startCoords[0] = coords[0] - _length / 2.0f - _lipidBilayerWidth / 2.0f;
	startCoords[1] = coords[1];
	startCoords[2] = coords[2];

	stopCoords[0] = coords[0] + _length / 2.0f + _lipidBilayerWidth / 2.0f;
	stopCoords[1] = coords[1];
	stopCoords[2] = coords[2];

	discPoint[0] = x0 + _length / 2.0f - _lipidBilayerWidth / 2.0f;
	discPoint[1] = y0;
	discPoint[2] = z0;

	// add inner layer
	length -= _lipidBilayerWidth;
	radius = innerRadius;
	addCircle(innerLayerVertices);
	addCircle(innerLayerVertices);
	addDisc(innerLayerIndices, slices + 2);
	innerLayerVAO = generateLayer(innerLayerVertices, innerLayerIndices);

	// add outer layer
	length += 2.0f * _lipidBilayerWidth;
	radius = outerRadius;
	addCircle(outerLayerVertices);
	addCircle(outerLayerVertices);
	addDisc(outerLayerIndices, slices + 2);
	outerLayerVAO = generateLayer(outerLayerVertices, outerLayerIndices);

	length = _length;
	radius = _radius;
}

collision::Type Axon::checkCollision(const float newCoords[3], const float oldCoords[3]) const
{
	// not in barrier length
	if (oldCoords[0] < startCoords[0] && newCoords[0] < startCoords[0] ||
		oldCoords[0] > stopCoords[0] && newCoords[0] > stopCoords[0])
		return collision::NONE;

	float oldD = getPointLineDistance(oldCoords, startCoords, stopCoords);
	float newD = getPointLineDistance(newCoords, startCoords, stopCoords);

	if (oldD < innerRadius && newD > innerRadius)
		return collision::INSIDE;
	if (newD < outerRadius && oldD > outerRadius)
		return collision::OUTSIDE;
	if (oldD < innerRadius &&
		newD < innerRadius &&
		oldCoords[0] < stopCoords[0] - lipidBilayerWidth && 
		newCoords[0] > stopCoords[0] - lipidBilayerWidth)
		return collision::DISC;

	// TODO better inside bilayer collision
	if (oldD < radius && newD > radius)
		return collision::INSIDE;

	return collision::NONE;
}

int Axon::getCollisionPoint(float* point, float newCoords[3], float oldCoords[3], collision::Type collisionType) const
{
	glm::vec3 currPoint = glm::vec3(oldCoords[0], oldCoords[1], oldCoords[2]);
	glm::vec3 dVec = glm::vec3(newCoords[0], newCoords[1], newCoords[2]) - currPoint;
	dVec *= precision;
	int iterations = 1 / precision;

	if (collisionType == collision::INSIDE) {
		for (int i = 0; i < iterations; ++i, currPoint += dVec) {
			float d = getPointLineDistance(&currPoint[0], startCoords, stopCoords);
			if (fabs(d - innerRadius) < precision) {
				point[0] = currPoint[0];
				point[1] = currPoint[1];
				point[2] = currPoint[2];
				return i;
			}
		}
	}
	else if (collisionType == collision::OUTSIDE) {
		for (int i = 0; i < iterations; ++i, currPoint += dVec) {
			float d = getPointLineDistance(&currPoint[0], startCoords, stopCoords);
			if (fabs(d - outerRadius) < precision) {
				point[0] = currPoint[0];
				point[1] = currPoint[1];
				point[2] = currPoint[2];
				return i;
			}
		}
	}
	else if (collisionType == collision::DISC) {
		point[0] = discPoint[0];
		point[1] = newCoords[1];
		point[2] = newCoords[2];
		return -3;
	}
	else
		return -2;

	point[0] = (newCoords[0] - oldCoords[0]) / 2.0f;
	point[1] = (newCoords[1] - oldCoords[1]) / 2.0f;
	point[2] = (newCoords[2] - oldCoords[2]) / 2.0f;
	return -1;
}

bool Axon::getCollisionNormalVec(float collisionPoint[3], glm::vec3& n, collision::Type collisionType) const
{
	if (collisionType == collision::DISC) {
		n = glm::vec3(-1.0f, 0.0f, 0.0f);
		return true;
	}

	n = glm::normalize(glm::vec3(0.0f, collisionPoint[1], collisionPoint[2]));

	if (collisionType == collision::OUTSIDE)
		n *= -1.0f;
	return true;
}

bool Axon::getRandPointOnInnerLayer(float* point, glm::vec3& inOutVec) const
{
	// draw x and angle
	float synapse = getRandDouble(0.0, 1.0);
	float angle = getRandDouble(0.0, 2 * phy::pi);
	float x, r;
	if (synapse < synapseProbability) {
		x = stopCoords[0] - lipidBilayerWidth;
		r = getRandDouble(0.0, 1.0) * innerRadius;
		inOutVec = glm::normalize(glm::vec3(1.0f, 0.0f, 0.0f));
	}
	else {
		x = getRandDouble(startCoords[0] + lipidBilayerWidth, stopCoords[0] - lipidBilayerWidth);;
		r = innerRadius;
		inOutVec = glm::normalize(glm::vec3(0, point[1] - y0, point[2] - z0));
	}

	point[0] = x;
	point[1] = r * sin(angle);
	point[2] = r * cos(angle);

	return true;
}

void Axon::setSynapseProbability(float probability)
{
	synapseProbability = probability;
}
