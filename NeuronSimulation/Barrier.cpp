#include "Barrier.h"

GLuint Barrier::generateCircles()
{
	addCircle();
	addCircle();
	return generateVAO(&circlesVertices, &circlesIndices);
}

GLuint Barrier::generateBilayer()
{
	GLuint offset = slices + 2; // 2 midpoints
	for (GLuint i = 1; i < slices; ++i) {
		bilayerIndices.push_back(UCoord3({ i, i + 1, i + offset }));
		bilayerIndices.push_back(UCoord3({ i + 1, i + offset, i + offset + 1 }));
	}
	bilayerIndices.push_back(UCoord3({ slices, slices + 1, slices + offset + 1 }));
	bilayerIndices.push_back(UCoord3({ slices, slices + offset, slices + offset + 1 }));

	return generateVAO(&bilayerVertices, &bilayerIndices);
}

void Barrier::addCircle()
{
	GLuint offset = circlesVertices.size();
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
	circlesVertices.push_back(FCoord8({ x + x0, y0, z0, 0.5f, 0.5f, normX, 0.0f, 0.0f })); //midpoint
	bilayerVertices.push_back(FCoord8({ x + x0, y0, z0, 0.5f, 0.5f, normX, 0.0f, 0.0f }));

	for (GLuint i = 1; i < slices + 2; ++i, alfa += theta) {
		sa = sin(alfa);
		ca = cos(alfa);
		z = radius * ca;
		y = radius * sa;
		circlesVertices.push_back(FCoord8({ x + x0, y + y0, z + z0, multi * (y + radius), multi * (z + radius), normX, sa, ca }));
		bilayerVertices.push_back(FCoord8({ x + x0, y + y0, z + z0, radius * alfa, corner, normX, sa, ca }));

		circlesIndices.push_back(UCoord3({ offset, i + offset, i + offset + 1 }));
	}
	circlesIndices.pop_back();
}

Barrier::Barrier(float coords[3], float _radius, float _length) :
	x0(coords[0]), y0(coords[1]), z0(coords[2]), radius(_radius), length(_length), slices(30)
{
	circleVAO = generateCircles();
	layerVAO = generateBilayer();

	circlesTexture = loadMipmapTexture(GL_TEXTURE0, "circlesBilayer.png");
	bilayerTexture = loadMipmapTexture(GL_TEXTURE0, "bilayer.png");
}

bool Barrier::checkCollision(float newCoords[3], float oldCoords[3])
{
	// TODO
	return false;
}

bool Barrier::getCollisionPoint(float* point, float newCoords[3], float oldCoords[3])
{
	// TODO
	point[0] = point[1] = point[2] = 0.0f;

	return true;
}

void Barrier::render()
{
	// only open barriers for now

	/*glBindTexture(GL_TEXTURE_2D, circlesTexture);
	glBindVertexArray(circleVAO);
	glDrawElements(GL_TRIANGLES, circlesIndices.size() * 3, GL_UNSIGNED_INT, 0);*/

	glBindTexture(GL_TEXTURE_2D, bilayerTexture);
	glBindVertexArray(layerVAO);
	glDrawElements(GL_TRIANGLES, bilayerIndices.size() * 3, GL_UNSIGNED_INT, 0);
}
