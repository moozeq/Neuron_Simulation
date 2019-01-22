#pragma once
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>
#include <algorithm>

#include <vector>

#include <glew.h>
#include <glfw3.h>

#include <ext.hpp>

#include "STBImage.h"
#include "Camera.h"
#include "ShaderProgram.h"
#include "Particle.h"
#include "Channel.h"

namespace phy {
	/* physics consts:										*/
	/* angstream a:					10^(-10)				*/
	/* epsilonR, water in 36.6:		~70.0					*/
	/* epsilon0:					8.8542 * 10^(-12) F/m	*/
	/* Na+ radius:					116 * 10^(-12) m		*/
	/* K+ radius:					152 * 10^(-12) m		*/
	/* Na+ channel width:			0.45 nm					*/
	/* K+ channel width:			1.5 nm					*/
	/* Lipid bilayer width:			6.0 nm					*/
	/* Synaptic gap width:			20.0 nm					*/
	/* Na+ channel open treshold:   -40 mV					*/
	/* Na+ channel open time:       1 ms					*/
	/* Na+ channel inactive time:   1 ms					*/

	constexpr double A = 1e-10;
	constexpr double e = 1.60218e-19;
	constexpr double k = 1.2839342030e+8;
	constexpr double u1 = 1.6605389274e-27;
	constexpr double pi = 3.14159265359;

	// scales
	constexpr double tempIonScale = 30;
	constexpr double tempChannelScale = 2;
	constexpr double tempTimeScale = 1e-3;

	constexpr double lipidBilayerWidth = 60 * A;
	constexpr double synapticGapWidth = 200 * A;

	constexpr double NapOpenTreshold = -40e-3;
	constexpr double NapRepolarizationTreshold = -70e-3;
	constexpr double NapOpenTime = 1e-3;

	constexpr double KpOpenTreshold = 40e-3;
	constexpr double KpRepolarizationTreshold = -40e-3;
	constexpr double KpCloseTime = 3e-3;

	constexpr double NtrR = 1.16 * A;
	constexpr double NtrM = 38.1754326758e-27;
	constexpr double NtrC = +1.0 * e;
	constexpr double NtrA = k * NtrC / NtrM;

	constexpr double NapR = 1.16 * A;
	constexpr double NapM = 38.1754326758e-27;
	constexpr double NapC = +1.0 * e;
	constexpr double NapA = k * NapC / NapM;

	constexpr double KpR = 1.52 * A;
	constexpr double KpM = 64.7007990034e-27;
	constexpr double KpC = +1.0 * e;
	constexpr double KpA = k * KpC / KpM;

	constexpr double ClmR = 1.67 * A;
	constexpr double ClmM = 58.0671408979e-27;
	constexpr double ClmC = -1.0 * e;
	constexpr double ClmA = k * ClmC / ClmM;

	constexpr double OanR = 2 * A;
	constexpr double OanM = 1000e-27;
	constexpr double OanC = -1.0 * e;
	constexpr double OanA = k * OanC / OanM;

	constexpr double NapChR = 20 * A;
	constexpr double KpChR = 34 * A;
}

struct UCoord3 {
	unsigned coords[3];
};

struct FCoord8 {
	float coords[8];
};

static inline double getRandDouble(double min, double max) {
	return ((max - min) * ((double)rand() / (double)RAND_MAX) + min);
}

static inline double getGap(double sphereRadius, double cylinderRadius) {
	return sphereRadius - sqrt(sphereRadius * sphereRadius - cylinderRadius * cylinderRadius);
}

static inline float getPointLineDistance(const float point[3], const float startPoint[3], const float stopPoint[3]) {
	glm::vec3 X0X1 = glm::vec3(point[0] - startPoint[0], point[1] - startPoint[1], point[2] - startPoint[2]);
	glm::vec3 X0X2 = glm::vec3(point[0] - stopPoint[0], point[1] - stopPoint[1], point[2] - stopPoint[2]);
	glm::vec3 X2X1 = glm::vec3(stopPoint[0] - startPoint[0], stopPoint[1] - startPoint[1], stopPoint[2] - startPoint[2]);
	return glm::length(glm::cross(X0X1, X0X2)) / glm::length(X2X1);
}

static GLuint generateVAO(std::vector<FCoord8>* vertices, std::vector<UCoord3>* indices) {
	GLuint VAO, VBO, EBO;
	glGenBuffers(1, &VBO);
	glGenBuffers(1, &EBO);
	glGenVertexArrays(1, &VAO);

	glBindVertexArray(VAO);
	glBindBuffer(GL_ARRAY_BUFFER, VBO);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);

	glBufferData(GL_ARRAY_BUFFER, vertices->size() * 8 * sizeof(GLfloat), &(*vertices)[0], GL_STATIC_DRAW);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices->size() * 3 * sizeof(GLuint), &(*indices)[0], GL_STATIC_DRAW);

	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(GLfloat), (GLvoid*)0);
	glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 8 * sizeof(GLfloat), (GLvoid*)(3 * sizeof(GLfloat)));
	glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(GLfloat), (GLvoid*)(5 * sizeof(GLfloat)));
	glEnableVertexAttribArray(0);
	glEnableVertexAttribArray(1);
	glEnableVertexAttribArray(2);

	glBindVertexArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

	return VAO;
}

static void setParticleProperties(Particle* particle, particle::Type type) {
	switch (type) {
	case particle::NEUROTRANSMITTER:
		particle->charge = phy::NtrC;
		particle->mass = phy::NtrM;
		particle->r = phy::NtrR;
		break;

	case particle::NAP:
		particle->charge = phy::NapC;
		particle->mass = phy::NapM;
		particle->r = phy::NapR;
		break;

	case particle::KP:
		particle->charge = phy::KpC;
		particle->mass = phy::KpM;
		particle->r = phy::KpR;
		break;

	case particle::CLM:
		particle->charge = phy::ClmC;
		particle->mass = phy::ClmM;
		particle->r = phy::ClmR;
		break;

	case particle::ORGANIC_ANION:
		particle->charge = phy::OanC;
		particle->mass = phy::OanM;
		particle->r = phy::OanR;
		break;

	default:
		particle->charge = 0.0f;
		particle->mass = phy::u1;
		particle->r = phy::A;
		break;
	}
}

static Particle* newParticle(float x, float y, float z, particle::Type type) {
	Particle* particle = new Particle();
	particle->vx = 0.0;
	particle->vy = 0.0;
	particle->vz = 0.0;
	particle->x = x;
	particle->y = y;
	particle->z = z;

	setParticleProperties(particle, type);

	return particle;
}

static Particle* newParticle(float coords[3], float velocities[3], particle::Type type) {
	Particle* particle = new Particle();
	particle->vx = velocities[0];
	particle->vy = velocities[1];
	particle->vz = velocities[2];
	particle->x = coords[0];
	particle->y = coords[1];
	particle->z = coords[2];

	setParticleProperties(particle, type);

	return particle;
}

static Particle* newParticle(double boundaries[3][2], particle::Type type) {
	Particle* particle = new Particle();
	particle->vx = 0.0;
	particle->vy = 0.0;
	particle->vz = 0.0;
	particle->x = getRandDouble(boundaries[0][0], boundaries[0][1]);
	particle->y = getRandDouble(boundaries[1][0], boundaries[1][1]);
	particle->z = getRandDouble(boundaries[2][0], boundaries[2][1]);

	setParticleProperties(particle, type);

	return particle;
}

static Particle* newParticle(float coords[3], float velocities[3], double charge, double mass, double r) {
	Particle* particle = new Particle();
	particle->vx = velocities[0];
	particle->vy = velocities[1];
	particle->vz = velocities[2];
	particle->x = coords[0];
	particle->y = coords[1];
	particle->z = coords[2];

	particle->charge = charge;
	particle->mass = mass;
	particle->r = r;

	return particle;
}

static void log(std::ofstream& stream, std::string prompt)
{
	time_t t = time(0);
	struct tm now;
	localtime_s(&now, &t);
	std::string currentDate = std::to_string(now.tm_mday) + '.' + std::to_string(now.tm_mon + 1) + '.' + std::to_string(now.tm_year + 1900);
	std::string currentTime = std::to_string(now.tm_hour) + ':' + std::to_string(now.tm_min) + ':' + std::to_string(now.tm_sec);

	// print to log file
	stream << currentDate << ' ' << currentTime << '\t' << prompt << std::endl;

	// print to debug console
	std::cout << currentDate << ' ' << currentTime << '\t' << prompt << std::endl;
}

static void clearAlphaChannel(unsigned char* image, int width, int height) {
	for (int i = 0; i < width * height * 4; i += 4) {
		if (image[i] == 0xff && image[i + 1] == 0xff && image[i + 2] == 0xff)
			image[i + 3] = 0x00;
	}
}

static GLuint loadMipmapTexture(GLuint textureId, const char* fname)
{
	int width, height, comp;

	unsigned char* image = stbi_load(fname, &width, &height, &comp, STBI_rgb_alpha);
	if (image == nullptr)
		throw std::exception("Failed to load texture file");
	clearAlphaChannel(image, width, height);

	GLuint texture;
	glCreateTextures(GL_TEXTURE_2D, 1, &texture);

	glActiveTexture(textureId);
	glBindTexture(GL_TEXTURE_2D, texture);

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, image);
	glGenerateMipmap(GL_TEXTURE_2D);

	glBindTexture(GL_TEXTURE_2D, 0);
	stbi_image_free(image);
	return texture;
}