#pragma once
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>

#include <vector>

#include <glew.h>
#include <glfw3.h>

#include <ext.hpp>

#include "STBImage.h"
#include "ShaderProgram.h"

#define GET_RAND_DOUBLE(min, max) ((max - min) * ( (double)rand() / (double)RAND_MAX ) + min)

namespace phy {
	/* physics consts:										*/
	/* epsilonR, water in 36.6:		~70.0					*/
	/* epsilon0:					8.8542 * 10^(-12) F/m	*/
	/* Na+ radius:					116 * 10^(-12) m		*/
	/* K+ radius:					152 * 10^(-12) m		*/
	/* Na+ channel width:			0.45 nm					*/
	/* K+ channel width:			1.5 nm					*/
	constexpr double e = 1.60218e-19;

	constexpr double NapR = 116e-12;
	constexpr double NapM = 38.1754326758e-27;
	constexpr double NapC = +1.0 * e;

	constexpr double KpR = 152e-12;
	constexpr double KpM = 64.7007990034e-27;
	constexpr double KpC = +1.0 * e;

	constexpr double k = 1.2839342030e+8;
	//constexpr double ENap = 1.5287557234539444552960119488831e+9;
	//constexpr double EKp = 8.9036257854900781641547510319302e+8;

	/* 1.0 in space is 1000nm */
	/* 1.0 sec real time is 1ms */
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
	delete(image);
	return texture;
}

static GLFWcursor* createCursor(std::string cursorImage)
{
	int cursorWidth, cursorHeight, comp;
	unsigned char* pixels = stbi_load(cursorImage.c_str(), &cursorWidth, &cursorHeight, &comp, STBI_rgb_alpha);

	// TODO make pics with alpha channel active, here's some workaround
	clearAlphaChannel(pixels, cursorWidth, cursorHeight);
	GLFWimage image;
	image.width = cursorWidth;
	image.height = cursorHeight;
	image.pixels = pixels;
	GLFWcursor* cursor = glfwCreateCursor(&image, 0, 0);
	delete(pixels);
	return cursor;
}

