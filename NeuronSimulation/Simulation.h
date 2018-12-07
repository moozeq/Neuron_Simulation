#pragma once
#include "Utilities.h"
#include "Config.h"
#include "Camera.h"

class Simulation
{
	Config config;
	Camera camera;
	std::ofstream logfile;

	int width;
	int height;
	double currentFrame;
	double lastFrame;
	double deltaTime;
	GLFWwindow* window;
	ShaderProgram* ionsRenderProgram;

	GLuint NapIonsVAO;
	GLuint NapIonsPos;
	GLuint NapIonTexture;

	float* NapIons;

	void loadConfig(const Config& _config);
	void setupOpenGL();
	void setupInput();
	void setupPrograms();
	void setupBuffers();
	static void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods);
	static void mouseButtonCallback(GLFWwindow* window, int button, int action, int mods);
	static void cursorPosCallback(GLFWwindow* window, double xpos, double ypos);
	static void scrollCallback(GLFWwindow* window, double xoffset, double yoffset);
	static void framebufferSizeCallback(GLFWwindow* window, int width, int height);

	bool updateFramebufferSize(int width, int height);

public:
	Simulation(const Config& config);
	~Simulation();
	void render(double time);
	void start(void);
};

