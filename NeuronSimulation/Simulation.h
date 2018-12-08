#pragma once
#include "Utilities.h"
#include "Config.h"
#include "Camera.h"
#include "Particle.h"

class Simulation
{
	Config config;
	Camera camera;
	std::ofstream logfile;

	double metricFactor;
	double timeFactor;
	double inversedTimeFactor;
	int width;
	int height;
	double currentFrame;
	double lastFrame;
	double deltaTime;
	GLFWwindow* window;
	ShaderProgram* ionsRenderProgram;

	GLuint NapIonsVAO;
	GLuint NapIonsPosBuf;
	GLuint NapIonTexture;

	std::vector<Particle> NapIons[2];
	std::vector<double> accels;
	float* NapIonsPos;
	unsigned short bufferNum;

	void loadConfig(const Config& _config);
	void setupOpenGL();
	void setupInput();
	void setupPrograms();
	void setupStructures();
	void setupBuffers();
	static void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods);
	static void mouseButtonCallback(GLFWwindow* window, int button, int action, int mods);
	static void cursorPosCallback(GLFWwindow* window, double xpos, double ypos);
	static void scrollCallback(GLFWwindow* window, double xoffset, double yoffset);
	static void framebufferSizeCallback(GLFWwindow* window, int width, int height);

	bool updateFramebufferSize(int width, int height);

	inline void updateIons();
	inline void update();
	inline void render();

public:
	Simulation(const Config& config);
	~Simulation();
	void start(void);

	double getDeltaTime(void) const;
};

