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

	double metricFactorSq;
	double timeFactor;
	double inversedTimeFactor;
	int width;
	int height;
	bool ice;
	double currentFrame;
	double lastFrame;
	double deltaTime;
	GLFWwindow* window;

	ShaderProgram* ionsRenderProgram;

	GLuint NapIonsVAO;
	GLuint KpIonsVAO;
	GLuint ClmIonsVAO;
	GLuint otherParticlesVAO;

	GLuint NapIonTexture;
	GLuint KpIonTexture;
	GLuint ClmIonTexture;
	GLuint otherParticlesTexture;

	GLuint particlesPosBuf;

	std::vector<Particle> particles[2];
	std::vector<double> accels;
	std::vector<double> partAccOrigin;
	float* particlesPos;
	unsigned short bufferNum;
	size_t particlesBufferSize;

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

	inline void updateParticles();
	inline void update();
	inline void render();

public:
	Simulation(const Config& config);
	~Simulation();
	void start(void);
	void freeze(void);

	double getDeltaTime(void) const;
};