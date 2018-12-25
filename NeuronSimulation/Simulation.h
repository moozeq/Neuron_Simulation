#pragma once
#include "Utilities.h"
#include "Config.h"
#include "Camera.h"
#include "Particle.h"
#include "Channel.h"
#include "Neuron.h"

class Simulation
{

/* -----------------SIMULATION CONFIG----------------- */

	// current camera
	Camera camera;

	// const config variables
	Config config;
	std::ofstream logfile;
	double metricFactorSq;
	double inversedTimeFactor;
	unsigned short bufferNum;
	int width;
	int height;

	// runtime config variables
	bool ice;
	bool rewind;
	double timeFactor;
	double currentFrame;
	double lastFrame;
	double deltaTime;

	// opengl window struct
	GLFWwindow* window;

/* -----------------PARTICLES----------------- */

	// particles rendering program
	ShaderProgram* ionsRenderProgram;

	// particles VAOs
	GLuint NapIonsVAO;
	GLuint KpIonsVAO;
	GLuint ClmIonsVAO;
	GLuint otherParticlesVAO;

	// particles textures
	GLuint NapIonTexture;
	GLuint KpIonTexture;
	GLuint ClmIonTexture;
	GLuint otherParticlesTexture;

	// particles structs
	std::vector<Particle> particles[2];
	std::vector<double> partAccOrigin;
	long particlesBufferSize;
	float* particlesPos;

/* -----------------CHANNELS----------------- */

	// channels rendering program
	ShaderProgram* channelsRenderProgram;

	// channels VAOs
	GLuint NapIonsChannelsVAO;
	GLuint KpIonsChannelsVAO;

	// channels textures
	GLuint NapIonChannelTexture[channel::STATES_COUNT];
	GLuint KpIonChannelTexture[channel::STATES_COUNT];

	// channels structs
	long channelsBufferSize;
	float* channelsStates;

/* -----------------NEURON----------------- */

	// neuron struct
	Neuron* neuron;

/* -----------------PRIVATE METHODS----------------- */

	void loadConfig(const Config& _config);
	void setupOpenGL();
	void setupInput();
	void setupPrograms();
	void setupStructures();
	void setupNeuronStructures();
	void setupParticlesStructures();
	void setupTextures();
	void setupBuffers();
	void setupParticlesBuffers();
	void setupChannelsBuffers();

	static void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods);
	static void mouseButtonCallback(GLFWwindow* window, int button, int action, int mods);
	static void cursorPosCallback(GLFWwindow* window, double xpos, double ypos);
	static void scrollCallback(GLFWwindow* window, double xoffset, double yoffset);
	static void framebufferSizeCallback(GLFWwindow* window, int width, int height);

	bool updateFramebufferSize(int width, int height);

	inline void updateChannelsStates();
	inline void updateParticlesPositions();
	inline void update();
	inline void render();

public:

/* -----------------PUBLIC METHODS----------------- */

	Simulation(const Config& config);
	~Simulation();
	void start(void);
	void freeze(void);
	void reverse(void);

	double getDeltaTime(void) const;
};