#include "Simulation.h"

void Simulation::loadConfig(const Config& _config)
{
	// copy config struct
	config = _config;

	// set variables based on config struct
	metricFactorSq = config.metricFactor * config.metricFactor;
	metricFactor = config.metricFactor;
	if (metricFactor <= 0.0)
		throw std::exception("Wrong metric factor");

	timeFactor = config.timeFactor;
	if (timeFactor <= 0.0)
		throw std::exception("Wrong time factor");
	bufferNum = 0;

	inversedTimeFactor = 1.0 / config.timeFactor;
	particlesBufferSize = 0;
	for (int type = 0; type < particle::TYPES_COUNT; ++type) {
		particlesOffset[type] = particlesBufferSize;
		particlesBufferSize += config.particlesCount[type];
		activeParticlesCount[type] = 0;
	}
	if (particlesBufferSize < 0)
		throw std::exception("Wrong particles number");

	ice = false;
	rewind = false;
	cursor = false;
	renderParticles = true;
	renderChannels = true;
}

void Simulation::setupOpenGL()
{
	// log info about client
	logfile.open(config.logPath, std::ofstream::app | std::ofstream::binary);
	if (!logfile.good())
		throw std::exception("Couldn't open log file");

	// setup client's openGL components
	if (glfwInit() != GL_TRUE)
		throw std::exception("GLFW initialization's failed");

	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 5);

	width = config.width;
	height = config.height;

	window = glfwCreateWindow(width, height, "ProjectN", nullptr, nullptr);
	if (window == nullptr)
		throw std::exception("GLFW window couldn't be created");
	glfwMakeContextCurrent(window);

	// new glew functions usage
	glewExperimental = GL_TRUE;
	if (glewInit() != GLEW_OK)
		throw std::exception("GLEW initialization's failed");

	// enable alpha channel in textures
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_DEPTH_TEST);
	glDisable(GL_CULL_FACE);

	glViewport(0, 0, width, height);

	// set user pointer to this simulation (need for input)
	glfwSetWindowUserPointer(window, this);
}

void Simulation::setupInput()
{
	// input functions
	glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
	glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_TRUE);
	glfwSetInputMode(window, GLFW_STICKY_MOUSE_BUTTONS, GL_FALSE);
	glfwSetKeyCallback(window, keyCallback);
	glfwSetMouseButtonCallback(window, mouseButtonCallback);
	glfwSetCursorPosCallback(window, cursorPosCallback);
	glfwSetScrollCallback(window, scrollCallback);
	glfwSetFramebufferSizeCallback(window, framebufferSizeCallback);
}

void Simulation::setupPrograms()
{
	std::vector<const GLchar*> ionsPaths({
		   "Ions.vert",
		   "Ions.geom",
		   "Ions.frag"
		});
	ionsRenderProgram = new ShaderProgram(shader::VFG, ionsPaths);

	std::vector<const GLchar*> channelsPaths({
		   "Channels.vert",
		   "Channels.geom",
		   "Channels.frag"
		});
	channelsRenderProgram = new ShaderProgram(shader::VFG, channelsPaths);
}

void Simulation::setupStructures()
{
	setupNeuronStructures();
	setupParticlesStructures();
}

void Simulation::setupNeuronStructures()
{
	neuron = new Neuron(config);

	channelsBufferSize = 0;
	for (int type = 0; type < channel::TYPES_COUNT; ++type) {
		channelsOffset[type] = channelsBufferSize;
		activeChannelsCount[type] = neuron->getChannelsCount((channel::Type)type);
		channelsBufferSize += activeChannelsCount[type];
	}

	synapsePosition = neuron->getSynapsePosition();
}

void Simulation::setupParticlesStructures()
{
	particles[0].reserve(particlesBufferSize);
	particles[1].reserve(particlesBufferSize);
	double boundaries[3][2] = { { 0, 0 }, { 0, 0 }, { 0, 0 } };

	for (int type = 0; type < particle::TYPES_COUNT; ++type) {
		particle::Type currentType = (particle::Type)type;
		for (int i = 0; i < config.particlesCount[currentType]; ++i) {
			Particle* particle = newParticle(boundaries, currentType);
			particles[0].push_back(*particle);
			particles[1].push_back(*particle);

			// for all particle types add one partAccOrigin
			if (i == 0)
				partAccOrigin[currentType] = (phy::k * particle->charge / particle->mass) / metricFactorSq;
			delete particle;
		}
	}
}

void Simulation::setupUniforms()
{
	particleRadius[particle::NAP] = phy::NapR / metricFactor;
	particleRadius[particle::KP] = phy::KpR / metricFactor;
	particleRadius[particle::CLM] = phy::ClmR / metricFactor;
	particleRadius[particle::ORGANIC_ANION] = phy::OanR / metricFactor;
	particleRadius[particle::NEUROTRANSMITTER] = phy::NtrR / metricFactor;

	channelRadius[channel::NAP] = phy::NapChR / metricFactor;
	channelRadius[channel::KP] = phy::KpChR / metricFactor;

	// set const values as uniforms in shader program
	uniforms.channelWidth = phy::lipidBilayerWidth / metricFactor;
	uniforms.opacity = 0.5f;

	channelsRenderProgram->use();
	channelsRenderProgram->setUniforms(uniforms);
	glUseProgram(0);
}

void Simulation::setupTextures()
{
	// load particles textures
	for (int type = 0; type < particle::TYPES_COUNT; ++type)
		particleTexture[(particle::Type)type] = loadMipmapTexture(GL_TEXTURE0, config.particlesTextures[(particle::Type)type].c_str());

	// load channels textures
	for (int type = 0; type < channel::TYPES_COUNT; ++type)
		channelTexture[(channel::Type)type] = loadMipmapTexture(GL_TEXTURE0, config.channelsTextures[(channel::Type)type].c_str());
}

void Simulation::setupBuffers()
{
	setupParticlesBuffers();
	setupChannelsBuffers();
}

void Simulation::setupParticlesBuffers()
{
	if (particlesBufferSize == 0)
		return;

	GLuint particlesPosBuf;

	// create vaos
	glCreateVertexArrays(particle::TYPES_COUNT, particleVAO);

	// creating buffers
	glCreateBuffers(1, &particlesPosBuf);

	GLbitfield flags = GL_MAP_WRITE_BIT | GL_MAP_PERSISTENT_BIT | GL_MAP_COHERENT_BIT;
	GLsizei bufferSize = particlesBufferSize * 3 * sizeof(GLfloat);
	glNamedBufferStorage(particlesPosBuf, bufferSize, nullptr, flags);

	for (int type = 0; type < particle::TYPES_COUNT; ++type) {
		glVertexArrayVertexBuffer(particleVAO[(particle::Type)type], 0, particlesPosBuf, 0, 3 * sizeof(GLfloat));
		glVertexArrayAttribFormat(particleVAO[(particle::Type)type], 0, 3, GL_FLOAT, GL_FALSE, 0);
		glVertexArrayAttribBinding(particleVAO[(particle::Type)type], 0, 0);
		glEnableVertexArrayAttrib(particleVAO[(particle::Type)type], 0);
	}

	// initialize buffers in GPU and get pointers to them
	particlesPos = (float*)glMapNamedBuffer(particlesPosBuf, GL_WRITE_ONLY);
	if (!particlesPos)
		throw std::exception("Buffer mapping failed");

	glUnmapNamedBuffer(particlesPosBuf);
}

void Simulation::setupChannelsBuffers()
{
	if (channelsBufferSize == 0)
		return;

	GLuint channelsBuf;

	// create vaos
	glCreateVertexArrays(channel::TYPES_COUNT, channelVAO);

	// creating buffers
	glCreateBuffers(1, &channelsBuf);

	GLbitfield flags = GL_MAP_WRITE_BIT | GL_MAP_PERSISTENT_BIT | GL_MAP_COHERENT_BIT;
	std::vector<float> channels = neuron->getChannels();
	glNamedBufferStorage(channelsBuf, channelsBufferSize * 4 * sizeof(GLfloat), &channels[0], flags);

	// positions
	for (int type = 0; type < channel::TYPES_COUNT; ++type) {
		channel::Type currentType = (channel::Type)type;
		glVertexArrayVertexBuffer(channelVAO[currentType], 0, channelsBuf, 0, 4 * sizeof(GLfloat));
		glVertexArrayAttribFormat(channelVAO[currentType], 0, 4, GL_FLOAT, GL_FALSE, 0);
		glVertexArrayAttribBinding(channelVAO[currentType], 0, 0);
		glEnableVertexArrayAttrib(channelVAO[currentType], 0);
	}

	// initialize buffers in GPU and get pointers to them
	channelsAttribs = (float*)glMapNamedBuffer(channelsBuf, GL_WRITE_ONLY);
	if (!channelsAttribs)
		throw std::exception("Buffer mapping failed");

	glUnmapNamedBuffer(channelsBuf);

	for (long i = 0; i < neuron->channels.size(); ++i) {
		const Channel& channel = neuron->channels[i];
		float color;
		switch (channel.state) {
		case channel::OPEN:
			color = 0.0f;
			break;
		case channel::INACTIVE:
			color = 0.5f;
			break;
		case channel::CLOSED:
			color = 1.0f;
			break;
		}
		channelsAttribs[i * 4 + 3] = color;
	}
}

void Simulation::keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	// exit
	if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
		glfwSetWindowShouldClose(window, GL_TRUE);

	Simulation* simulation = static_cast<Simulation*>(glfwGetWindowUserPointer(window));

	// F - freeze
	if (key == GLFW_KEY_F && action == GLFW_PRESS) {
		simulation->freeze();
		log(simulation->logfile, "[*] Simulation freezed = " + std::to_string(simulation->ice));
	}

	// X - reverse time for particles
	if (key == GLFW_KEY_X && action == GLFW_PRESS) {
		simulation->reverse();
		log(simulation->logfile, "[~] Simulation rewind = " + std::to_string(simulation->rewind));
	}

	// UP - camera speed + 20%
	if (key == GLFW_KEY_UP && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
		simulation->camera.movementSpeed *= 1.2f;
		log(simulation->logfile, "[+] Camera movement speed = " + std::to_string(simulation->camera.movementSpeed));
	}

	// DOWN - camera speed - 20%
	if (key == GLFW_KEY_DOWN && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
		simulation->camera.movementSpeed *= 0.8f;
		log(simulation->logfile, "[-] Camera movement speed = " + std::to_string(simulation->camera.movementSpeed));
	}

	// O - opacity -0.125
	if (key == GLFW_KEY_O && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
		if (simulation->uniforms.opacity >= 0.125f)
			simulation->uniforms.opacity -= 0.125f;
		log(simulation->logfile, "[-] Lipid bilayer opacity changed, opacity = " + std::to_string(simulation->uniforms.opacity));
	}

	// P - opacity +0.05
	if (key == GLFW_KEY_P && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
		if (simulation->uniforms.opacity < 1.0f)
			simulation->uniforms.opacity += 0.125f;
		log(simulation->logfile, "[+] Lipid bilayer opacity changed, opacity = " + std::to_string(simulation->uniforms.opacity));
	}
	
	// 1 - increase neurotransmitters count
	if (key == GLFW_KEY_1 && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
		simulation->increaseNeurotransmitters(10);
		log(simulation->logfile, "[+] Neurotransmitters count increased = " + std::to_string(simulation->activeParticlesCount[particle::NEUROTRANSMITTER]));
	}
	
	// 2 - decrease neurotransmitters count
	if (key == GLFW_KEY_2 && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
		simulation->decreaseNeurotransmitters(10);
		log(simulation->logfile, "[-] Neurotransmitters count decreased = " + std::to_string(simulation->activeParticlesCount[particle::NEUROTRANSMITTER]));
	}
	
	// 3 - decrease neurotransmitters count to 0
	if (key == GLFW_KEY_3 && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
		simulation->decreaseNeurotransmitters(simulation->activeParticlesCount[particle::NEUROTRANSMITTER]);
		log(simulation->logfile, "[-] Neurotransmitters count decreased = " + std::to_string(simulation->activeParticlesCount[particle::NEUROTRANSMITTER]));
	}

	// N - turn on/off rendering particles
	if (key == GLFW_KEY_N && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
		simulation->renderParticles = !simulation->renderParticles;
		log(simulation->logfile, "[*] Rendering particles = " + std::to_string(simulation->renderParticles));
	}
	
	// M - turn on/off rendering channels
	if (key == GLFW_KEY_M && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
		simulation->renderChannels = !simulation->renderChannels;
		log(simulation->logfile, "[*] Rendering channels = " + std::to_string(simulation->renderChannels));
	}

	// R - reset simulation
	if (key == GLFW_KEY_R && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
		simulation->reset();
		log(simulation->logfile, "[*] Simulation's been reset");
	}

	// TODO lock current simulation state before changing timeFactor
	// TODO better camera movement adjustment
	// RIGHT - time speed + 20%
	if (key == GLFW_KEY_RIGHT && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
		simulation->timeFactor *= 1.2f;
		simulation->inversedTimeFactor = 1.0f / simulation->timeFactor;
		std::ostringstream timeFactor;
		timeFactor << simulation->timeFactor;
		log(simulation->logfile, "[+] Time factor = " + timeFactor.str());

		// after changing time factor needs to change movement speed
		simulation->camera.movementSpeed *= 0.8333f;
		log(simulation->logfile, "[-] Camera movement speed = " + std::to_string(simulation->camera.movementSpeed));
	}

	// LEFT - time speed - 20%
	if (key == GLFW_KEY_LEFT && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
		simulation->timeFactor *= 0.8f;
		simulation->inversedTimeFactor = 1.0f / simulation->timeFactor;
		std::ostringstream timeFactor;
		timeFactor << simulation->timeFactor;
		log(simulation->logfile, "[-] Time factor = " + timeFactor.str());

		// after changing time factor needs to change movement speed
		simulation->camera.movementSpeed *= 1.25f;
		log(simulation->logfile, "[+] Camera movement speed = " + std::to_string(simulation->camera.movementSpeed));
	}

	if (key == GLFW_KEY_W && (action == GLFW_PRESS || action == GLFW_REPEAT))
		simulation->camera.processKeyboard(cam::FORWARD, simulation->getDeltaTime());

	if (key == GLFW_KEY_S && (action == GLFW_PRESS || action == GLFW_REPEAT))
		simulation->camera.processKeyboard(cam::BACKWARD, simulation->getDeltaTime());

	if (key == GLFW_KEY_A && (action == GLFW_PRESS || action == GLFW_REPEAT))
		simulation->camera.processKeyboard(cam::LEFT, simulation->getDeltaTime());

	if (key == GLFW_KEY_D && (action == GLFW_PRESS || action == GLFW_REPEAT))
		simulation->camera.processKeyboard(cam::RIGHT, simulation->getDeltaTime());
}

void Simulation::mouseButtonCallback(GLFWwindow* window, int button, int action, int mods)
{
	Simulation* simulation = static_cast<Simulation*>(glfwGetWindowUserPointer(window));

	if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {
		simulation->cursor = !simulation->cursor;
		if (simulation->cursor)
			glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
		else
			glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
	}
}

void Simulation::cursorPosCallback(GLFWwindow* window, double xpos, double ypos)
{
	Simulation* simulation = static_cast<Simulation*>(glfwGetWindowUserPointer(window));

	// if cursor is visible, don't callback to position changes
	if (simulation->cursor)
		return;

	float xoffset = (float)(xpos - simulation->camera.lastX);
	float yoffset = (float)(simulation->camera.lastY - ypos);

	simulation->camera.lastX = (float)xpos;
	simulation->camera.lastY = (float)ypos;

	simulation->camera.processMouseMovement(xoffset, yoffset);
}

void Simulation::scrollCallback(GLFWwindow * window, double xoffset, double yoffset)
{
	Simulation* simulation = static_cast<Simulation*>(glfwGetWindowUserPointer(window));
	simulation->camera.processMouseScroll((float)yoffset);
}

void Simulation::framebufferSizeCallback(GLFWwindow* window, int width, int height)
{
	Simulation* simulation = static_cast<Simulation*>(glfwGetWindowUserPointer(window));
	if (simulation->updateFramebufferSize(width, height))
		glViewport(0, 0, width, height);
}

bool Simulation::updateFramebufferSize(int width, int height)
{
	if (width <= 0 || height <= 0) {
		log(logfile, "Wrong framebuffer size! Width = " + std::to_string(width) + ", Height = " + std::to_string(height));
		return false;
	}
	this->width = width;
	this->height = height;
	return true;
}

inline void Simulation::updateChannelsStates()
{
	const int channelsBufferSizeInLoop = channelsBufferSize;

	// neurotransmitters
	const int neurotransmittersOffset = particlesOffset[particle::NEUROTRANSMITTER];
	const int activeNeurotransmittersInLoopIndex = particlesOffset[particle::NEUROTRANSMITTER] + activeParticlesCount[particle::NEUROTRANSMITTER];

	// nap ions
	const int NapIonsOffset = particlesOffset[particle::NAP];
	const int activeNapIonsInLoopIndex = particlesOffset[particle::NAP] + activeParticlesCount[particle::NAP];
	
	// kp ions
	const int KpIonsOffset = particlesOffset[particle::KP];
	const int activeKpIonsInLoopIndex = particlesOffset[particle::KP] + activeParticlesCount[particle::KP];

	// parallelization can be done due to no influence from other channels
#pragma loop(hint_parallel(0))
#pragma loop(ivdep)
	for (int i = 0; i < channelsBufferSizeInLoop; ++i) {
		Channel& currChannel = neuron->channels[i];
		
		// TODO probability to open instead of threshold (hidden markov model)
		// TODO add relative refraction
		// TODO channels open/close/inactive when deltaTime < 0

		// voltage gated channels need voltage inside neuron computing
		if (currChannel.gating == channel::VOLTAGE_GATED) {
			Particle* particle;
			float dx, dy, dz, d;
			float Ein = -0.065f;
			float Eout = 0.0f;
			float U;

			// calc voltage inside neuron but only from Nap ions inside neuron
			for (int j = NapIonsOffset; j < activeNapIonsInLoopIndex; ++j) {
				particle = &particles[bufferNum][j];

				dx = particle->x - currChannel.xIn;
				dy = particle->y - currChannel.yIn;
				dz = particle->z - currChannel.zIn;

				d = sqrt(dx * dx + dy * dy + dz * dz);
				if (d == 0.0)
					continue;

				Ein += phy::k * particle->charge / (metricFactor * d);
			}

			// calc voltage outside neuron but only from Kp ions outside neuron
			for (int j = KpIonsOffset; j < activeKpIonsInLoopIndex; ++j) {
				particle = &particles[bufferNum][j];

				dx = particle->x - currChannel.xOut;
				dy = particle->y - currChannel.yOut;
				dz = particle->z - currChannel.zOut;

				d = sqrt(dx * dx + dy * dy + dz * dz);
				if (d == 0.0)
					continue;

				Eout += phy::k * particle->charge / (metricFactor * d);
			}

			currChannel.U = U = Ein - Eout;
			
			if (currChannel.type == channel::NAP) {
				if (currChannel.state == channel::CLOSED) {
					if (U > phy::NapOpenTreshold) {
						currChannel.state = channel::OPEN;
						currChannel.timeLeft = phy::NapOpenTime;
						channelsAttribs[i * 4 + 3] = 0.0f;
					}
				}
				else if (currChannel.state == channel::OPEN) {
					currChannel.timeLeft -= fabs(deltaTime);
					if (currChannel.timeLeft < 0.0f) {
						currChannel.state = channel::INACTIVE;
						channelsAttribs[i * 4 + 3] = 0.5f;
					}
				}
				else if (currChannel.state == channel::INACTIVE) {
					if (U < phy::NapRepolarizationTreshold) {
						currChannel.state = channel::CLOSED;
						channelsAttribs[i * 4 + 3] = 1.0f;
					}
				}
			}
			else if (currChannel.type == channel::KP) {
				if (currChannel.state == channel::CLOSED) {
					if (U > phy::KpOpenTreshold) {
						currChannel.state = channel::OPEN;
						channelsAttribs[i * 4 + 3] = 0.0f;
					}
				}
				else if (currChannel.state == channel::OPEN) {
					if (U < phy::KpRepolarizationTreshold) {
						currChannel.state = channel::CLOSED;
						channelsAttribs[i * 4 + 3] = 1.0f;
					}
				}
			}
		}
		// ligand-gated channel's state opening in collisions calculating, need to be closed when no neurotransmitters
		else if (currChannel.gating == channel::LIGAND_GATED) {
			// open when neurotransmitter in its coordsOut position
			bool open = false;
			for (int k = neurotransmittersOffset; k < activeNeurotransmittersInLoopIndex && !open; ++k) {
				const Particle* neurotransmitter = &particles[bufferNum][k];
				// neurotransmitter nearby opening channel
				if (neurotransmitter->x == currChannel.xOut && neurotransmitter->y == currChannel.yOut && neurotransmitter->z == currChannel.zOut)
					open = true;
			}
			if (open) {
				currChannel.state = channel::OPEN;
				channelsAttribs[i * 4 + 3] = 0.0f;
			}
			else {
				currChannel.state = channel::CLOSED;
				channelsAttribs[i * 4 + 3] = 1.0f;
			}
		}
	}
}

inline void Simulation::calculateParticlesPositions()
{
	const unsigned short nextBufferNum = (bufferNum + 1) % 2;

	const int particlesTypes = particle::TYPES_COUNT;

	int offsets[particle::TYPES_COUNT];
#pragma loop(ivdep)
	for (int i = 0; i < particle::TYPES_COUNT; ++i)
		offsets[i] = particlesOffset[i];

	// for all particles types [ft] = first type
#pragma loop(ivdep)
	for (int ft = 0; ft < particlesTypes; ++ft) {
		const int firstParticlesOffset = offsets[ft];
		const int firstParticlesLastIndex = firstParticlesOffset + activeParticlesCount[ft];
		const double partAccOriginDt = partAccOrigin[ft] * deltaTime;

		// for all particle of one type
#pragma loop(hint_parallel(0))
#pragma loop(ivdep)
		for (int i = firstParticlesOffset; i < firstParticlesLastIndex; ++i) {
			Particle& currParticle = particles[bufferNum][i];
			Particle& prevParticle = particles[nextBufferNum][i];

			Particle* particle;
			float dx, dy, dz, dSq;
			float partAccOther;
			float ax = 0.0, ay = 0.0, az = 0.0;

			// again for all particles types [st] = second type
#pragma loop(ivdep)
			for (int st = 0; st < particlesTypes; ++st) {
				const int secondParticlesOffset = offsets[st];
				const int secondParticlesLastIndex = secondParticlesOffset + activeParticlesCount[st];

				// again for all particle of one type
#pragma loop(ivdep)
				for (int j = secondParticlesOffset; j < secondParticlesLastIndex; ++j) {
					particle = &particles[nextBufferNum][j];

					dx = particle->x - currParticle.x;
					dy = particle->y - currParticle.y;
					dz = particle->z - currParticle.z;

					// distance between 2 particles squared (d - distance, Sq - squared)
					dSq = dx * dx + dy * dy + dz * dz;
					if (dSq == 0.0)
						continue;

					// part of acceleration which depends on other particle's charge and distance
					partAccOther = particle->charge / dSq;

					ax -= partAccOther * dx;
					ay -= partAccOther * dy;
					az -= partAccOther * dz;
				}
			}

			// muliply by part of acceleration which depends on current particle
			double axdt = ax * partAccOriginDt;
			double aydt = ay * partAccOriginDt;
			double azdt = az * partAccOriginDt;

			// new coordinates
			// new coordinates: x = x0 + vx0 * dt + ax * dt^2 / 2
			// 'axdt = ax * dt' so optimized equation: x += dt * (vx0 + axdt / 2)
			currParticle.x = prevParticle.x + deltaTime * (currParticle.vx + axdt / 2);
			currParticle.y = prevParticle.y + deltaTime * (currParticle.vy + aydt / 2);
			currParticle.z = prevParticle.z + deltaTime * (currParticle.vz + azdt / 2);

			// update velocities: vx = vx0 + ax * dt
			currParticle.vx += axdt;
			currParticle.vy += aydt;
			currParticle.vz += azdt;
		}
	}
}

inline void Simulation::calculateCollisions()
{
	const unsigned short nextBufferNum = (bufferNum + 1) % 2;
	const int particlesTypes = particle::TYPES_COUNT;

	int offsets[particle::TYPES_COUNT];
	for (int i = 0; i < particle::TYPES_COUNT; ++i)
		offsets[i] = particlesOffset[i];

	// Nap and Kp ions
#pragma loop(hint_parallel(0))
#pragma loop(ivdep)
	for (int ft = 0; ft < particlesTypes; ++ft) {
		const int firstParticlesOffset = offsets[ft];
		const int firstParticlesLastIndex = firstParticlesOffset + activeParticlesCount[ft];

		for (int i = firstParticlesOffset; i < firstParticlesLastIndex; ++i) {
			Particle& currParticle = particles[bufferNum][i];
			Particle& prevParticle = particles[nextBufferNum][i];

			// check if collide with lipid bilayer and if, then check if bounced or pass through channel
			neuron->checkCollision(currParticle, prevParticle, (particle::Type)ft);
		}
	}
}

inline void Simulation::updateNapIonsFromChannels()
{
	const unsigned short nextBufferNum = (bufferNum + 1) % 2;

	if (activeParticlesCount[particle::NAP] == config.particlesCount[particle::NAP])
		return;

	for (int i = 0; i < channelsBufferSize; ++i) {
		Channel& currChannel = neuron->channels[i];

		if (currChannel.type == channel::NAP && currChannel.state == channel::OPEN) {
			float NapV0 = 5000;
			int ionsNum = 0;
			int iterations = 1;
			double NapIonsOccuredNum = config.particlesFlow[channel::NAP] * deltaTime;

			if (currChannel.gating == channel::LIGAND_GATED)
				iterations = 200;

			do {
				if (NapIonsOccuredNum < 1.0) {
					if (getRandDouble(0.0, 1.0) < NapIonsOccuredNum)
						++ionsNum;
				}
				else
					ionsNum = NapIonsOccuredNum;

				--iterations;
			} while (iterations > 0);

			for (int j = 0; j < ionsNum && activeParticlesCount[particle::NAP] < config.particlesCount[particle::NAP]; ++j) {
				float NapCoord0 = getRandDouble(0.001, 0.002);
				Particle& currParticle = particles[bufferNum][particlesOffset[particle::NAP] + activeParticlesCount[particle::NAP]];
				Particle& prevParticle = particles[nextBufferNum][particlesOffset[particle::NAP] + activeParticlesCount[particle::NAP]];

				++activeParticlesCount[particle::NAP];
				glm::vec3 n = glm::normalize(glm::vec3(currChannel.xIn - currChannel.xOut, currChannel.yIn - currChannel.yOut, currChannel.zIn - currChannel.zOut));
				glm::vec3 v = n * NapV0;
				glm::vec3 coords = n * NapCoord0;

				prevParticle.x = currParticle.x = currChannel.xIn + coords[0];
				prevParticle.y = currParticle.y = currChannel.yIn + coords[1];
				prevParticle.z = currParticle.z = currChannel.zIn + coords[2];

				prevParticle.vx = currParticle.vx = v[0];
				prevParticle.vy = currParticle.vy = v[1];
				prevParticle.vz = currParticle.vz = v[2];
			}
		}
	}
}

inline void Simulation::updateKpIonsFromChannels()
{
	const unsigned short nextBufferNum = (bufferNum + 1) % 2;
	const int channelsBufferSizeInLoop = channelsBufferSize;

	if (activeParticlesCount[particle::KP] == config.particlesCount[particle::KP])
		return;

	for (int i = 0; i < channelsBufferSizeInLoop; ++i) {
		Channel& currChannel = neuron->channels[i];

		if (currChannel.type == channel::KP && currChannel.state == channel::OPEN) {
			float KpV0 = 5000;
			int ionsNum = 0;
			int iterations = 1;
			double NapIonsOccuredNum = config.particlesFlow[channel::KP] * deltaTime;

			do {
				if (NapIonsOccuredNum < 1.0) {
					if (getRandDouble(0.0, 1.0) < NapIonsOccuredNum)
						++ionsNum;
				}
				else
					ionsNum = NapIonsOccuredNum;

				--iterations;
			} while (iterations > 0);

			for (int j = 0; j < ionsNum && activeParticlesCount[particle::KP] < config.particlesCount[particle::KP]; ++j) {
				float KpCoord0 = getRandDouble(0.001, 0.002);
				Particle& currParticle = particles[bufferNum][particlesOffset[particle::KP] + activeParticlesCount[particle::KP]];
				Particle& prevParticle = particles[nextBufferNum][particlesOffset[particle::KP] + activeParticlesCount[particle::KP]];

				++activeParticlesCount[particle::KP];
				glm::vec3 n = glm::normalize(glm::vec3(currChannel.xOut - currChannel.xIn, currChannel.yOut - currChannel.yIn, currChannel.zOut - currChannel.zIn));
				glm::vec3 v = n * KpV0;
				glm::vec3 coords = n * KpCoord0;

				prevParticle.x = currParticle.x = currChannel.xOut + coords[0];
				prevParticle.y = currParticle.y = currChannel.yOut + coords[1];
				prevParticle.z = currParticle.z = currChannel.zOut + coords[2];

				prevParticle.vx = currParticle.vx = v[0];
				prevParticle.vy = currParticle.vy = v[1];
				prevParticle.vz = currParticle.vz = v[2];
			}
		}
	}
}

inline void Simulation::updateParticlesPositions()
{
	const unsigned short nextBufferNum = (bufferNum + 1) % 2;

	const int particlesTypes = particle::TYPES_COUNT;

	int offsets[particle::TYPES_COUNT];
	for (int i = 0; i < particle::TYPES_COUNT; ++i)
		offsets[i] = particlesOffset[i];

	for (int ft = 0; ft < particlesTypes; ++ft) {
		const int particlesOffset = offsets[ft];
		const int particlesLastIndex = particlesOffset + activeParticlesCount[ft];
		// Nap and Kp ions
#pragma loop(hint_parallel(0))
#pragma loop(ivdep)
		for (int i = particlesOffset; i < particlesLastIndex; ++i) {
			Particle& currParticle = particles[bufferNum][i];
			Particle& prevParticle = particles[nextBufferNum][i];

			// update positions, all particles have 3 coords, x, y, z, that's why 'index * 3' in particlesPos[]
			particlesPos[i * 3 + 0] = currParticle.x;
			particlesPos[i * 3 + 1] = currParticle.y;
			particlesPos[i * 3 + 2] = currParticle.z;
		}
	}

}

inline void Simulation::update()
{
	if (ice)
		return;

	// update channles states
	updateChannelsStates();

	// calculate particles positions based on eletric forces
	calculateParticlesPositions();

	// calculate particles positions when collide
	calculateCollisions();

	// update Nap ions inflow from open channels
	updateNapIonsFromChannels();

	// update Kp ions outflow from open channels
	updateKpIonsFromChannels();

	// update particles positions buffer with their calculated positions
	updateParticlesPositions();

	// swap to next particles vector
	++bufferNum %= 2;
}

inline void Simulation::render()
{
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	uniforms.viewMatrix = glm::perspective(glm::radians(camera.zoom), (float)width / (float)height, 0.1f, 100.0f) * camera.GetViewMatrix();

	if (renderChannels) {
		// render channels
		channelsRenderProgram->use();

		for (int type = 0; type < channel::TYPES_COUNT; ++type) {
			channel::Type currentType = (channel::Type)type;
			uniforms.channelRadius = channelRadius[currentType];
			channelsRenderProgram->setUniforms(uniforms);
			glBindTexture(GL_TEXTURE_2D, channelTexture[currentType]);
			glBindVertexArray(channelVAO[currentType]);
			glDrawArrays(GL_POINTS, channelsOffset[currentType], activeChannelsCount[currentType]);
		}
	}
	
	if (renderParticles) {
		// render particles
		ionsRenderProgram->use();

		for (int type = 0; type < particle::TYPES_COUNT; ++type) {
			particle::Type currentType = (particle::Type)type;
			uniforms.ionRadius = particleRadius[currentType];
			ionsRenderProgram->setUniforms(uniforms);
			glBindTexture(GL_TEXTURE_2D, particleTexture[currentType]);
			glBindVertexArray(particleVAO[currentType]);
			glDrawArrays(GL_POINTS, particlesOffset[currentType], activeParticlesCount[currentType]);
		}
	}
	
	// render neuron
	neuron->render(uniforms);

	glFlush();
	glfwSwapBuffers(window);
	glBindVertexArray(0);
	glUseProgram(0);
	glfwPollEvents();
}

Simulation::Simulation(const Config& config)
{
	log(logfile, "--- Simulation initialization has started ---");

	log(logfile, "--- Loading configuration... ---");
	loadConfig(config);
	log(logfile, "--- Completed ---");

	log(logfile, "--- Setup OpenGL context... ---");
	setupOpenGL();
	log(logfile, "--- Completed ---");

	log(logfile, "--- Setup GLFW input methods... ---");
	setupInput();
	log(logfile, "--- Completed ---");

	log(logfile, "--- Setup shader programs... ---");
	setupPrograms();
	log(logfile, "--- Completed ---");

	log(logfile, "--- Setup simulation structures... ---");
	setupStructures();
	log(logfile, "--- Completed ---");

	log(logfile, "--- Setup uniforms... ---");
	setupUniforms();
	log(logfile, "--- Completed ---");

	log(logfile, "--- Setup textures... ---");
	setupTextures();
	log(logfile, "--- Completed ---");

	log(logfile, "--- Setup particles and channels buffers... ---");
	setupBuffers();
	log(logfile, "--- Completed ---");

	log(logfile, "--- Simulation initialization's ended successfully ---");
}

Simulation::~Simulation()
{
	log(logfile, "--- Simulation terminating has started... ---");
	glfwTerminate();
	delete neuron;
	delete ionsRenderProgram;
	delete channelsRenderProgram;
	log(logfile, "--- Simulation ends ---");
	logfile.close();
}

void Simulation::start(void)
{
	log(logfile, "--- Simulation's started ---");
	currentFrame = lastFrame = deltaTime = 0.0;

	// main simulation loop
	while (!glfwWindowShouldClose(window))
	{
		currentFrame = glfwGetTime();
		deltaTime = timeFactor * (currentFrame - lastFrame);
		lastFrame = currentFrame;

		// simulation logic
		update();

		// rendering
		render();
	}
}

void Simulation::freeze(void)
{
	ice = !ice;
}

void Simulation::reverse(void)
{
	rewind = !rewind;
	timeFactor = -timeFactor;
}

void Simulation::reset(void)
{
	for (int i = 0; i < particle::TYPES_COUNT; ++i)
		activeParticlesCount[i] = 0;
	bufferNum = 0;
	const long neurotransmittersOffset = particlesBufferSize - config.particlesCount[particle::NEUROTRANSMITTER];
	for (long i = 0; i < channelsBufferSize; ++i) {
		Channel& channel = neuron->channels[i];
		channel.state = channel::CLOSED;
		channelsAttribs[4 * i + 3] = 1.0f;
	}
	for (long i = 0; i < config.particlesCount[particle::NEUROTRANSMITTER]; ++i) {
		float offset[3] = { getRandDouble(-0.05, -0.01), getRandDouble(-0.03, 0.03), getRandDouble(-0.03, 0.03) };
		float coords[3] = { synapsePosition[0] + offset[0], synapsePosition[1] + offset[1], synapsePosition[2] + offset[2] };
		float velocities[3] = { 10000.0f, 0.0f, 0.0f };
		
		particles[0][neurotransmittersOffset + i].x = particles[1][neurotransmittersOffset + i].x = coords[0];
		particles[0][neurotransmittersOffset + i].y = particles[1][neurotransmittersOffset + i].y = coords[1];
		particles[0][neurotransmittersOffset + i].z = particles[1][neurotransmittersOffset + i].z = coords[2];

		particles[0][neurotransmittersOffset + i].vx = particles[1][neurotransmittersOffset + i].vx = velocities[0];
		particles[0][neurotransmittersOffset + i].vy = particles[1][neurotransmittersOffset + i].vy = velocities[1];
		particles[0][neurotransmittersOffset + i].vz = particles[1][neurotransmittersOffset + i].vz = velocities[2];
	}
}

void Simulation::decreaseNeurotransmitters(const unsigned n)
{
	if (n > activeParticlesCount[particle::NEUROTRANSMITTER])
		activeParticlesCount[particle::NEUROTRANSMITTER] = 0;
	else
		activeParticlesCount[particle::NEUROTRANSMITTER] -= n;
}

void Simulation::increaseNeurotransmitters(const unsigned n)
{
	const long neurotransmittersOffset = particlesBufferSize - config.particlesCount[particle::NEUROTRANSMITTER] + activeParticlesCount[particle::NEUROTRANSMITTER];

	if (activeParticlesCount[particle::NEUROTRANSMITTER] < config.particlesCount[particle::NEUROTRANSMITTER] - n) {
		activeParticlesCount[particle::NEUROTRANSMITTER] += n;

		for (long i = 0; i < n; ++i) {
			float offset[3] = { -phy::synapticGapWidth / metricFactor, getRandDouble(-0.03, 0.03), getRandDouble(-0.03, 0.03) };
			float coords[3] = { synapsePosition[0] + offset[0], synapsePosition[1] + offset[1], synapsePosition[2] + offset[2] };
			float velocities[3] = { 100000.0f, 0.0f, 0.0f };

			particles[0][neurotransmittersOffset + i].x = particles[1][neurotransmittersOffset + i].x = coords[0];
			particles[0][neurotransmittersOffset + i].y = particles[1][neurotransmittersOffset + i].y = coords[1];
			particles[0][neurotransmittersOffset + i].z = particles[1][neurotransmittersOffset + i].z = coords[2];

			particles[0][neurotransmittersOffset + i].vx = particles[1][neurotransmittersOffset + i].vx = velocities[0];
			particles[0][neurotransmittersOffset + i].vy = particles[1][neurotransmittersOffset + i].vy = velocities[1];
			particles[0][neurotransmittersOffset + i].vz = particles[1][neurotransmittersOffset + i].vz = velocities[2];
		}
	}
}

double Simulation::getDeltaTime(void) const
{
	return inversedTimeFactor * fabs(deltaTime);
}
