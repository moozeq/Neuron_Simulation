#include "Simulation.h"

void Simulation::loadConfig(const Config& _config)
{
	// copy config struct
	config = _config;

	// set variables based on config struct
	metricFactorSq = config.metricFactor * config.metricFactor;
	metricFactor = config.metricFactor;
	timeFactor = config.timeFactor;
	bufferNum = 0;

	inversedTimeFactor = 1.0 / config.timeFactor;
	particlesBufferSize = 0;
	for (int i = 0; i < particle::TYPES_COUNT; ++i) {
		particlesOffsets[i] = particlesBufferSize;
		particlesBufferSize += config.particlesCount[i];
	}

	activeParticlesCount = 0;
	activeNeurotransmittersCount = 0;
	activeNapsCount = 0;
	activeKpsCount = 0;

	ice = false;
	rewind = false;
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
	log(logfile, "--- Simulation initialization's started ---");

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

	log(logfile, "--- Simulation initialization's ended successfully ---");
}

void Simulation::setupInput()
{
	// input functions
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
	neuron = new Neuron(metricFactor, timeFactor, config);
	channelsBufferSize = neuron->channels.size();
	NapChannelsCount = neuron->allNapChannelsCount;
	KpChannelsCount = neuron->allKpChannelsCount;
	synapsePosition = neuron->getSynapsePosition();
}

void Simulation::setupParticlesStructures()
{
	particles[0].reserve(particlesBufferSize);
	particles[1].reserve(particlesBufferSize);
	partAccOrigin.reserve(particlesBufferSize);
	size_t i = 0;
	size_t offset = 0;
	double boundaries[3][2] = { { -1.4f, 12.4f }, { -0.03f, 0.03f }, { -0.03f, 0.03f } };

	offset += config.particlesCount[particle::NAP];
	while (i < offset) {
		Particle* particle = newParticle(boundaries, particle::NAP);
		particles[0].push_back(*particle);
		particles[1].push_back(*particle);
		partAccOrigin.push_back(phy::NapA / metricFactorSq);
		delete particle;
		++i;
	}

	offset += config.particlesCount[particle::KP];
	while (i < offset) {
		Particle* particle = newParticle(boundaries, particle::KP);
		particles[0].push_back(*particle);
		particles[1].push_back(*particle);
		partAccOrigin.push_back(phy::KpA / metricFactorSq);
		delete particle;
		++i;
	}

	offset += config.particlesCount[particle::CLM];
	while (i < offset) {
		Particle* particle = newParticle(boundaries, particle::CLM);
		particles[0].push_back(*particle);
		particles[1].push_back(*particle);
		partAccOrigin.push_back(phy::ClmA / metricFactorSq);
		delete particle;
		++i;
	}

	offset += config.particlesCount[particle::ORGANIC_ANION];
	while (i < offset) {
		Particle* particle = newParticle(boundaries, particle::ORGANIC_ANION);
		particles[0].push_back(*particle);
		particles[1].push_back(*particle);
		partAccOrigin.push_back((phy::k * particle->charge / particle->mass) / metricFactorSq);
		delete particle;
		++i;
	}

	offset += config.particlesCount[particle::NEUROTRANSMITTER];
	while (i < offset) {
		float offset[3] = { getRandDouble(-0.05, -0.01), getRandDouble(-0.003, 0.003), getRandDouble(-0.003, 0.003) };
		float coords[3] = { synapsePosition[0] + offset[0], synapsePosition[1] + offset[1], synapsePosition[2] + offset[2] };
		float velocities[3] = { 10000.0f, 0.0f, 0.0f };
		Particle* particle = newParticle(coords, velocities, particle::NEUROTRANSMITTER);
		particles[0].push_back(*particle);
		particles[1].push_back(*particle);
		partAccOrigin.push_back((phy::k * particle->charge / particle->mass) / metricFactorSq);
		delete particle;
		++i;
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
	NapIonTexture = loadMipmapTexture(GL_TEXTURE0, config.particlesTextures[particle::NAP].c_str());
	KpIonTexture = loadMipmapTexture(GL_TEXTURE0, config.particlesTextures[particle::KP].c_str());
	ClmIonTexture = loadMipmapTexture(GL_TEXTURE0, config.particlesTextures[particle::CLM].c_str());
	otherParticlesTexture = loadMipmapTexture(GL_TEXTURE0, config.particlesTextures[particle::ORGANIC_ANION].c_str());
	neurotransmittersTexture = loadMipmapTexture(GL_TEXTURE0, config.particlesTextures[particle::NEUROTRANSMITTER].c_str());

	// load channels textures
	NapIonChannelTexture = loadMipmapTexture(GL_TEXTURE0, config.NapChannelTexturePath.c_str());
	KpIonChannelTexture = loadMipmapTexture(GL_TEXTURE0, config.KpChannelTexturePath.c_str());
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
	glCreateVertexArrays(1, &NapIonsVAO);
	glCreateVertexArrays(1, &KpIonsVAO);
	glCreateVertexArrays(1, &ClmIonsVAO);
	glCreateVertexArrays(1, &otherParticlesVAO);
	glCreateVertexArrays(1, &neurotransmittersVAO);

	// creating buffers
	glCreateBuffers(1, &particlesPosBuf);

	GLbitfield flags = GL_MAP_WRITE_BIT | GL_MAP_PERSISTENT_BIT | GL_MAP_COHERENT_BIT;
	GLsizei bufferSize = particlesBufferSize * 3 * sizeof(GLfloat);
	glNamedBufferStorage(particlesPosBuf, bufferSize, nullptr, flags);

	glVertexArrayVertexBuffer(NapIonsVAO, 0, particlesPosBuf, 0, 3 * sizeof(GLfloat));
	glVertexArrayAttribFormat(NapIonsVAO, 0, 3, GL_FLOAT, GL_FALSE, 0);
	glVertexArrayAttribBinding(NapIonsVAO, 0, 0);
	glEnableVertexArrayAttrib(NapIonsVAO, 0);

	glVertexArrayVertexBuffer(KpIonsVAO, 0, particlesPosBuf, 0, 3 * sizeof(GLfloat));
	glVertexArrayAttribFormat(KpIonsVAO, 0, 3, GL_FLOAT, GL_FALSE, 0);
	glVertexArrayAttribBinding(KpIonsVAO, 0, 0);
	glEnableVertexArrayAttrib(KpIonsVAO, 0);

	glVertexArrayVertexBuffer(ClmIonsVAO, 0, particlesPosBuf, 0, 3 * sizeof(GLfloat));
	glVertexArrayAttribFormat(ClmIonsVAO, 0, 3, GL_FLOAT, GL_FALSE, 0);
	glVertexArrayAttribBinding(ClmIonsVAO, 0, 0);
	glEnableVertexArrayAttrib(ClmIonsVAO, 0);

	glVertexArrayVertexBuffer(otherParticlesVAO, 0, particlesPosBuf, 0, 3 * sizeof(GLfloat));
	glVertexArrayAttribFormat(otherParticlesVAO, 0, 3, GL_FLOAT, GL_FALSE, 0);
	glVertexArrayAttribBinding(otherParticlesVAO, 0, 0);
	glEnableVertexArrayAttrib(otherParticlesVAO, 0);

	glVertexArrayVertexBuffer(neurotransmittersVAO, 0, particlesPosBuf, 0, 3 * sizeof(GLfloat));
	glVertexArrayAttribFormat(neurotransmittersVAO, 0, 3, GL_FLOAT, GL_FALSE, 0);
	glVertexArrayAttribBinding(neurotransmittersVAO, 0, 0);
	glEnableVertexArrayAttrib(neurotransmittersVAO, 0);

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
	glCreateVertexArrays(1, &NapIonsChannelsVAO);
	glCreateVertexArrays(1, &KpIonsChannelsVAO);

	// creating buffers
	glCreateBuffers(1, &channelsBuf);

	GLbitfield flags = GL_MAP_WRITE_BIT | GL_MAP_PERSISTENT_BIT | GL_MAP_COHERENT_BIT;
	std::vector<float> channels = neuron->getChannels();

	glNamedBufferStorage(channelsBuf, channelsBufferSize * 4 * sizeof(GLfloat), &channels[0], flags);

	// positions
	glVertexArrayVertexBuffer(NapIonsChannelsVAO, 0, channelsBuf, 0, 4 * sizeof(GLfloat));
	glVertexArrayAttribFormat(NapIonsChannelsVAO, 0, 4, GL_FLOAT, GL_FALSE, 0);
	glVertexArrayAttribBinding(NapIonsChannelsVAO, 0, 0);
	glEnableVertexArrayAttrib(NapIonsChannelsVAO, 0);

	glVertexArrayVertexBuffer(KpIonsChannelsVAO, 0, channelsBuf, 0, 4 * sizeof(GLfloat));
	glVertexArrayAttribFormat(KpIonsChannelsVAO, 0, 4, GL_FLOAT, GL_FALSE, 0);
	glVertexArrayAttribBinding(KpIonsChannelsVAO, 0, 0);
	glEnableVertexArrayAttrib(KpIonsChannelsVAO, 0);

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

	// R - reverse
	if (key == GLFW_KEY_R && action == GLFW_PRESS) {
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

	// O - opacity -0.05
	if (key == GLFW_KEY_O && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
		if (simulation->uniforms.opacity > 0.05f)
			simulation->uniforms.opacity -= 0.05f;
		log(simulation->logfile, "[-] Lipid bilayer opacity changed, opacity = " + std::to_string(simulation->uniforms.opacity));
	}

	// P - opacity +0.05
	if (key == GLFW_KEY_P && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
		if (simulation->uniforms.opacity < 1.0f)
			simulation->uniforms.opacity += 0.05f;
		log(simulation->logfile, "[+] Lipid bilayer opacity changed, opacity = " + std::to_string(simulation->uniforms.opacity));
	}
	
	// 1 - increase neurotransmitters count
	if (key == GLFW_KEY_1 && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
		simulation->increaseNeurotransmitters(10);
		log(simulation->logfile, "[+] Neurotransmitters count increased = " + std::to_string(simulation->activeNeurotransmittersCount));
	}
	
	// 2 - decrease neurotransmitters count
	if (key == GLFW_KEY_2 && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
		simulation->decreaseNeurotransmitters(10);
		log(simulation->logfile, "[-] Neurotransmitters count decreased = " + std::to_string(simulation->activeNeurotransmittersCount));
	}
	
	// 3 - decrease neurotransmitters count to 0
	if (key == GLFW_KEY_3 && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
		simulation->decreaseNeurotransmitters(simulation->activeNeurotransmittersCount);
		log(simulation->logfile, "[-] Neurotransmitters count decreased = " + std::to_string(simulation->activeNeurotransmittersCount));
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

	// E - reset simulation
	if (key == GLFW_KEY_E && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
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
	if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {
		double x, y;
		glfwGetCursorPos(window, &x, &y);

		Simulation* simulation = static_cast<Simulation*>(glfwGetWindowUserPointer(window));

		// translate from [0, WIDTH/HEIGHT] to [-1.0, 1.0]
		x /= simulation->width / 2;
		y /= simulation->height / 2;
		x -= 1.0;
		y -= 1.0;
		y *= -1.0;
	}
}

void Simulation::cursorPosCallback(GLFWwindow* window, double xpos, double ypos)
{
	Simulation* simulation = static_cast<Simulation*>(glfwGetWindowUserPointer(window));

	float xoffset = (float)(xpos - simulation->camera.lastX);
	float yoffset = (float)(simulation->camera.lastY - ypos);

	simulation->camera.lastX = (float)xpos;
	simulation->camera.lastY = (float)ypos;

	int state = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT);
	if (state == GLFW_PRESS)
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
	const long channelsBufferSizeInLoop = channelsBufferSize;

	// neurotransmitters from first to last
	const long neurotransmittersOffset = particlesOffsets[particle::NEUROTRANSMITTER];
	const long activeNeurotransmittersInLoopIndex = particlesOffsets[particle::NEUROTRANSMITTER] + activeNeurotransmittersCount;

	// nap ions counted from last to first
	const long NapIonsOffset = config.particlesCount[particle::NAP] - activeNapsCount;
	const long activeNapIonsInLoopIndex = config.particlesCount[particle::NAP];
	
	// kp ions counted from first to last
	const long KpIonsOffset = particlesOffsets[particle::KP];
	const long activeKpIonsInLoopIndex = particlesOffsets[particle::KP] + activeKpsCount;
	
	const long particlesOffset = config.particlesCount[particle::NAP] - activeNapsCount;
	const long particlesBufferSizeInLoopIndex = particlesOffset + activeNapsCount + activeKpsCount;

	// parallelization can be done due to no influence from other channels
#pragma loop(hint_parallel(0))
#pragma loop(ivdep)
	for (long i = 0; i < channelsBufferSizeInLoop; ++i) {
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
			for (long j = NapIonsOffset; j < activeNapIonsInLoopIndex; ++j) {
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
			for (long j = KpIonsOffset; j < activeKpIonsInLoopIndex; ++j) {
				particle = &particles[bufferNum][j];

				dx = particle->x - currChannel.xOut;
				dy = particle->y - currChannel.yOut;
				dz = particle->z - currChannel.zOut;

				d = sqrt(dx * dx + dy * dy + dz * dz);
				if (d == 0.0)
					continue;

				Ein -= phy::k * particle->charge / (metricFactor * d);
			}

			currChannel.U = U = Ein - Eout;
			if (i < 10)
				std::cout << "\n[Channel] index = " + std::to_string(i) + " Ein = " + std::to_string(Ein) + ", Eout = " + std::to_string(Eout) + ", U = " + std::to_string(U);

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
						currChannel.timeLeft = phy::NapOpenTime;
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
			if (currChannel.state == channel::CLOSED)
				channelsAttribs[i * 4 + 3] = 1.0f;
			else if (currChannel.state == channel::OPEN) {
				// open when neurotransmitter in its coordsOut position
				bool open = false;
				for (long k = neurotransmittersOffset; k < activeNeurotransmittersInLoopIndex && !open; ++k) {
					const Particle* neurotransmitter = &particles[bufferNum][k];
					if (neurotransmitter->x == currChannel.xOut && neurotransmitter->y == currChannel.yOut && neurotransmitter->z == currChannel.zOut)
						open = true;
				}
				if (!open) {
					currChannel.state = channel::CLOSED;
					channelsAttribs[i * 4 + 3] = 1.0f;
				}
				else
					channelsAttribs[i * 4 + 3] = 0.0f;
			}
		}
	}
}

inline void Simulation::calculateParticlesPositions()
{
	const unsigned short nextBufferNum = (bufferNum + 1) % 2;

	// neurotransmitters from first to last
	const long neurotransmittersOffset = particlesOffsets[particle::NEUROTRANSMITTER];
	const long activeNeurotransmittersInLoopIndex = particlesOffsets[particle::NEUROTRANSMITTER] + activeNeurotransmittersCount;

	const long particlesOffset = config.particlesCount[particle::NAP] - activeNapsCount;
	const long particlesBufferSizeInLoopIndex = particlesOffset + activeNapsCount + activeKpsCount;

	// parallelization can be done due to double buffering particles vector
#pragma loop(hint_parallel(0))
#pragma loop(ivdep)
	for (long i = particlesOffset; i < particlesBufferSizeInLoopIndex; ++i) {
		Particle& currParticle = particles[bufferNum][i];
		Particle& prevParticle = particles[nextBufferNum][i];

		//const size_t index = currParticle.index;
		const double partAccOriginDt = partAccOrigin[i] * deltaTime;

		Particle* particle;
		float dx, dy, dz, dSq;
		float partAccOther;
		float ax = 0.0, ay = 0.0, az = 0.0;

		// caluclating forces with omiting checking if (i == j) with for loops: [0; current - 1] and [current + 1; last]
		// calculate forces from particles 0 to current - 1
		for (long j = particlesOffset; j < i; ++j) {
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

		// calculate forces from particles current + 1 to last
		for (long j = i + 1; j < particlesBufferSizeInLoopIndex; ++j) {
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

	// parallelization can be done due to double buffering particles vector
#pragma loop(hint_parallel(0))
#pragma loop(ivdep)
	for (long i = neurotransmittersOffset; i < activeNeurotransmittersInLoopIndex; ++i) {
		Particle& currParticle = particles[bufferNum][i];
		Particle& prevParticle = particles[nextBufferNum][i];

		//const size_t index = currParticle.index;
		const double partAccOriginDt = partAccOrigin[i] * deltaTime;

		Particle* particle;
		float dx, dy, dz, dSq;
		float partAccOther;
		float ax = 0.0, ay = 0.0, az = 0.0;

		// caluclating forces with omiting checking if (i == j) with for loops: [0; current - 1] and [current + 1; last]
		// calculate forces from particles 0 to current - 1
		for (long j = neurotransmittersOffset; j < i; ++j) {
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

		// calculate forces from particles current + 1 to last
		for (long j = i + 1; j < activeNeurotransmittersInLoopIndex; ++j) {
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

inline void Simulation::calculateCollisions()
{
	const unsigned short nextBufferNum = (bufferNum + 1) % 2;
	const long neurotransmittersOffset = particlesOffsets[particle::NEUROTRANSMITTER];
	const long activeNeurotransmittersCountInLoopIndex = particlesOffsets[particle::NEUROTRANSMITTER] + activeNeurotransmittersCount;
	const long particlesOffset = config.particlesCount[particle::NAP] - activeNapsCount;
	const long particlesBufferSizeInLoopIndex = particlesOffset + activeNapsCount + activeKpsCount;
	
	long offsets[particle::TYPES_COUNT];
	for (int i = 0; i < particle::TYPES_COUNT; ++i)
		offsets[i] = particlesOffsets[i];

	// Nap and Kp ions
#pragma loop(hint_parallel(0))
#pragma loop(ivdep)
	for (long i = particlesOffset; i < particlesBufferSizeInLoopIndex; ++i) {
		Particle& currParticle = particles[bufferNum][i];
		Particle& prevParticle = particles[nextBufferNum][i];

		// check if collide with lipid bilayer and if, then check if bounced or pass through channel
		if (i < offsets[particle::NAP])
			neuron->checkCollision(currParticle, prevParticle, particle::NAP);
		else if (i < offsets[particle::KP])
			neuron->checkCollision(currParticle, prevParticle, particle::KP);
		else if (i < offsets[particle::CLM])
			neuron->checkCollision(currParticle, prevParticle, particle::CLM);
		else if (i < offsets[particle::ORGANIC_ANION])
			neuron->checkCollision(currParticle, prevParticle, particle::ORGANIC_ANION);
	}

	// neurotransmitters
#pragma loop(hint_parallel(0))
#pragma loop(ivdep)
	for (long i = neurotransmittersOffset; i < activeNeurotransmittersCountInLoopIndex; ++i) {
		Particle& currParticle = particles[bufferNum][i];
		Particle& prevParticle = particles[nextBufferNum][i];

		// check if collide with lipid bilayer and if, then check if bounced or pass through channel
		neuron->checkCollision(currParticle, prevParticle, particle::NEUROTRANSMITTER);
	}
}

inline void Simulation::updateNapIonsFromChannels()
{
	const unsigned short nextBufferNum = (bufferNum + 1) % 2;
	const long channelsBufferSizeInLoop = channelsBufferSize;

	if (activeNapsCount == config.particlesCount[particle::NAP])
		return;

	for (long i = 0; i < channelsBufferSizeInLoop; ++i) {
		Channel& currChannel = neuron->channels[i];

		if (currChannel.type == channel::NAP && currChannel.state == channel::OPEN) {
			float NapV0 = 5000;
			unsigned ionsNum = 0;// timeFactor * 1e7 * getRandDouble(0.0, config.NapInflow);

			if (currChannel.gating == channel::LIGAND_GATED) {
				for (long k = 0; k < 20; ++k) {
					if (getRandDouble(0.0, 1.0) < 0.005 * (1.0 / config.timeFactor) * timeFactor)
						++ionsNum;
				}
			}

			if (getRandDouble(0.0, 1.0) < 0.005 * (1.0 / config.timeFactor) * timeFactor)
				++ionsNum;

			for (long j = 0; j < ionsNum; ++j) {
				float NapCoord0 = getRandDouble(0.001, 0.002);
				++activeParticlesCount;
				++activeNapsCount;
				Particle& currParticle = particles[bufferNum][config.particlesCount[particle::NAP] - activeNapsCount];
				Particle& prevParticle = particles[nextBufferNum][config.particlesCount[particle::NAP] - activeNapsCount];

				glm::vec3 n = glm::normalize(glm::vec3(currChannel.xIn - currChannel.xOut, currChannel.yIn - currChannel.yOut, currChannel.zIn - currChannel.zOut));
				glm::vec3 v = n * NapV0;
				glm::vec3 coords = n * NapCoord0;

				prevParticle.x = currParticle.x = currChannel.xIn + coords[0];
				prevParticle.y = currParticle.y = currChannel.yIn + coords[1];
				prevParticle.z = currParticle.z = currChannel.zIn + coords[2];

				prevParticle.vx = currParticle.vx = v[0];
				prevParticle.vy = currParticle.vy = v[1];
				prevParticle.vz = currParticle.vz = v[2];

				if (activeNapsCount == config.particlesCount[particle::NAP])
					return;
			}
		}
	}
}

inline void Simulation::updateKpIonsFromChannels()
{
	const unsigned short nextBufferNum = (bufferNum + 1) % 2;
	const long channelsBufferSizeInLoop = channelsBufferSize;

	if (activeKpsCount == config.particlesCount[particle::KP])
		return;

	for (long i = 0; i < channelsBufferSizeInLoop; ++i) {
		Channel& currChannel = neuron->channels[i];

		if (currChannel.type == channel::KP && currChannel.state == channel::OPEN) {
			float KpV0 = 5000;
			unsigned ionsNum = 0;// timeFactor * 1e7 * getRandDouble(0.0, config.NapInflow);

			if (getRandDouble(0.0, 1.0) < 0.01 * (1.0 / config.timeFactor) * timeFactor)
				++ionsNum;

			for (long j = 0; j < ionsNum; ++j) {
				float KpCoord0 = getRandDouble(0.001, 0.002);
				++activeParticlesCount;
				++activeKpsCount;
				Particle& currParticle = particles[bufferNum][particlesOffsets[particle::KP] + activeKpsCount];
				Particle& prevParticle = particles[nextBufferNum][particlesOffsets[particle::KP] + activeKpsCount];

				glm::vec3 n = glm::normalize(glm::vec3(currChannel.xOut - currChannel.xIn, currChannel.yOut - currChannel.yIn, currChannel.zOut - currChannel.zIn));
				glm::vec3 v = n * KpV0;
				glm::vec3 coords = n * KpCoord0;

				prevParticle.x = currParticle.x = currChannel.xOut + coords[0];
				prevParticle.y = currParticle.y = currChannel.yOut + coords[1];
				prevParticle.z = currParticle.z = currChannel.zOut + coords[2];

				prevParticle.vx = currParticle.vx = v[0];
				prevParticle.vy = currParticle.vy = v[1];
				prevParticle.vz = currParticle.vz = v[2];

				if (activeKpsCount == config.particlesCount[particle::KP])
					return;
			}
		}
	}
}

inline void Simulation::updateParticlesPositions()
{
	const unsigned short nextBufferNum = (bufferNum + 1) % 2;
	const long neurotransmittersOffset = particlesOffsets[particle::NEUROTRANSMITTER];
	const long activeNeurotransmittersCountInLoopIndex = particlesOffsets[particle::NEUROTRANSMITTER] + activeNeurotransmittersCount;
	const long particlesOffset = config.particlesCount[particle::NAP] - activeNapsCount;
	const long particlesBufferSizeInLoop = particlesOffset + activeNapsCount + activeKpsCount;

	// Nap and Kp ions
#pragma loop(hint_parallel(0))
#pragma loop(ivdep)
	for (long i = particlesOffset; i < particlesBufferSizeInLoop; ++i) {
		Particle& currParticle = particles[bufferNum][i];
		Particle& prevParticle = particles[nextBufferNum][i];

		// update positions, all particles have 3 coords, x, y, z, that's why 'index * 3' in particlesPos[]
		particlesPos[i * 3 + 0] = currParticle.x;
		particlesPos[i * 3 + 1] = currParticle.y;
		particlesPos[i * 3 + 2] = currParticle.z;
	}

	// neurotransmitters
#pragma loop(hint_parallel(0))
#pragma loop(ivdep)
	for (long i = neurotransmittersOffset; i < activeNeurotransmittersCountInLoopIndex; ++i) {
		Particle& currParticle = particles[bufferNum][i];
		Particle& prevParticle = particles[nextBufferNum][i];

		// update positions, all particles have 3 coords, x, y, z, that's why 'index * 3' in particlesPos[]
		particlesPos[i * 3 + 0] = currParticle.x;
		particlesPos[i * 3 + 1] = currParticle.y;
		particlesPos[i * 3 + 2] = currParticle.z;
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

	// update particles positions buffer with their calculated positions
	updateParticlesPositions();

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
		size_t channelsOffset = 0;

		channelsRenderProgram->use();

		uniforms.channelRadius = channelRadius[channel::NAP];
		channelsRenderProgram->setUniforms(uniforms);

		glBindTexture(GL_TEXTURE_2D, NapIonChannelTexture);
		glBindVertexArray(NapIonsChannelsVAO);
		glDrawArrays(GL_POINTS, channelsOffset, NapChannelsCount);
		channelsOffset += NapChannelsCount;

		uniforms.channelRadius = channelRadius[channel::KP];
		channelsRenderProgram->setUniforms(uniforms);

		glBindTexture(GL_TEXTURE_2D, KpIonChannelTexture);
		glBindVertexArray(KpIonsChannelsVAO);
		glDrawArrays(GL_POINTS, channelsOffset, KpChannelsCount);
		channelsOffset += KpChannelsCount;
	}
	
	if (renderParticles) {
		// render particles
		size_t ionsOffset = config.particlesCount[particle::NAP] - activeNapsCount;

		ionsRenderProgram->use();

		uniforms.ionRadius = particleRadius[particle::NAP];
		ionsRenderProgram->setUniforms(uniforms);
		glBindTexture(GL_TEXTURE_2D, NapIonTexture);
		glBindVertexArray(NapIonsVAO);
		glDrawArrays(GL_POINTS, ionsOffset, activeNapsCount);

		ionsOffset = particlesOffsets[particle::KP];
		uniforms.ionRadius = particleRadius[particle::KP];
		ionsRenderProgram->setUniforms(uniforms);
		glBindTexture(GL_TEXTURE_2D, KpIonTexture);
		glBindVertexArray(KpIonsVAO);
		glDrawArrays(GL_POINTS, ionsOffset, activeKpsCount);

		ionsOffset = particlesOffsets[particle::KP];
		uniforms.ionRadius = particleRadius[particle::CLM];
		ionsRenderProgram->setUniforms(uniforms);
		glBindTexture(GL_TEXTURE_2D, ClmIonTexture);
		glBindVertexArray(ClmIonsVAO);
		glDrawArrays(GL_POINTS, ionsOffset, config.particlesCount[particle::CLM]);

		ionsOffset = particlesOffsets[particle::CLM];
		uniforms.ionRadius = particleRadius[particle::ORGANIC_ANION];
		ionsRenderProgram->setUniforms(uniforms);
		glBindTexture(GL_TEXTURE_2D, otherParticlesTexture);
		glBindVertexArray(otherParticlesVAO);
		glDrawArrays(GL_POINTS, ionsOffset, config.particlesCount[particle::ORGANIC_ANION]);

		ionsOffset = particlesOffsets[particle::ORGANIC_ANION];
		uniforms.ionRadius = particleRadius[particle::NEUROTRANSMITTER];
		ionsRenderProgram->setUniforms(uniforms);
		glBindTexture(GL_TEXTURE_2D, neurotransmittersTexture);
		glBindVertexArray(neurotransmittersVAO);
		glDrawArrays(GL_POINTS, ionsOffset, activeNeurotransmittersCount);
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
	loadConfig(config);
	setupOpenGL();
	setupInput();
	setupPrograms();
	setupStructures();
	setupUniforms();
	setupTextures();
	setupBuffers();
}

Simulation::~Simulation()
{
	glfwTerminate();
	delete neuron;
	delete ionsRenderProgram;
	delete channelsRenderProgram;
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
	activeParticlesCount = 0;
	bufferNum = 0;
	activeNeurotransmittersCount = 0;
	activeNapsCount = 0;
	const long neurotransmittersOffset = particlesBufferSize - config.particlesCount[particle::NEUROTRANSMITTER];
	for (long i = 0; i < neuron->channels.size(); ++i) {
		Channel& channel = neuron->channels[i];
		if (channel.gating == channel::CONST_OPEN)
			continue;
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
	if (n > activeParticlesCount)
		return;

	if (n > activeNeurotransmittersCount) {
		activeParticlesCount -= activeNeurotransmittersCount;
		activeNeurotransmittersCount = 0;
	}
	else {
		activeParticlesCount -= n;
		activeNeurotransmittersCount -= n;
	}
}

void Simulation::increaseNeurotransmitters(const unsigned n)
{
	const long neurotransmittersOffset = particlesBufferSize - config.particlesCount[particle::NEUROTRANSMITTER] + activeNeurotransmittersCount;

	if (activeNeurotransmittersCount < config.particlesCount[particle::NEUROTRANSMITTER] - n) {
		activeNeurotransmittersCount += n;
		activeParticlesCount += n;

		for (long i = 0; i < n; ++i) {
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
}

double Simulation::getDeltaTime(void) const
{
	return inversedTimeFactor * fabs(deltaTime);
}
