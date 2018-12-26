#include "Simulation.h"

void Simulation::loadConfig(const Config& _config)
{
	// copy config struct
	config = _config;

	// set variables based on config struct
	metricFactorSq = config.metricFactor * config.metricFactor;
	metricFactor = config.metricFactor;
	timeFactor = config.timeFactor;

	inversedTimeFactor = 1.0 / config.timeFactor;
	particlesBufferSize = config.NapIonsNum + config.KpIonsNum + config.ClmIonsNum + config.otherParticlesNum;
	channelsBufferSize = config.NapIonsChannelsNum + config.KpIonsChannelsNum;
	bufferNum = 0;
	ice = false;
	rewind = false;
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
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_BLEND);
	glEnable(GL_DEPTH_TEST);

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
	neuron = new Neuron();
}

void Simulation::setupParticlesStructures()
{
	particles[0].reserve(particlesBufferSize);
	particles[1].reserve(particlesBufferSize);
	partAccOrigin.reserve(particlesBufferSize);
	size_t i = 0;
	size_t offset = 0;

	offset += config.NapIonsNum;
	while (i < offset) {
		Particle* particle = newParticle(0.1, particle::NAP);
		particles[0].push_back(*particle);
		particles[1].push_back(*particle);
		partAccOrigin.push_back(phy::NapA / metricFactorSq);
		delete particle;
		++i;
	}

	offset += config.KpIonsNum;
	while (i < offset) {
		Particle* particle = newParticle(0.1, particle::KP);
		particles[0].push_back(*particle);
		particles[1].push_back(*particle);
		partAccOrigin.push_back(phy::KpA / metricFactorSq);
		delete particle;
		++i;
	}

	offset += config.ClmIonsNum;
	while (i < offset) {
		Particle* particle = newParticle(0.1, particle::CLM);
		particles[0].push_back(*particle);
		particles[1].push_back(*particle);
		partAccOrigin.push_back(phy::ClmA / metricFactorSq);
		delete particle;
		++i;
	}

	offset += config.otherParticlesNum;
	while (i < offset) {
		Particle* particle = newParticle(0.1, particle::MASSIVEION);
		particles[0].push_back(*particle);
		particles[1].push_back(*particle);
		partAccOrigin.push_back((phy::k * particle->charge / particle->mass) / metricFactorSq);
		delete particle;
		++i;
	}
}

void Simulation::setupTextures()
{
	// set const values as uniforms in shader program
	shader::Uniforms uniforms;
	uniforms.ionRadius = config.ionRadius;
	uniforms.channelRadius = config.channelRadius;

	// ions
	ionsRenderProgram->use();

	// load textures
	NapIonTexture = loadMipmapTexture(GL_TEXTURE0, config.NapIonTexturePath.c_str());
	KpIonTexture = loadMipmapTexture(GL_TEXTURE0, config.KpIonTexturePath.c_str());
	ClmIonTexture = loadMipmapTexture(GL_TEXTURE0, config.ClmIonTexturePath.c_str());
	otherParticlesTexture = loadMipmapTexture(GL_TEXTURE0, config.otherParticlesTexturePath.c_str());

	// set uniforms
	ionsRenderProgram->setUniforms(uniforms);
	glUseProgram(0);

	// channels
	channelsRenderProgram->use();

	// load textures
	for (int i = 0; i < channel::STATES_COUNT; ++i)
		NapIonChannelTexture[i] = loadMipmapTexture(GL_TEXTURE0, config.NapChannelTexturePath[i].c_str());
	for (int i = 0; i < channel::STATES_COUNT; ++i)
		KpIonChannelTexture[i] = loadMipmapTexture(GL_TEXTURE0, config.KpChannelTexturePath[i].c_str());
	
	// set uniforms
	channelsRenderProgram->setUniforms(uniforms);
	glUseProgram(0);
}

void Simulation::setupBuffers()
{
	setupParticlesBuffers();
	setupChannelsBuffers();
}

void Simulation::setupParticlesBuffers()
{
	GLuint particlesPosBuf;

	// create vaos
	glCreateVertexArrays(1, &NapIonsVAO);
	glCreateVertexArrays(1, &KpIonsVAO);
	glCreateVertexArrays(1, &ClmIonsVAO);
	glCreateVertexArrays(1, &otherParticlesVAO);

	// creating buffers
	glCreateBuffers(1, &particlesPosBuf);

	GLbitfield flags = GL_MAP_WRITE_BIT | GL_MAP_PERSISTENT_BIT | GL_MAP_COHERENT_BIT;
	GLsizei bufferSize = (GLsizei)particlesBufferSize * 3 * sizeof(GLfloat);
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

	// initialize buffers in GPU and get pointers to them
	particlesPos = (float*)glMapNamedBuffer(particlesPosBuf, GL_WRITE_ONLY);
	if (!particlesPos)
		throw std::exception("Buffer mapping failed");

	glUnmapNamedBuffer(particlesPosBuf);
}

void Simulation::setupChannelsBuffers()
{
	GLuint channelsBuf;

	// create vaos
	glCreateVertexArrays(1, &NapIonsChannelsVAO);
	glCreateVertexArrays(1, &KpIonsChannelsVAO);

	// creating buffers
	glCreateBuffers(1, &channelsBuf);

	GLbitfield flags = GL_MAP_WRITE_BIT | GL_MAP_PERSISTENT_BIT | GL_MAP_COHERENT_BIT;
	GLsizei channelBufferSize = (GLsizei)neuron->channels.size() * 4 * sizeof(GLfloat);
	std::vector<float> channels = neuron->getChannels();

	glNamedBufferStorage(channelsBuf, channelBufferSize, &channels[0], flags);

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

	// TODO lock current simulation state before changing timeFactor
	// RIGHT - time speed + 20%
	if (key == GLFW_KEY_RIGHT && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
		simulation->timeFactor *= 1.2f;
		std::ostringstream timeFactor;
		timeFactor << simulation->timeFactor;
		log(simulation->logfile, "[+] Time factor = " + timeFactor.str());

		// after changing time factor needs to change movement speed
		simulation->camera.movementSpeed *= 0.8f;
		log(simulation->logfile, "[-] Camera movement speed = " + std::to_string(simulation->camera.movementSpeed));
	}

	// LEFT - time speed - 20%
	if (key == GLFW_KEY_LEFT && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
		simulation->timeFactor *= 0.8f;
		std::ostringstream timeFactor;
		timeFactor << simulation->timeFactor;
		log(simulation->logfile, "[-] Time factor = " + timeFactor.str());

		// after changing time factor needs to change movement speed
		simulation->camera.movementSpeed *= 1.2f;
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

	// parallelization can be done due to no influence from other channels
#pragma loop(hint_parallel(0))
#pragma loop(ivdep)
	for (long i = 0; i < channelsBufferSizeInLoop; ++i) {
		Channel& currChannel = neuron->channels[i];

		Particle* particle;
		float dx, dy, dz, d;
		float Ein = 0.0f;
		float Eout = 0.0f;
		float U = 0.0f;

		// calc electric field inside neuron
		for (long j = 0; j < particlesBufferSize; ++j) {
			particle = &particles[bufferNum][j];

			dx = particle->x - currChannel.xIn;
			dy = particle->y - currChannel.yIn;
			dz = particle->z - currChannel.zIn;

			d = sqrt(dx * dx + dy * dy + dz * dz);
			if (d == 0.0)
				continue;

			Ein += phy::k * particle->charge / (metricFactor * d);
		}

		// calc electric field outside neuron
		for (long j = 0; j < particlesBufferSize; ++j) {
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

		// TODO probability to open instead of threshold
		// TODO add relative refraction 
		// TODO refactor, delete logs

		std::ostringstream streamObj;
		switch (currChannel.state) {
		case channel::CLOSED:
			if (U > phy::NapOpenTreshold) {

				log(logfile, "[C] Channel opened, U = " + std::to_string(U));

				currChannel.state = channel::OPEN;
				currChannel.timeLeft = config.timeFactor / phy::NapOpenTime;
				channelsAttribs[i * 4 + 3] = 1.0f;
			}
			break;
		case channel::OPEN:
			
			currChannel.timeLeft -= deltaTime;

			log(logfile, "[O] Channel open, U = " + std::to_string(U));
			streamObj << currChannel.timeLeft;
			log(logfile, "[O] Channel open, time left = " + streamObj.str());

			if (currChannel.timeLeft < 0.0f) {

				log(logfile, "[I] Channel inactivated, U = " + std::to_string(U));

				currChannel.state = channel::INACTIVE;
				currChannel.timeLeft = config.timeFactor / phy::NapInactiveTime;
				channelsAttribs[i * 4 + 3] = 0.5f;
			}
			break;
		case channel::INACTIVE:
			currChannel.timeLeft -= deltaTime;

			log(logfile, "[I] Channel inactive, U = " + std::to_string(U));
			streamObj << currChannel.timeLeft;
			log(logfile, "[I] Channel inactive, time left = " + streamObj.str());

			if (currChannel.timeLeft < 0.0f) {

				log(logfile, "[C] Channel closed, U = " + std::to_string(U));

				currChannel.state = channel::CLOSED;
				channelsAttribs[i * 4 + 3] = 0.0f;
			}
			break;
		}
	}
}

inline void Simulation::updateParticlesPositions()
{
	const unsigned short nextBufferNum = (bufferNum + 1) % 2;
	const long particlesBufferSizeInLoop = particlesBufferSize;
	
	// parallelization can be done due to double buffering particles vector
#pragma loop(hint_parallel(0))
#pragma loop(ivdep)
	for (long i = 0; i < particlesBufferSizeInLoop; ++i) {
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
		for (long j = 0; j < i; ++j) {
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
		for (long j = i + 1; j < particlesBufferSize; ++j) {
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

		// old coordinates
		float ox = prevParticle.x;
		float oy = prevParticle.y;
		float oz = prevParticle.z;

		// muliply by part of acceleration which depends on current particle
		double axdt = ax * partAccOriginDt;
		double aydt = ay * partAccOriginDt;
		double azdt = az * partAccOriginDt;

		// new coordinates
		// new coordinates: x = x0 + vx0 * dt + ax * dt^2 / 2
		// 'axdt = ax * dt' so optimized equation: x += dt * (vx0 + axdt / 2)
		float nx = ox + deltaTime * (currParticle.vx + axdt / 2);
		float ny = oy + deltaTime * (currParticle.vy + aydt / 2);
		float nz = oz + deltaTime * (currParticle.vz + azdt / 2);

		// all particles have 3 coords, x, y, z, that's why 'index * 3' in particlesPos[]
		particlesPos[i * 3 + 0] = currParticle.x = nx;
		particlesPos[i * 3 + 1] = currParticle.y = ny;
		particlesPos[i * 3 + 2] = currParticle.z = nz;

		// update velocities: vx = vx0 + ax * dt
		currParticle.vx += axdt;
		currParticle.vy += aydt;
		currParticle.vz += azdt;

		// collision detection and reaction
		if (nx <= -1.0f || nx >= 1.0f) {
			if (nx <= -1.0f)
				particlesPos[i * 3 + 0] = currParticle.x = nx = -2.0f - nx;
			else
				particlesPos[i * 3 + 0] = currParticle.x = nx = 2.0f - nx;

			currParticle.vx -= 2 * currParticle.vx;
		}

		if (ny <= -1.0f || ny >= 1.0f) {
			if (ny <= -1.0f)
				particlesPos[i * 3 + 1] = currParticle.y = ny = -2.0f - ny;
			else
				particlesPos[i * 3 + 1] = currParticle.y = ny = 2.0f - ny;

			currParticle.vy -= 2 * currParticle.vy;
		}

		if (nz <= -1.0f || nz >= 1.0f) {
			if (nz <= -1.0f)
				particlesPos[i * 3 + 2] = currParticle.z = nz = -2.0f - nz;
			else
				particlesPos[i * 3 + 2] = currParticle.z = nz = 2.0f - nz;

			currParticle.vz -= 2 * currParticle.vz;
		}
	}
}

inline void Simulation::update()
{
	if (ice)
		return;

	// update channles states
	updateChannelsStates();

	// update particles in current vector and positions buffer
	updateParticlesPositions();

	// swap to next particles vector
	++bufferNum %= 2;
}

inline void Simulation::render()
{
	glClearColor(0.3f, 0.3f, 0.3f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glm::mat4 viewMatrix = glm::perspective(glm::radians(camera.zoom), (float)width / (float)height, 0.1f, 100.0f) * camera.GetViewMatrix();

	// render neuron
	neuron->render();

	// render channels
	channelsRenderProgram->use();
	channelsRenderProgram->setViewMatrix(viewMatrix);

	size_t channelsOffset = 0;

	glBindTexture(GL_TEXTURE_2D, NapIonChannelTexture[channel::NONE]);
	glBindVertexArray(NapIonsChannelsVAO);
	glDrawArrays(GL_POINTS, channelsOffset, config.NapIonsChannelsNum);
	channelsOffset += config.NapIonsChannelsNum;

	glBindTexture(GL_TEXTURE_2D, KpIonChannelTexture[channel::NONE]);
	glBindVertexArray(KpIonsChannelsVAO);
	glDrawArrays(GL_POINTS, channelsOffset, config.KpIonsChannelsNum);
	channelsOffset += config.KpIonsChannelsNum;

	// render ions
	ionsRenderProgram->use();
	ionsRenderProgram->setViewMatrix(viewMatrix);

	size_t ionsOffset = 0;

	glBindTexture(GL_TEXTURE_2D, NapIonTexture);
	glBindVertexArray(NapIonsVAO);
	glDrawArrays(GL_POINTS, ionsOffset, config.NapIonsNum);
	ionsOffset += config.NapIonsNum;

	glBindTexture(GL_TEXTURE_2D, KpIonTexture);
	glBindVertexArray(KpIonsVAO);
	glDrawArrays(GL_POINTS, ionsOffset, config.KpIonsNum);
	ionsOffset += config.KpIonsNum;

	glBindTexture(GL_TEXTURE_2D, ClmIonTexture);
	glBindVertexArray(ClmIonsVAO);
	glDrawArrays(GL_POINTS, ionsOffset, config.ClmIonsNum);
	ionsOffset += config.ClmIonsNum;

	glBindTexture(GL_TEXTURE_2D, otherParticlesTexture);
	glBindVertexArray(otherParticlesVAO);
	glDrawArrays(GL_POINTS, ionsOffset, config.otherParticlesNum);
	ionsOffset += config.otherParticlesNum;

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
	setupTextures();
	setupBuffers();
}

Simulation::~Simulation()
{
	glfwTerminate();
	delete neuron;
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

double Simulation::getDeltaTime(void) const
{
	return inversedTimeFactor * fabs(deltaTime);
}
