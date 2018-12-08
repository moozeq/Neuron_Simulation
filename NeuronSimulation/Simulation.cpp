#include "Simulation.h"

void Simulation::loadConfig(const Config& _config)
{
	config = _config;
	metricFactor = config.metricFactor;
	timeFactor = config.timeFactor;

	inversedTimeFactor = 1.0 / config.timeFactor;
	particlesBufferSize = config.KpIonsNum + config.NapIonsNum + config.ClmIonsNum + config.otherParticlesNum;
	bufferNum = 0;
	ice = false;
}

void Simulation::setupOpenGL()
{
	// log info about client
	logfile.open(config.logPath, std::ofstream::app | std::ofstream::binary);
	if (!logfile.good())
		throw("Couldn't open log file");

	// setup client's openGL components
	log(logfile, "--- Simulation initialization's started ---");

	if (glfwInit() != GL_TRUE)
		throw("GLFW initialization's failed");

	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 5);

	width = config.width;
	height = config.height;

	window = glfwCreateWindow(width, height, "ProjectN", nullptr, nullptr);
	if (window == nullptr)
		throw("GLFW window couldn't be created");
	glfwMakeContextCurrent(window);

	// new glew functions usage
	glewExperimental = GL_TRUE;
	if (glewInit() != GLEW_OK)
		throw("GLEW initialization's failed");

	// enable alpha channel in textures
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_BLEND);
	glEnable(GL_DEPTH_TEST);

	int integ;
	glGetIntegerv(GL_MAX_GEOMETRY_OUTPUT_VERTICES, &integ);

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
	ionsRenderProgram = new ShaderProgram(sh::VFG, ionsPaths);
}

void Simulation::setupStructures()
{
	particles[0].reserve(particlesBufferSize);
	particles[1].reserve(particlesBufferSize);
	size_t i = 0;
	size_t offset = 0;

	offset += config.NapIonsNum;
	while (i < offset) {
		Particle* particle = newParticle(0.1, par::NAP, i++);
		particles[0].push_back(*particle);
		particles[1].push_back(*particle);
		delete particle;
	}

	offset += config.KpIonsNum;
	while (i < offset) {
		Particle* particle = newParticle(0.1, par::KP, i++);
		particles[0].push_back(*particle);
		particles[1].push_back(*particle);
		delete particle;
	}

	offset += config.ClmIonsNum;
	while (i < offset) {
		Particle* particle = newParticle(0.1, par::CLM, i++);
		particles[0].push_back(*particle);
		particles[1].push_back(*particle);
		delete particle;
	}

	offset += config.otherParticlesNum;
	while (i < offset) {
		Particle* particle = newParticle(0.1, par::MASSIVEION, i++);
		particles[0].push_back(*particle);
		particles[1].push_back(*particle);
		delete particle;
	}

	accels.reserve(particlesBufferSize * 3);
	accels.resize(particlesBufferSize * 3);
}

void Simulation::setupBuffers()
{
	// set const values as uniforms in shader program
	ionsRenderProgram->use();
	NapIonTexture = loadMipmapTexture(GL_TEXTURE0, config.NapIonTexturePath.c_str());
	KpIonTexture = loadMipmapTexture(GL_TEXTURE0, config.KpIonTexturePath.c_str());
	ClmIonTexture = loadMipmapTexture(GL_TEXTURE0, config.ClmIonTexturePath.c_str());
	otherParticlesTexture = loadMipmapTexture(GL_TEXTURE0, config.otherParticlesTexturePath.c_str());
	sh::Uniforms uniforms;
	uniforms.ionRadius = config.ionRadius;
	ionsRenderProgram->setUniforms(uniforms);
	glUseProgram(0);

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
		throw("Buffer mapping failed");

	glUnmapNamedBuffer(particlesPosBuf);
}

void Simulation::keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	// exit
	if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
		glfwSetWindowShouldClose(window, GL_TRUE);

	Simulation* simulation = static_cast<Simulation*>(glfwGetWindowUserPointer(window));

	if (key == GLFW_KEY_SPACE && action == GLFW_PRESS)
		simulation->freeze();

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

inline void Simulation::updateIons()
{
	const unsigned short nextBufferNum = (++bufferNum) % 2;
	const long ionsBufferSize = (long)particlesBufferSize;
	
#pragma loop(hint_parallel(0))
#pragma loop(ivdep)
	for (long i = 0; i < ionsBufferSize; ++i) {
		Particle& currParticle = particles[bufferNum][i];

		// all particles have 3 accels, ax, ay, az, that's why index * 3
		const size_t accelIndex = currParticle.index * 3;

		for (long j = 0; j < ionsBufferSize; ++j) {
			const Particle& particle = particles[bufferNum][j];
			double dx = metricFactor * (particle.x - currParticle.x);
			double dy = metricFactor * (particle.y - currParticle.y);
			double dz = metricFactor * (particle.z - currParticle.z);

			double d = cbrt(dx * dx + dy * dy + dz * dz);
			if (d == 0.0)
				continue;
			double F = phy::k * currParticle.charge * particle.charge / d;
			double a = F / currParticle.mass;

			accels[accelIndex] -= a * dx / d;
			accels[accelIndex + 1] -= a * dy / d;
			accels[accelIndex + 2] -= a * dz / d;
		}
		currParticle.x += currParticle.vx * deltaTime + accels[accelIndex] * deltaTime * deltaTime / 2;
		currParticle.y += currParticle.vy * deltaTime + accels[accelIndex + 1] * deltaTime * deltaTime / 2;
		currParticle.z += currParticle.vz * deltaTime + accels[accelIndex + 2] * deltaTime * deltaTime / 2;

		currParticle.vx += accels[accelIndex] * deltaTime;
		currParticle.vy += accels[accelIndex + 1] * deltaTime;
		currParticle.vz += accels[accelIndex + 2] * deltaTime;

		particlesPos[accelIndex] = (float)currParticle.x;
		particlesPos[accelIndex + 1] = (float)currParticle.y;
		particlesPos[accelIndex + 2] = (float)currParticle.z;
	}
}

inline void Simulation::update()
{
	if (ice)
		return;

	// erase accels vector
	std::fill(accels.begin(), accels.end(), 0);
	updateIons();

	++bufferNum %= 2;
}

void Simulation::render()
{
	glClearColor(0.3f, 0.3f, 0.3f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	ionsRenderProgram->use();
	glm::mat4 viewMatrix = glm::perspective(glm::radians(camera.zoom), (float)width / (float)height, 0.1f, 100.0f) * camera.GetViewMatrix();
	ionsRenderProgram->setViewMatrix(viewMatrix);

	size_t offset = 0;

	glBindTexture(GL_TEXTURE_2D, NapIonTexture);
	glBindVertexArray(NapIonsVAO);
	glDrawArrays(GL_POINTS, offset, config.NapIonsNum);
	offset += config.NapIonsNum;

	glBindTexture(GL_TEXTURE_2D, KpIonTexture);
	glBindVertexArray(KpIonsVAO);
	glDrawArrays(GL_POINTS, offset, config.KpIonsNum);
	offset += config.KpIonsNum;

	glBindTexture(GL_TEXTURE_2D, ClmIonTexture);
	glBindVertexArray(ClmIonsVAO);
	glDrawArrays(GL_POINTS, offset, config.ClmIonsNum);
	offset += config.ClmIonsNum;

	glBindTexture(GL_TEXTURE_2D, otherParticlesTexture);
	glBindVertexArray(otherParticlesVAO);
	glDrawArrays(GL_POINTS, offset, config.otherParticlesNum);
	offset += config.otherParticlesNum;

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
	setupBuffers();
}

Simulation::~Simulation()
{
	glfwDestroyWindow(window);
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

double Simulation::getDeltaTime(void) const
{
	return inversedTimeFactor * deltaTime;
}
