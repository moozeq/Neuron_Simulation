#include "Simulation.h"

void Simulation::loadConfig(const Config& _config)
{
	config = _config;
	metricFactor = config.metricFactor;
	timeFactor = config.timeFactor;

	inversedTimeFactor = 1.0 / config.timeFactor;
	particlesBufferSize = config.KpIonsNum + config.NapIonsNum + config.ClmIonsNum;
	bufferNum = 0;
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

	//camera = Camera(glm::vec3(0.0f, 0.0f, 10.0f));
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
	while (i++ < offset) {
		Particle particle = Particle({ GET_RAND_DOUBLE(-0.10, 0.10), GET_RAND_DOUBLE(-0.10, 0.10), GET_RAND_DOUBLE(-0.10, 0.10), 0.0, 0.0, 0.0, phy::NapC, phy::NapM, i });
		particles[0].push_back(particle);
		particles[1].push_back(particle);
	}

	offset += config.KpIonsNum;
	while (i++ < offset) {
		Particle particle = Particle({ GET_RAND_DOUBLE(-0.10, 0.10), GET_RAND_DOUBLE(-0.10, 0.10), GET_RAND_DOUBLE(-0.10, 0.10), 0.0, 0.0, 0.0, phy::KpC, phy::KpM, i });
		particles[0].push_back(particle);
		particles[1].push_back(particle);
	}

	offset += config.ClmIonsNum;
	while (i++ < offset) {
		Particle particle = Particle({ GET_RAND_DOUBLE(-0.10, 0.10), GET_RAND_DOUBLE(-0.10, 0.10), GET_RAND_DOUBLE(-0.10, 0.10), 0.0, 0.0, 0.0, phy::ClmC, phy::ClmM, i });
		particles[0].push_back(particle);
		particles[1].push_back(particle);
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
	sh::Uniforms uniforms;
	uniforms.ionRadius = config.ionRadius;
	ionsRenderProgram->setUniforms(uniforms);
	glUseProgram(0);

	// create vaos
	glCreateVertexArrays(1, &NapIonsVAO);
	glCreateVertexArrays(1, &KpIonsVAO);
	glCreateVertexArrays(1, &ClmIonsVAO);

	// creating buffers
	glCreateBuffers(1, &particlesPosBuf);

	GLbitfield flags = GL_MAP_WRITE_BIT | GL_MAP_PERSISTENT_BIT | GL_MAP_COHERENT_BIT;
	GLsizei bufferSize = particlesBufferSize * 3 * sizeof(GLfloat);
	size_t offset = 0;
	glNamedBufferStorage(particlesPosBuf, bufferSize, nullptr, flags);

	glVertexArrayVertexBuffer(NapIonsVAO, 0, particlesPosBuf, offset, 3 * sizeof(GL_FLOAT));
	glVertexArrayAttribFormat(NapIonsVAO, 0, 3, GL_FLOAT, GL_FALSE, 0);
	glVertexArrayAttribBinding(NapIonsVAO, 0, 0);
	glEnableVertexArrayAttrib(NapIonsVAO, 0);
	offset += config.NapIonsNum;

	glVertexArrayVertexBuffer(KpIonsVAO, 0, particlesPosBuf, offset, 3 * sizeof(GL_FLOAT));
	glVertexArrayAttribFormat(KpIonsVAO, 0, 3, GL_FLOAT, GL_FALSE, 0);
	glVertexArrayAttribBinding(KpIonsVAO, 0, 0);
	glEnableVertexArrayAttrib(KpIonsVAO, 0);
	offset += config.KpIonsNum;

	glVertexArrayVertexBuffer(ClmIonsVAO, 0, particlesPosBuf, offset, 3 * sizeof(GL_FLOAT));
	glVertexArrayAttribFormat(ClmIonsVAO, 0, 3, GL_FLOAT, GL_FALSE, 0);
	glVertexArrayAttribBinding(ClmIonsVAO, 0, 0);
	glEnableVertexArrayAttrib(ClmIonsVAO, 0);
	offset += config.ClmIonsNum;

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
		//mouseButtonCallback(window, GLFW_MOUSE_BUTTON_LEFT, GLFW_PRESS, 0);
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
	const long ionsBufferSize = particlesBufferSize;
	
#pragma loop(hint_parallel(0))
#pragma loop(ivdep)
	for (long i = 0; i < ionsBufferSize; ++i) {
		Particle& currParticle = particles[bufferNum][i];
		const size_t index = currParticle.index;
		accels[index] = 0.0;
		accels[index + 1] = 0.0;
		accels[index + 2] = 0.0;

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

			accels[index] -= a * dx / d;
			accels[index + 1] -= a * dy / d;
			accels[index + 2] -= a * dz / d;
		}
		currParticle.x += currParticle.vx * deltaTime + accels[index] * deltaTime * deltaTime / 2;
		currParticle.y += currParticle.vy * deltaTime + accels[index + 1] * deltaTime * deltaTime / 2;
		currParticle.z += currParticle.vz * deltaTime + accels[index + 2] * deltaTime * deltaTime / 2;

		currParticle.vx += accels[index] * deltaTime;
		currParticle.vy += accels[index + 1] * deltaTime;
		currParticle.vz += accels[index + 2] * deltaTime;
	}
}

inline void Simulation::update()
{
	updateIons();
	const long ionsBufferSize = particlesBufferSize;

	for (long i = 0; i < ionsBufferSize; ++i) {
		Particle& currParticle = particles[bufferNum][i];
		particlesPos[i] = currParticle.x;
		particlesPos[i + 1] = currParticle.y;
		particlesPos[i + 2] = currParticle.z;
	}

	++bufferNum %= 2;
}

void Simulation::render()
{
	glClearColor(0.3f, 0.3f, 0.3f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	ionsRenderProgram->use();
	glm::mat4 viewMatrix = glm::perspective(glm::radians(camera.zoom), (float)width / (float)height, 0.1f, 100.0f) * camera.GetViewMatrix();
	ionsRenderProgram->setViewMatrix(viewMatrix);

	glBindTexture(GL_TEXTURE_2D, NapIonTexture);
	glBindVertexArray(NapIonsVAO);
	glDrawArrays(GL_POINTS, 0, config.NapIonsNum);

	glBindTexture(GL_TEXTURE_2D, KpIonTexture);
	glBindVertexArray(KpIonsVAO);
	glDrawArrays(GL_POINTS, 0, config.KpIonsNum);

	glBindTexture(GL_TEXTURE_2D, ClmIonTexture);
	glBindVertexArray(ClmIonsVAO);
	glDrawArrays(GL_POINTS, 0, config.ClmIonsNum);

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

double Simulation::getDeltaTime(void) const
{
	return inversedTimeFactor * deltaTime;
}
