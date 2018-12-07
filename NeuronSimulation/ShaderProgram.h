#pragma once
#include "Utilities.h"

namespace sh {
	/*
	VF = vert + frag
	VFG = vert + frag + geo
	*/
	enum ShadersTypes {
		VF, VFG
	};

	enum Shaders {
		VERT, GEO, FRAG
	};
	
/*
to add new uniform, change:
	Uniforms name
	ShaderProgram location (& getloc)
	setUniforms glUniform*(location, name)
*/

	struct Uniforms {
		GLfloat ionRadius;
		GLfloat channelRadius;
		GLint NapIonTexture;
		GLint NApChannelTexture;
		GLint KpIonTexture;
		GLint KpChannelTexture;
	};
}

class ShaderProgram
{
	GLuint id;

	// uniforms
	GLuint ionRadiusLoc;
	GLuint channelRadiusLoc;

	GLuint NapIonTextureLoc;
	GLuint NApChannelTextureLoc;

	GLuint KpIonTextureLoc;
	GLuint KpChannelTextureLoc;

	// view uniforms
	GLuint viewMatrixLoc;

	std::string readShaderCode(const GLchar* shaderPath)
	{
		std::ifstream shaderFile;
		shaderFile.exceptions(std::ifstream::badbit);
		shaderFile.open(shaderPath);

		// check for errors while opening file and throw if any
		if (!shaderFile.is_open())
			throw std::exception("Shader source code not found");
		std::stringstream shaderStream;
		shaderStream << shaderFile.rdbuf();
		shaderFile.close();
		return shaderStream.str();
	}

	GLuint compileShader(const GLchar* shaderCode, GLenum shaderType)
	{
		GLuint shaderId = glCreateShader(shaderType);
		glShaderSource(shaderId, 1, &shaderCode, NULL);
		glCompileShader(shaderId);

		// check for errors and throw if any
		GLint success = 0;
		glGetShaderiv(shaderId, GL_COMPILE_STATUS, &success);
		if (!success)
		{
			GLchar infoLog[512];
			glGetShaderInfoLog(shaderId, sizeof(infoLog), NULL, infoLog);
			std::string msg = std::string("Shader compilation error: ") + infoLog;
			throw std::exception(msg.c_str());
		}
		return shaderId;
	}

	void addShader(const GLchar* shaderPath, GLenum shaderType) {
		std::string shaderCode = readShaderCode(shaderPath);
		GLuint shaderId = compileShader(shaderCode.c_str(), shaderType);
		glAttachShader(id, shaderId);

		// delete shaders as they're linked into our program now and no longer necessery
		glDeleteShader(shaderId);
	}

public:
	ShaderProgram(const sh::ShadersTypes types, const std::vector<const GLchar*>& shadersPaths) {
		id = glCreateProgram();

		switch (types) {
		case sh::VF:
			addShader(shadersPaths[sh::VERT], GL_VERTEX_SHADER);
			addShader(shadersPaths[sh::FRAG], GL_FRAGMENT_SHADER);
			break;

		case sh::VFG:
			addShader(shadersPaths[sh::VERT], GL_VERTEX_SHADER);
			addShader(shadersPaths[sh::GEO], GL_GEOMETRY_SHADER);
			addShader(shadersPaths[sh::FRAG], GL_FRAGMENT_SHADER);
			break;
		}

		// link shader program
		glLinkProgram(id);

		// check for errors in linking and throw if any
		GLint success;
		glGetProgramiv(id, GL_LINK_STATUS, &success);
		if (!success)
		{
			GLchar infoLog[512];
			glGetProgramInfoLog(id, sizeof(infoLog), NULL, infoLog);
			std::string msg = std::string("Shader program linking error: ") + infoLog;
			throw std::exception(msg.c_str());
		}

		// uniforms here
		ionRadiusLoc = glGetUniformLocation(id, "ionRadius");
		channelRadiusLoc = glGetUniformLocation(id, "channelRadius");
		NapIonTextureLoc = glGetUniformLocation(id, "NapIonTexture");
		NApChannelTextureLoc = glGetUniformLocation(id, "NApChannelTexture");
		KpIonTextureLoc = glGetUniformLocation(id, "KpIonTexture");
		KpChannelTextureLoc = glGetUniformLocation(id, "KpChannelTexture");

		// view uniform here
		viewMatrixLoc = glGetUniformLocation(id, "viewMatrix");
	}
	~ShaderProgram() {
		glDeleteProgram(id);
	}

	// setting uniforms
	void setUniforms(sh::Uniforms& uniforms) {
		glUniform1f(ionRadiusLoc, uniforms.ionRadius);
		glUniform1f(channelRadiusLoc, uniforms.channelRadius);
		glUniform1i(NapIonTextureLoc, uniforms.NapIonTexture);
		glUniform1i(NApChannelTextureLoc, uniforms.NApChannelTexture);
		glUniform1i(KpIonTextureLoc, uniforms.KpIonTexture);
		glUniform1i(KpChannelTextureLoc, uniforms.KpChannelTexture);
	}

	void setViewMatrix(glm::mat4 viewMatrix) {
		glUniformMatrix4fv(viewMatrixLoc, 1, GL_FALSE, glm::value_ptr(viewMatrix));
	}

	void use(void) const {
		glUseProgram(id);
	}
};
