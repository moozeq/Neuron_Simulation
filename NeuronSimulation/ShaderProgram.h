#pragma once
#include "Utilities.h"

namespace shader {
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
		GLfloat channelWidth;
		GLfloat opacity;
		glm::mat4 viewMatrix;
	};
}

class ShaderProgram
{
	GLuint id;

	// uniforms
	GLuint ionRadiusLoc;
	GLuint channelRadiusLoc;
	GLuint channelWidthLoc;
	GLuint opacityLoc;

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
	ShaderProgram(const shader::ShadersTypes types, const std::vector<const GLchar*>& shadersPaths) {
		id = glCreateProgram();
		
		switch (types) {
		case shader::VF:
			if (shadersPaths.size() != 2)
				throw std::exception("Invalid shadersPaths (VF)");
			addShader(shadersPaths[0], GL_VERTEX_SHADER);
			addShader(shadersPaths[1], GL_FRAGMENT_SHADER);
			break;

		case shader::VFG:
			if (shadersPaths.size() != 3)
				throw std::exception("Invalid shadersPaths (VFG)");
			addShader(shadersPaths[0], GL_VERTEX_SHADER);
			addShader(shadersPaths[1], GL_GEOMETRY_SHADER);
			addShader(shadersPaths[2], GL_FRAGMENT_SHADER);
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
		channelWidthLoc = glGetUniformLocation(id, "channelWidth");
		opacityLoc = glGetUniformLocation(id, "opacity");

		// view uniform here
		viewMatrixLoc = glGetUniformLocation(id, "viewMatrix");
	}
	~ShaderProgram() {
		glDeleteProgram(id);
	}

	// setting uniforms
	void setUniforms(shader::Uniforms& uniforms) {
		glUniform1f(ionRadiusLoc, uniforms.ionRadius);
		glUniform1f(channelRadiusLoc, uniforms.channelRadius);
		glUniform1f(channelWidthLoc, uniforms.channelWidth);
		glUniform1f(opacityLoc, uniforms.opacity);
		glUniformMatrix4fv(viewMatrixLoc, 1, GL_FALSE, glm::value_ptr(uniforms.viewMatrix));
	}

	void use(void) const {
		glUseProgram(id);
	}
};
