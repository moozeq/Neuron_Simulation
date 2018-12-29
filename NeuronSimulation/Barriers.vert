#version 450 core

layout (location = 0) in vec3 position;
layout (location = 1) in vec2 texCoord;
layout (location = 2) in vec3 normal;

out vec3 vertFragPos;
out vec2 vertTexCoord;
out vec3 vertNormal;

uniform mat4 viewMatrix;

void main()
{
    gl_Position = viewMatrix * vec4(position, 1.0f);

	vertFragPos = position;
	vertTexCoord = texCoord;
	vertNormal = normal;
} 