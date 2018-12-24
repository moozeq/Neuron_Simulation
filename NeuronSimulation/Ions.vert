#version 450 core

layout (location = 0) in vec3 ionPosition;

void main(void) {
	gl_Position = vec4(ionPosition.xyz, 1.0f);
}