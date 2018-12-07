#version 450 core

layout (location = 0) in vec3 ion;

void main(void) {
	gl_Position = vec4(ion.xyz, 1.0f);
}