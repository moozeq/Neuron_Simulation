#version 450 core

layout (location = 0) in vec3 ions;

void main(void) {
	gl_Position = vec4(ions.xyz, 1.0f);
}