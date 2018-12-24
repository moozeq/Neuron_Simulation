#version 450 core

layout (location = 0) in vec3 channelPosition;

void main(void) {
	gl_Position = vec4(channelPosition.xyz, 1.0f);
}