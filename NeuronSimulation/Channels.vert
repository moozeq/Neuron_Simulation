#version 450 core

layout (location = 0) in vec4 channel;

out float vertChannelState;

void main(void) {
	gl_Position = vec4(channel.xyz, 1.0f);
	vertChannelState = channel[3];
}