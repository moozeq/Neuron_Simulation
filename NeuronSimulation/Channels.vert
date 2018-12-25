#version 450 core

layout (location = 0) in vec3 channelPosition;
layout (location = 1) in float channelState;

//out float vertChannelState;

void main(void) {
	gl_Position = vec4(channelPosition.xyz, 1.0f);
	//vertChannelState = channelState;
}