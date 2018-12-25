#version 450 core

uniform sampler2D channelTexture;

in vec2 texCoord;
//in float geomChannelState;

out vec4 color;

void main(void) {
	color = texture(channelTexture, texCoord);
	//color.a = geomChannelState;
}