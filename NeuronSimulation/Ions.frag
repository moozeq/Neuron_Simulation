#version 450 core

uniform sampler2D ionTexture;

in vec2 texCoord;

out vec4 color;

void main(void) {
	color = texture(ionTexture, texCoord);
}