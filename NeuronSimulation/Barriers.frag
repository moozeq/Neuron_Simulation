#version 450 core

uniform sampler2D bilayerTexture;
uniform float opacity;

in vec3 vertFragPos;
in vec2 vertTexCoord;
in vec3 vertNormal;

out vec4 color;

void main()
{
	color = texture(bilayerTexture, vertTexCoord);
	color.a = opacity;
} 