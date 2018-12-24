#version 450 core

out vec4 color;

uniform float opacity;
uniform uint red;
uniform uint green;
uniform uint blue;

void main()
{
	color = vec4(red, green, blue, opacity);
} 