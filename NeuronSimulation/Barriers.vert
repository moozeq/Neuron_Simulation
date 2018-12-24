#version 450 core

layout (location = 0) in vec3 position;

uniform mat4 viewMatrix;

void main()
{
    gl_Position = viewMatrix * vec4(position, 1.0f);
} 