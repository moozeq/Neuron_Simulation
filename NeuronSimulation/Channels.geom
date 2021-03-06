#version 450 core

layout (points) in;
layout (triangle_strip, max_vertices = 24) out;

uniform float channelRadius;
uniform float channelWidth;
uniform mat4 viewMatrix;

in float[1] vertChannelState;

out vec2 texCoord;
out float geomChannelState;

void main(void) {
	const vec3 offset[] = vec3[] (
							vec3(-channelRadius, channelRadius, channelRadius), vec3(channelRadius, channelRadius, channelRadius), vec3(-channelRadius, channelRadius, -channelRadius), vec3(channelRadius, channelRadius, -channelRadius),
							vec3(-channelRadius, -channelRadius, channelRadius), vec3(channelRadius, -channelRadius, channelRadius), vec3(-channelRadius, -channelRadius, -channelRadius), vec3(channelRadius, -channelRadius, -channelRadius),
							vec3(-channelRadius, channelRadius, channelRadius), vec3(channelRadius, channelRadius, channelRadius), vec3(-channelRadius, -channelRadius, channelRadius), vec3(channelRadius, -channelRadius, channelRadius),
							vec3(channelRadius, channelRadius, -channelRadius), vec3(channelRadius, channelRadius, channelRadius), vec3(channelRadius, -channelRadius, -channelRadius), vec3(channelRadius, -channelRadius, channelRadius),
							vec3(channelRadius, channelRadius, -channelRadius), vec3(-channelRadius, channelRadius, -channelRadius), vec3(channelRadius, -channelRadius, -channelRadius), vec3(-channelRadius, -channelRadius, -channelRadius),
							vec3(-channelRadius, channelRadius, channelRadius), vec3(-channelRadius, channelRadius, -channelRadius), vec3(-channelRadius, -channelRadius, channelRadius), vec3(-channelRadius, -channelRadius, -channelRadius)
	);
	const vec2 texCoords[] = vec2[] (vec2(0.0f, 1.0f), vec2(1.0f, 1.0f), vec2(0.0f, 0.0f), vec2(1.0f, 0.0f));
	
	for (int i = 0; i < 6; ++i) {
		for (int j = 0; j < 4; ++j) {
			geomChannelState = vertChannelState[0];
			gl_Position = viewMatrix * (gl_in[0].gl_Position + vec4(offset[i * 4 + j], 0.0f));
			texCoord = texCoords[j];
			EmitVertex();
		}
		EndPrimitive();
	}
}