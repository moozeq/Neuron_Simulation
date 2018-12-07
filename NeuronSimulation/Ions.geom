#version 450 core

layout (points) in;
layout (triangle_strip, max_vertices = 4) out;

uniform float ionRadius;

out vec2 texCoord;

void main(void) {
	const vec2 offset[] = vec2[] (vec2(-ionRadius, ionRadius), vec2(ionRadius, ionRadius), vec2(-ionRadius, -ionRadius), vec2(ionRadius, -ionRadius));
	const vec2 texCoords[] = vec2[] (vec2(0.0f, 1.0f), vec2(1.0f, 1.0f), vec2(0.0f, 0.0f), vec2(1.0f, 0.0f));
	
	for (int i = 0; i < 4; ++i) {
		gl_Position = gl_in[0].gl_Position + vec4(offset[i], 0.0f, 0.0f);
		texCoord = texCoords[i];
		EmitVertex();
	}
}