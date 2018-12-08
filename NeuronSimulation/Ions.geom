#version 450 core

layout (points) in;
layout (triangle_strip, max_vertices = 24) out;

uniform float ionRadius;
uniform mat4 viewMatrix;

out vec2 texCoord;

void main(void) {
	const vec3 offset[] = vec3[] (
							vec3(-ionRadius, ionRadius, ionRadius), vec3(ionRadius, ionRadius, ionRadius), vec3(-ionRadius, ionRadius, -ionRadius), vec3(ionRadius, ionRadius, -ionRadius),
							vec3(-ionRadius, -ionRadius, ionRadius), vec3(ionRadius, -ionRadius, ionRadius), vec3(-ionRadius, -ionRadius, -ionRadius), vec3(ionRadius, -ionRadius, -ionRadius),
							vec3(-ionRadius, ionRadius, ionRadius), vec3(ionRadius, ionRadius, ionRadius), vec3(-ionRadius, -ionRadius, ionRadius), vec3(ionRadius, -ionRadius, ionRadius),
							vec3(ionRadius, ionRadius, -ionRadius), vec3(ionRadius, ionRadius, ionRadius), vec3(ionRadius, -ionRadius, -ionRadius), vec3(ionRadius, -ionRadius, ionRadius),
							vec3(ionRadius, ionRadius, -ionRadius), vec3(-ionRadius, ionRadius, -ionRadius), vec3(ionRadius, -ionRadius, -ionRadius), vec3(-ionRadius, -ionRadius, -ionRadius),
							vec3(-ionRadius, ionRadius, ionRadius), vec3(-ionRadius, ionRadius, -ionRadius), vec3(-ionRadius, -ionRadius, ionRadius), vec3(-ionRadius, -ionRadius, -ionRadius)
	);
	const vec2 texCoords[] = vec2[] (vec2(0.0f, 1.0f), vec2(1.0f, 1.0f), vec2(0.0f, 0.0f), vec2(1.0f, 0.0f));
	
	for (int i = 0; i < 6; ++i) {
		for (int j = 0; j < 4; ++j) {
			gl_Position = viewMatrix * (gl_in[0].gl_Position + vec4(offset[i * 4 + j], 0.0f));
			texCoord = texCoords[j];
			EmitVertex();
		}
		EndPrimitive();
	}
}