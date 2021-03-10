#version 450

// a single quad (as a triangle strip) covering all of the xy clip plane
vec2 positions[4] = vec2[] (
	vec2(-1.0, 1.0),
	vec2(1.0, 1.0),
	vec2(-1.0, -1.0),
	vec2(1.0, -1.0)
);

void main() {
	gl_Position = vec4(positions[gl_VertexIndex], 0.0, 1.0);
}
