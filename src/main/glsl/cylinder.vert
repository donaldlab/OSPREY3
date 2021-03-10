#version 450

layout(location = 0) in vec3 inPosWorld[2];
layout(location = 2) in float inRadiusWorld[2];
layout(location = 4) in vec4 inColor[2];
layout(location = 6) in int inIndex[2];
layout(location = 8) in uvec4 inEffect[2];

layout(location = 0) out vec3 outPosCamera[2];
layout(location = 2) out float outRadiusCamera[2];
layout(location = 4) out vec4 outColor[2];
layout(location = 6) out int outIndex[2];
layout(location = 8) out uvec4 outEffect[2];

#include "view.glsl"


void main() {

	// NOTE: vertex shaders convert from world space to clip space
	// via world -> camera -> NDC -> clip

	// transform from world to camera space
	vec3[] posCamera = {
		worldToCamera(inPosWorld[0]),
		worldToCamera(inPosWorld[1]),
	};
	float[] radiusCamera = {
		worldToCamera(inRadiusWorld[0]),
		worldToCamera(inRadiusWorld[1])
	};

	// transform from camera to clip space
	vec4[] posClip = {
		cameraToClip(posCamera[0]),
		cameraToClip(posCamera[1])
	};
	vec2[] radiusClip = {
		cameraPerpendicularToClip(radiusCamera[0], posCamera[0].z),
		cameraPerpendicularToClip(radiusCamera[1], posCamera[1].z)
	};

	// get the line vector
	vec2 dir = normalize(posClip[1].xy - posClip[0].xy);

	// do a 90 degree ccw rotation to get a perpendicular direction
	vec2 pdir = vec2(-dir[1], dir[0]);

	// TODO: could make this tighter by doing some elliptical math
	// dunno if it's worth it though
	float r = max(
		max(radiusClip[0].x, radiusClip[1].x),
		max(radiusClip[0].y, radiusClip[1].y)
	);

	// offset vertices by the radius to make a billboard quad, in triangle strip order, with ccw facing
	int i = gl_VertexIndex % 4;
	if (i == 0) {
		gl_Position = posClip[0] + vec4(0 - dir*r - pdir*r, 0, 0);
	} else if (i == 1) {
		gl_Position = posClip[0] + vec4(0 - dir*r + pdir*r, 0, 0);
	} else if (i == 2) {
		gl_Position = posClip[1] + vec4(0 + dir*r - pdir*r, 0, 0);
	} else {
		gl_Position = posClip[1] + vec4(0 + dir*r + pdir*r, 0, 0);
	}

	// send outputs to next shader stages
	outPosCamera = posCamera;
	outRadiusCamera = inRadiusWorld;
	outColor = inColor;
	outIndex = inIndex;
	outEffect = inEffect;
}
