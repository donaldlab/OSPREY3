#version 450

// in int gl_VertexIndex;

layout(location = 0) out vec3 outPosCamera;
layout(location = 1) out vec2 outRadiusClip;

#include "view.glsl"
#include "light.glsl" // NOTE: bindings 1,2 aren't used, so will get compiled out
#include "ambientOcclusion.glsl"


void main() {

	// convert the vertex index into a world pos based on the sampling grid
	ivec3 id = ivec3(
		gl_VertexIndex % inOcclusionField.samples.x,
		gl_VertexIndex/inOcclusionField.samples.x % inOcclusionField.samples.y,
		gl_VertexIndex/inOcclusionField.samples.y/inOcclusionField.samples.x
	);
	vec3 posWorld = inOcclusionField.min + vec3(id)*(inOcclusionField.max - inOcclusionField.min)/vec3(inOcclusionField.samples - ivec3(1));

	// NOTE: vertex shaders convert from world space to clip space
	// via world -> camera -> NDC -> clip

	// transform from world to camera space
	vec3 posCamera = worldToCamera(posWorld);

	// transform from camera to clip space
	vec4 posClip = cameraToClip(posCamera);
	vec2 radiusClip = cameraPerpendicularToClip(RADIUS_CAMERA, posCamera.z);

	// send outputs to next shader stages
	gl_Position = posClip;
	outPosCamera = posCamera;
	outRadiusClip = radiusClip;
}
