#version 450

layout(location = 0) in vec3 inPosWorld;
layout(location = 1) in float inRadiusWorld;
layout(location = 2) in vec4 inColor;
layout(location = 3) in int inIndex;
layout(location = 4) in uvec4 inEffect;

layout(location = 0) out vec3 outPosCamera;
layout(location = 1) out float outRadiusCamera;
layout(location = 2) out vec4 outColor;
layout(location = 3) out int outIndex;
layout(location = 4) out uvec4 outEffect;

#include "view.glsl"


void main() {

	// NOTE: vertex shaders convert from world space to clip space
	// via world -> camera -> NDC -> clip

	// transform from world to camera space
	vec3 posCamera = worldToCamera(inPosWorld);
	float radiusCamera = worldToCamera(inRadiusWorld);

	// transform from camera to clip space
	vec4 posClip = cameraToClip(posCamera);
	vec2 radiusClip = cameraPerpendicularToClip(radiusCamera, posCamera.z);

	// offset vertices by the radius to make a billboard quad, in triangle strip order, with ccw facing
	if ((gl_VertexIndex & 0x02) != 0) {
		radiusClip.x = -radiusClip.x;
	}
	if ((gl_VertexIndex & 0x01) != 0) {
		radiusClip.y = -radiusClip.y;
	}
	posClip += vec4(radiusClip, 0, 0);

	// send outputs to next shader stages
	gl_Position = posClip;
	outPosCamera = posCamera;
	outRadiusCamera = inRadiusWorld;
	outColor = inColor;
	outIndex = inIndex;
	outEffect = inEffect;
}
