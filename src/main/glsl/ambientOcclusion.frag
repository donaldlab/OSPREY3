#version 450

layout(location = 0) flat in vec3 inPosCamera;

layout(location = 0) out vec4 outColor;
layout(location = 1) out ivec2 outIndex;
layout(location = 2) out uvec4 outEffects;
layout(depth_less) out; // float gl_FragDepth

layout(push_constant) uniform ViewIndex {
	int inViewIndex;
};

#include "view.glsl"
#include "light.glsl"
#include "ambientOcclusion.glsl"


float intersectSphereRay(vec3 p, float r, vec2 ray) {

	// see math/sphere.wxm for the derivation
	return p.z-sqrt(-ray.y*ray.y+2*p.y*ray.y-ray.x*ray.x+2*p.x*ray.x+r*r-p.y*p.y-p.x*p.x);
}


// list cardinal direction indices
// NOTE: negative,positive pairs must be consecutive and in that order
const uint DIR_XN = 0;
const uint DIR_XP = 1;
const uint DIR_YN = 2;
const uint DIR_YP = 3;
const uint DIR_ZN = 4;
const uint DIR_ZP = 5;
const uint NUM_DIRS = 6;

const vec3[] DIRS = {
	vec3(-1.0, 0.0, 0.0),
	vec3(+1.0, 0.0, 0.0),
	vec3(0.0, -1.0, 0.0),
	vec3(0.0, +1.0, 0.0),
	vec3(0.0, 0.0, -1.0),
	vec3(0.0, 0.0, +1.0)
};

uint getDir(vec3 v) {

	float bestd = dot(v, DIRS[0]);
	uint iDir = 0;

	// find the direction with the greatest dot product
	for (uint i=1; i<NUM_DIRS; i++) {
		float d = dot(v, DIRS[i]);
		if (d > bestd) {
			bestd = d;
			iDir = i;
		}
	}

	return iDir;
}


void main() {

	// transform from framebuf to camera space
	vec3 posPixelCamera = framebufToCamera(gl_FragCoord);

	// compute the pixel z pos on the sphere, in camera space
	float z = intersectSphereRay(inPosCamera, RADIUS_CAMERA, posPixelCamera.xy);
	if (!isnan(z)) {

		// pixel is on or inside the sphere, update the camera pos
		posPixelCamera.z = z;

		// get occlusion for this position
		float occlusion = nearestOcclusion(worldToOcclusionField(cameraToWorld(inPosCamera)));
		if (occlusion > 0) {

			// convert the occlusion to a color
			vec3 color;
			if (occlusion >= 0.5) {
				color = mix(vec3(1, 0, 0), vec3(0, 0, 0), (occlusion - 0.5)*2);
			} else {
				color = mix(vec3(1, 1, 1), vec3(1, 0, 0), occlusion*2);
			}

			// output the pixel
			outColor = vec4(color, 1);
			gl_FragDepth = cameraZToClip(posPixelCamera.z);

		} else {

			// drop the pixel
			outColor = vec4(0);
			gl_FragDepth = 1;
		}

	} else {

		// pixel is outside the sphere
		outColor = vec4(0);
		gl_FragDepth = 1;
	}

	outIndex = ivec2(-1, -1);
	outEffects = uvec4(0, 0, 0, 0);
}
