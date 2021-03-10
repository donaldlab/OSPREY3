#version 450

layout(location = 0) flat in vec3 inPosCamera;
layout(location = 1) flat in float inRadiusCamera;
layout(location = 2) flat in vec4 inColor;
layout(location = 3) flat in int inIndex;
layout(location = 4) flat in uvec4 inEffect;

layout(location = 0) out vec4 outColor;
layout(location = 1) out ivec2 outIndices;
layout(location = 2) out uvec4 outEffect;
layout(depth_less) out; // float gl_FragDepth

layout(push_constant) uniform ViewIndex {
	int inViewIndex;
};

#include "view.glsl"
#include "light.glsl"


float intersectSphereRay(vec3 p, float r, vec2 ray) {

	// see math/sphere.wxm for the derivation
	return p.z-sqrt(-ray.y*ray.y+2*p.y*ray.y-ray.x*ray.x+2*p.x*ray.x+r*r-p.y*p.y-p.x*p.x);
}


void main() {

	// transform from framebuf to camera space
	vec3 posPixelCamera = framebufToCamera(gl_FragCoord);

	// compute the pixel z pos on the sphere, in camera space
	float z = intersectSphereRay(inPosCamera, inRadiusCamera, posPixelCamera.xy);
	if (!isnan(z)) {

		// pixel is on or inside the sphere, update the camera pos
		posPixelCamera.z = z;

		// calc the sphere normal
		vec3 normal = normalize(posPixelCamera - inPosCamera);

		outColor = light(inColor, posPixelCamera, normal);
		outIndices = ivec2(inIndex, inViewIndex);
		outEffect = inEffect;
		gl_FragDepth = cameraZToClip(posPixelCamera.z);

	} else {

		// pixel is outside the sphere
		outColor = vec4(0);
		outIndices = ivec2(-1, -1);
		outEffect = ivec4(0, 0, 0, 0);
		gl_FragDepth = 1;
	}
}
