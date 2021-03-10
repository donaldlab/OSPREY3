#version 450

layout(location = 0) flat in vec3 inPosCamera[2];
layout(location = 2) flat in float inRadiusCamera[2];
layout(location = 4) flat in vec4 inColor[2];
layout(location = 6) flat in int inIndex[2];
layout(location = 8) flat in uvec4 inEffect[2];

layout(location = 0) out vec4 outColor;
layout(location = 1) out ivec2 outIndices;
layout(location = 2) out uvec4 outEffect;
layout(depth_any) out; // float gl_FragDepth

layout(push_constant) uniform ViewIndex {
	int inViewIndex;
};

#include "view.glsl"
#include "math.glsl"
#include "light.glsl"


float calcZEndcap(vec2 xy, uint i, vec3 normals[2]) {

	// intersect the ray with the plane
	float z = dot(inPosCamera[i].xy - xy, normals[i].xy)/normals[i].z + inPosCamera[i].z;

	// do a range check
	bool inRange = distance(vec3(xy, z), inPosCamera[i]) <= inRadiusCamera[i];
	if (inRange) {
		return z;
	}

	return NAN;
}

// p, n are the cylinder point and normalized direction
// r is the cylinder radius
float[2] intersectCylinderRay(vec3 p, vec3 n, float r, float len, vec2 ray) {

	// see math/cylinder.wxm for the derivation
	// this got manually optimized a bit more too
	float a1 = n.x*n.x;
	float a2 = n.y*n.y;
	float a3 = a2+a1;
	float a5 = (-a2-a1)*p.z;
	float a6 = sqrt(-a1*ray.y*ray.y+(2*n.x*n.y*ray.x+2*a1*p.y-2*n.x*n.y*p.x)*ray.y-a2*ray.x*ray.x+(2*a2*p.x-2*n.x*n.y*p.y)*ray.x+a3*r*r-a1*p.y*p.y+2*n.x*n.y*p.x*p.y-a2*p.x*p.x);
	float a7 = -n.z*(n.y*ray.y + n.x*ray.x - n.y*p.y - n.x*p.x);
	float z = -(a6+a7+a5)/a3;

	// get the normalized distance along the cylindrical axis
	float t = dot(vec3(ray, z) - p, n)/len;

	// if we're out of range, drop the intersection entirely
	if (t < 0 || t > 1) {
		z = NAN;
		t = NAN;
	}

	float result[2] = { z, t };
	return result;
}

void main() {

	// transform from framebuf to camera space
	vec3 posPixelCamera = framebufToCamera(gl_FragCoord);

	// convert to a pos,normal cylinder representation
	vec3 posCylinder = inPosCamera[0];
	vec3 axisCylinder = inPosCamera[1] - inPosCamera[0];
	float len = length(axisCylinder);
	axisCylinder /= len;

	// intersect with the cylindrical surface
	float intersection[2] = intersectCylinderRay(posCylinder, axisCylinder, inRadiusCamera[0], len, posPixelCamera.xy);
	float zCylinder = intersection[0];
	float tCylinder = intersection[1];

	// if we have a valid point on the cylinder, use that first
	if (!isnan(zCylinder)) {

		posPixelCamera.z = zCylinder;

		// compute the normal
		vec3 center = posCylinder + tCylinder*len*axisCylinder;
		vec3 normal = normalize(posPixelCamera - center);

		// pick which end to use based on the cylinder parameter
		uint iEnd = tCylinder <= 0.5 ? 0 : 1;

		outColor = light(inColor[iEnd], posPixelCamera, normal);
		outIndices = ivec2(inIndex[iEnd], inViewIndex);
		outEffect = inEffect[iEnd];
		gl_FragDepth = cameraZToClip(posPixelCamera.z);

	} else {

		// otherwise, find out which endcap (if any) intersects
		vec3 normalEndcap[2] = {
			-axisCylinder,
			axisCylinder
		};
		float zEndcap[2] = {
			calcZEndcap(posPixelCamera.xy, 0, normalEndcap),
			calcZEndcap(posPixelCamera.xy, 1, normalEndcap)
		};

		// did we hit any endcaps?
		int iEndcap = -1;
		if (!isnan(zEndcap[0]) && !isnan(zEndcap[1])) {

			// got both, pick the closer one
			if (zEndcap[0] < zEndcap[1]) {
				iEndcap = 0;
			} else {
				iEndcap = 1;
			}

		} else if (!isnan(zEndcap[0])) {
			iEndcap = 0;
		} else if (!isnan(zEndcap[1])) {
			iEndcap = 1;
		}

		if (iEndcap >= 0) {

			// yup, found an endcap
			posPixelCamera.z = zEndcap[iEndcap];

			outColor = light(inColor[iEndcap], posPixelCamera, normalEndcap[iEndcap]);
			outIndices = ivec2(inIndex[iEndcap], inViewIndex);
			outEffect = inEffect[iEndcap];
			gl_FragDepth = cameraZToClip(posPixelCamera.z);

		} else {

			// no intersection
			outColor = vec4(0);
			outIndices = ivec2(-1, -1);
			outEffect = ivec4(0, 0, 0, 0);
			gl_FragDepth = 1;
		}
	}
}
