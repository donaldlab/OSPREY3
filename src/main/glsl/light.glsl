
#ifndef _LIGHT_GLSL_
#define _LIGHT_GLSL_


#include "math.glsl"


layout(binding = 1) uniform isampler3D inOcclusion;

layout(binding = 2, std140) uniform readonly restrict OcclusionField {
	// (with padding to align vec3s at 16-byte boundaries)
	ivec3 samples;
	int maxOcclusion;
	vec3 min; // in world space
	float minOcclusion; // in [0,1], offset = 28
	vec3 max;
	float pad3;
} inOcclusionField;

layout(binding = 3, std140) uniform readonly restrict Settings {
	vec3 backgroundColor;
	float colorWeight;
	float shadingWeight;
	float lightWeight;
	float depthWeight;
	float depthZMin;
	float depthZMax;
	float ambientOcclusionWeight;
} inSettings;


const vec3 toLight = normalize(vec3(1, 1, -1));

vec3 pbrLambert(vec3 albedo, vec3 normal, float intensity) {

	vec3 a = albedo;
	vec3 n = normal;
	vec3 l = toLight;
	float i = intensity;

	// see math/lighting-lambert.wxm for derivation
	return (a*i*(2*PI*l.z*n.z+2*PI*l.y*n.y+2*PI*l.x*n.x+3*PI))/(6*PI);
}

float sampleOcclusion(vec3 posField) {

	// convert the position to texture coords
	// scale the field coords to the centers of the pixels
	vec3 pixelSize = vec3(1)/vec3(inOcclusionField.samples);
	vec3 uvw = posField*(1 - pixelSize)/(vec3(inOcclusionField.samples) - vec3(1)) + pixelSize/2;

	return float(texture(inOcclusion, uvw).r)/float(inOcclusionField.maxOcclusion);
}

float nearestOcclusion(vec3 posField) {
	return sampleOcclusion(vec3(ivec3(posField + vec3(0.5))));
}

vec3 worldToOcclusionField(vec3 posWorld) {
	vec3 posField = posWorld - inOcclusionField.min;
	posField *= ivec3(inOcclusionField.samples) - ivec3(1);
	posField /= inOcclusionField.max - inOcclusionField.min;
	return posField;
}

vec4 light(vec4 color, vec3 posCamera, vec3 normalCamera) {

	vec3 rgb = vec3(1, 1, 1);

	// apply color
	if (inSettings.colorWeight > 0) {
		rgb = mix(rgb, color.rgb, inSettings.colorWeight);
	}

	// apply shading
	if (inSettings.shadingWeight > 0) {
		rgb = mix(rgb, pbrLambert(rgb, normalCamera, inSettings.lightWeight), inSettings.shadingWeight);
	}

	// apply depth fading (ie far away things are more like the background)
	if (inSettings.depthWeight > 0) {
		float depth = cameraToNDC(posCamera).z;
		float width = inSettings.depthZMax - inSettings.depthZMin;
		if (width < 0) {
			width = 0;
		}
		depth = (depth - inSettings.depthZMin)/width;
		depth = clamp(depth, 0, 1);
		rgb = mix(rgb, inSettings.backgroundColor, inSettings.depthWeight*depth);
	}

	// apply ambient occlusion
	if (inSettings.ambientOcclusionWeight > 0) {
		float occlusion = sampleOcclusion(worldToOcclusionField(cameraToWorld(posCamera)));

		// normalize occlusion so the min is actually 0
		if (inOcclusionField.minOcclusion >= 1) {
			occlusion = 0;
		} else {
			occlusion = max(0, occlusion - inOcclusionField.minOcclusion)/(1 - inOcclusionField.minOcclusion);
		}

		rgb *= 1 - inSettings.ambientOcclusionWeight*occlusion;
	}
	
	return vec4(rgb, color.a);
}

vec4 showNormal(vec3 normal) {
	return vec4((normal.xy + vec2(1))/2, 1 - (normal.z + 1)/2, 1);
}

vec4 showZ(float zClip) {
	return vec4(vec3(1 - zClip), 1);
}


#endif // _LIGHT_GLSL
