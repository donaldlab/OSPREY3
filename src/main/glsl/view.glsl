
#ifndef _VIEW_GLSL_
#define _VIEW_GLSL_


layout(binding = 0, std140) uniform restrict readonly View {
	vec3 cameraWorld;
	// 4 bytes pad
	vec3 cameraSide;
	// 4 bytes pad
	vec3 cameraUp;
	// 4 bytes pad
	vec3 cameraLook;
	// 4 bytes pad
	vec2 windowSize;
	float zNearCamera;
	float zFarCamera;
	float magnification;
} view;


//=========================
// direction transforms
//=========================

// camera space:
// z points in the direction the camera is looking
// y points towards the up direction
// x is defined by applying the right-hand rule on y and z

vec3 dirWorldToCamera(vec3 v) {
	return vec3(
		dot(view.cameraSide, v),
		dot(view.cameraUp, v),
		dot(view.cameraLook, v)
	);
}

vec3 dirCameraToWorld(vec3 v) {
	return vec3(
		dot(vec3(view.cameraSide.x, view.cameraUp.x, view.cameraLook.x), v),
		dot(vec3(view.cameraSide.y, view.cameraUp.y, view.cameraLook.y), v),
		dot(vec3(view.cameraSide.z, view.cameraUp.z, view.cameraLook.z), v)
	);
}


//=========================
// position transforms
//=========================

vec3 worldToCamera(vec3 v) {
	return dirWorldToCamera(v - view.cameraWorld);
}

vec3 cameraToWorld(vec3 v) {
	return dirCameraToWorld(v) + view.cameraWorld;
}


// NDC space (normalized device coords) for Vulkan:
// x in [-1,1]   y in [-1,1]   z in [0,1]
// right-handed axis system: x+right, y+down, z+away
// camera -> NDC is perspective projection

// on world coords, if z is away, and y is up, RHR puts x pointing left
// going from camera to NDC coords means moving:
//   y from up to down
//   x from left to right
// so negate both x and y components together

vec3 cameraToNDC(vec3 v) {
	vec2 viewSize = view.windowSize/view.magnification;
	return vec3(
		-2.0*v.xy*view.zNearCamera/viewSize/v.z,
		(v.z - view.zNearCamera)/(view.zFarCamera - view.zNearCamera)
	);
}

vec3 NDCToCamera(vec3 v) {
	vec2 viewSize = view.windowSize/view.magnification;
	float zCamera = v.z*(view.zFarCamera - view.zNearCamera) + view.zNearCamera;
	return vec3(
		-v.xy*zCamera*viewSize/2/view.zNearCamera,
		zCamera
	);
}


// clip space:
// same as NDC space, but in homogeneous coords
vec4 NDCToClip(vec3 v) {
	return vec4(v, 1);
}

vec3 clipToNDC(vec4 v) {
	return v.xyz/v.w;
}


// framebuffer space:
// units are in pixels
// y axis is down!

vec3 framebufToNDC(vec3 v) {
	return vec3(
		v.xy*2/view.windowSize - vec2(1),
		v.z
	);
}


// combination transforms
vec4 cameraToClip(vec3 v) {
	return NDCToClip(cameraToNDC(v));
}

vec3 framebufToCamera(vec4 v) {
	return NDCToCamera(framebufToNDC(v.xyz));
}


//=========================
// distance transforms
//=========================


float worldToCamera(float d) {
	// camera transformation always preserves distances
	return d;
}

vec2 cameraPerpendicularToClip(float d, float z) {
	vec2 viewSize = view.windowSize/view.magnification;
	return 2*view.zNearCamera/viewSize/z*vec2(d);
}

float cameraZToClip(float z) {
	return (z - view.zNearCamera)/(view.zFarCamera - view.zNearCamera);
}


#endif // _VIEW_GLSL_
