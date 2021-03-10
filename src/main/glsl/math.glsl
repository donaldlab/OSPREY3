
#ifndef _MATH_GLSL_
#define _MATH_GLSL_


const float NAN = 0.0/0.0;
const float PI = 3.141592654;
const float TwoPI = 2*PI;


float length2(vec3 a) {
	return dot(a, a);
}

float distance2(vec3 a, vec3 b) {
	return length2(a - b);
}


#endif // _MATH_GLSL_
