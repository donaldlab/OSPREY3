
#ifdef cl_khr_fp64
    #pragma OPENCL EXTENSION cl_khr_fp64 : enable
#elif defined(cl_amd_fp64)
    #pragma OPENCL EXTENSION cl_amd_fp64 : enable
#else
    #error "Double precision floating point not supported by OpenCL implementation."
#endif


kernel void addFloat(global const float *a, global const float *b, global float *out) {
	int i = get_global_id(0);
	out[i] = a[i] + b[i];
}

kernel void addDouble(global const double *a, global const double *b, global double *out) {
	int i = get_global_id(0);
	out[i] = a[i] + b[i];
}
