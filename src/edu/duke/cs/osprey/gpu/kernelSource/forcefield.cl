
#ifdef cl_khr_fp64
    #pragma OPENCL EXTENSION cl_khr_fp64 : enable
#elif defined(cl_amd_fp64)
    #pragma OPENCL EXTENSION cl_amd_fp64 : enable
#else
    #error "Double precision floating point not supported by OpenCL implementation."
#endif

kernel void calc(global const double *in, global double *out) {
	int i = get_global_id(0);
	
	// simulate a force-field type workload
	// so I can test speed without having to implement the whole shebang just yet
	
	// do a few lookups
	float a = in[i];
	float b = in[0];
	float c = in[4];
	float d = in[5];
	float e = in[7];
	float f = in[10];
	
	// do some math
	float g = a*b + c*d + e*f;
	
	// do one sqrt
	g = sqrt(g);
	
	// do some more math
	g = (a + b)/g; 
	
	out[i] = g;
}
