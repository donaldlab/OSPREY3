
extern "C" __global__ void add(int n, float *a, float *b, float *out) {

	int i = blockIdx.x * blockDim.x + threadIdx.x;
	
	if (i < n) {
		out[i] = a[i] + b[i];
	}
}
