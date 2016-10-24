
// compile with
// nvcc -fatbin -arch=compute_20 "kernelSource/test.cu" -o "kernelBinaries/test.bin"

extern "C" __global__ void add(int n, double *a, double *b, double *out) {

	int i = blockIdx.x*blockDim.x + threadIdx.x;
	
	if (i < n) {
		out[i] = a[i] + b[i];
	}
}
