
/* compile with
	nvcc -cubin -dlink
	-gencode arch=compute_35,code=sm_35
	-gencode arch=compute_50,code=sm_50
	-gencode arch=compute_52,code=sm_52
	"kernelSource/testDP.cu" -o "kernelBinaries/testDP.bin"  -L /usr/lib/x86_64-linux-gnu
	
	see: http://docs.nvidia.com/cuda/maxwell-compatibility-guide/index.html#building-maxwell-compatible-apps-using-cuda-6-0
*/

const int NumRuns = 1000;

__global__ void add(int n, float *a, float *b, float *out) {

	int i = blockIdx.x*blockDim.x + threadIdx.x;
	
	if (i < n) {
		out[i] = a[i] + b[i];
	}
}

__device__ int calcNumBlocks(int numThreads, int blockThreads) {
	return (numThreads + blockThreads - 1)/blockThreads;
}

extern "C" __global__ void loop(int n, float *a, float *b, float *out) {
	
	// just in case
	if (threadIdx.x != 0) {
		return;
	}
	
	int blockThreads = 256;
	int gridBlocks = calcNumBlocks(n, blockThreads);
	
	for (int i=0; i<NumRuns; i++) {
		add<<<gridBlocks, blockThreads>>>(n, a, b, out);
		cudaDeviceSynchronize();
	}
}
