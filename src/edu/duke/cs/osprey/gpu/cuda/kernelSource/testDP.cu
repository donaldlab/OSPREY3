
/* compile with
	nvcc -cubin -dlink
	-gencode arch=compute_35,code=sm_35
	-gencode arch=compute_50,code=sm_50
	-gencode arch=compute_52,code=sm_52
	"kernelSource/testDP.cu" -o "kernelBinaries/testDP.bin"  -L /usr/lib/x86_64-linux-gnu
	
	see: http://docs.nvidia.com/cuda/maxwell-compatibility-guide/index.html#building-maxwell-compatible-apps-using-cuda-6-0
*/

const int NumElements = 10000;
const int NumRuns = 10;

const int BlockThreads = 256;
const int GridBlocks = 40;//(NumElements + BlockThreads - 1)/BlockThreads;

extern "C" __global__ void add(int n, float *a, float *b, float *out) {

	int i = blockIdx.x*blockDim.x + threadIdx.x;
	
	if (i < n) {
		out[i] = a[i] + b[i];
	}
}


extern "C" __global__ void loop(int n, float *a, float *b, float *out) {
	
	if (threadIdx.x == 0) {
	
		for (int i=0; i<NumRuns; i++) {
			add<<<GridBlocks, BlockThreads>>>(n, a, b, out);
			cudaDeviceSynchronize();
		}
	}
}
