
#include "global.h"
#include "cuda.h"


namespace cuda {

	__host__
	int get_arch() {

		int major;
		CUDACHECK(cudaDeviceGetAttribute(&major, cudaDevAttrComputeCapabilityMajor, 0));

		int minor;
		CUDACHECK(cudaDeviceGetAttribute(&minor, cudaDevAttrComputeCapabilityMinor, 0));

		return major*100 + minor*10;
	}

	__host__
	int get_max_shared_size(int device) {
		int val;
		cudaDeviceGetAttribute(&val, cudaDevAttrMaxSharedMemoryPerBlock, device);
		return val;
	}

	__host__
	int optimize_threads_void(const void * func, size_t shared_size_static, size_t shared_size_per_thread) {

		int max_block_size;
		CUDACHECK(cudaDeviceGetAttribute(&max_block_size, cudaDevAttrMaxThreadsPerBlock, 0));

		int best_block_size = 0;
		for (int block_size=1; block_size<=max_block_size; block_size*=2) {
			int num_blocks;
			size_t shared_size = shared_size_static + shared_size_per_thread*block_size;
			CUDACHECK(cudaOccupancyMaxActiveBlocksPerMultiprocessor(&num_blocks, func, block_size, shared_size));
			if (num_blocks <= 0) {
				break;
			} else {
				best_block_size = block_size;
			}
		}

		return best_block_size;
	}
}
