
#include "global.h"
#include "cuda.h"


namespace cuda {

	__host__
	void check_error() {
		cudaError_t error = cudaGetLastError();
		if (error != cudaSuccess) {
			throw std::runtime_error(cudaGetErrorString(error));
		}
	}

	__host__
	int optimize_threads_void(const void * func, size_t shared_size_static, size_t shared_size_per_thread) {

		// TODO: device selection

		int max_block_size;
		cudaDeviceGetAttribute(&max_block_size, cudaDevAttrMaxThreadsPerBlock, 0);
		check_error();

		int best_block_size = 0;
		for (int block_size=1; block_size<=max_block_size; block_size*=2) {
			int num_blocks;
			size_t shared_size = shared_size_static + shared_size_per_thread*block_size;
			cudaOccupancyMaxActiveBlocksPerMultiprocessor(&num_blocks, func, block_size, shared_size);
			check_error();
			if (num_blocks <= 0) {
				break;
			} else {
				best_block_size = block_size;
			}
		}

		return best_block_size;
	}
}
