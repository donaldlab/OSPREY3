
#include "global.h"
#include "cuda.h"
#include "array.h"
#include "formats.h"
#include "rotation.h"
#include "atoms.h"
#include "confspace.h"
#include "assignment.h"
#include "energy.h"
#include "energy_ambereef1.h"
#include "motions.h"
#include "minimization.h"
#include "confecalc.h"


namespace osprey {

	void exception_handler() {

		// get the stack trace
		const int SIZE = 32;
		void * trace_elems[SIZE];
		int num_frames = backtrace(trace_elems, SIZE);
		char ** symbols = backtrace_symbols(trace_elems, num_frames);

		// print it to the console
		std::cerr << "Unhandled exception in CUDA conformation energy calculator native code!" << std::endl;
		for (int i=0 ; i<num_frames; i++) {
			std::cerr << "\t" << symbols[i] << std::endl;
		}

		// cleanup
		free(symbols);

		// exit the program with an error code
		std::abort();
	}

	void free_conf_space(int device, void * p) {
		CUDACHECK(cudaSetDevice(device));
		CUDACHECK(cudaFree(p));
	}

	Stream * alloc_stream(int device, int64_t host_bytes, int64_t device_bytes) {
		CUDACHECK(cudaSetDevice(device));
		return new Stream(host_bytes, device_bytes);
	}

	void free_stream(int device, Stream * stream) {
		CUDACHECK(cudaSetDevice(device));
		delete stream;
	}


	// set kernel launch bounds for each real type

	template<>
	__global__
	__launch_bounds__(MINIMIZE_THREADS_F32, MINIMIZE_BLOCKS)
	void minimize_kernel<float32_t,osprey::ambereef1::calc_energy>(int64_t shared_size,
	                                                               const ConfSpace<float32_t> * conf_space,
	                                                               int64_t max_num_inters,
	                                                               int8_t * xfer_buf) {
		minimize_kernel_impl<float32_t,osprey::ambereef1::calc_energy>(shared_size, conf_space, max_num_inters, xfer_buf);
	}

	template<>
	__global__
	__launch_bounds__(MINIMIZE_THREADS_F64, MINIMIZE_BLOCKS)
	void minimize_kernel<float64_t,osprey::ambereef1::calc_energy>(int64_t shared_size,
	                                                               const ConfSpace<float64_t> * conf_space,
	                                                               int64_t max_num_inters,
	                                                               int8_t * xfer_buf) {
		minimize_kernel_impl<float64_t,osprey::ambereef1::calc_energy>(shared_size, conf_space, max_num_inters, xfer_buf);
	}
}
