
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

#define API extern "C" [[maybe_unused]]


API int version_major() {
	return CudaConfEcalc_VERSION_MAJOR;
}

API int version_minor() {
	return CudaConfEcalc_VERSION_MINOR;
}

API int cuda_version_driver() {
	int v;
	cudaDriverGetVersion(&v);
	return v;
}

API int cuda_version_runtime() {
	int v;
	cudaRuntimeGetVersion(&v);
	switch (cudaGetLastError()) {
		case cudaSuccess: break;
		case cudaErrorInsufficientDriver: return -1;
		case cudaErrorNoDevice: return -2;
		default: return INT_MIN;
	}
	return v;
}

API int cuda_version_required() {
	return CUDART_VERSION;
}

struct GpuInfo {
	char bus_id[16];
	char name[256];
	int integrated;
	int concurrent_kernels;
	int num_processors;
	int num_async_engines;
	int64_t mem_total;
	int64_t mem_free;
};
ASSERT_JAVA_COMPATIBLE(GpuInfo, 304);

API osprey::Array<GpuInfo> * alloc_gpu_infos() {

	int num_gpus;
	CUDACHECK(cudaGetDeviceCount(&num_gpus));

	osprey::Array<GpuInfo> * infos = osprey::Array<GpuInfo>::make(num_gpus);
	for (int device=0; device<num_gpus; device++) {
		GpuInfo & info = (*infos)[device];

		CUDACHECK(cudaDeviceGetPCIBusId(reinterpret_cast<char *>(&info.bus_id), 16, device));

		cudaDeviceProp props {};
		CUDACHECK(cudaGetDeviceProperties(&props, device));

		memcpy(&info.name, props.name, 256);
		info.integrated = props.integrated;
		info.concurrent_kernels = props.concurrentKernels;
		info.num_processors = props.multiProcessorCount;
		info.num_async_engines = props.asyncEngineCount;

		size_t mem_free;
		size_t mem_total;
		CUDACHECK(cudaMemGetInfo(&mem_free, &mem_total));

		info.mem_free = mem_free;
		info.mem_total = mem_total;
	}

	return infos;
}

API void free_gpu_infos(osprey::Array<GpuInfo> * p) {
	delete[] p;
}


namespace osprey {

	template<typename T>
	static void * alloc_conf_space(int device, const ConfSpace<T> & conf_space) {

		CUDACHECK(cudaSetDevice(device));

		// allocate the device memory
		void * p_device;
		CUDACHECK(cudaMalloc(&p_device, conf_space.size));

		// upload the conf space
		CUDACHECK(cudaMemcpy(p_device, &conf_space, conf_space.size, cudaMemcpyHostToDevice));

		return p_device;
	}
}

API void * alloc_conf_space_f32(int device, const osprey::ConfSpace<float32_t> & conf_space) {
	return osprey::alloc_conf_space(device, conf_space);
}
API void * alloc_conf_space_f64(int device, const osprey::ConfSpace<float64_t> & conf_space) {
	return osprey::alloc_conf_space(device, conf_space);
}

API void free_conf_space(int device, void * p) {

	CUDACHECK(cudaSetDevice(device));

	CUDACHECK(cudaFree(p));
}


class Stream {
	public:

		static const size_t d_buf_size = 2*1024*1024; // 8 MiB, big enough to hold all the coords
		static const size_t h_buf_size = 16; // just big enough to hold the resulting energy

		__host__
		inline Stream() {
			CUDACHECK(cudaStreamCreate(&stream));
			CUDACHECK(cudaMallocHost(&h_buf, h_buf_size));
			CUDACHECK(cudaMalloc(&d_buf, d_buf_size));
		}

		Stream(const Stream & other) = delete;

		__host__
		inline ~Stream() {
			cudaFreeHost(h_buf);
			cudaFree(d_buf);
			cudaStreamDestroy(stream);
		}

		inline cudaStream_t get_stream() const {
			return stream;
		}

		inline int8_t * get_h_buf() {
			return h_buf;
		}

		inline int8_t * get_d_buf() {
			return d_buf;
		}

	private:
		cudaStream_t stream;
		int8_t * h_buf;
		int8_t * d_buf;
};

API Stream * alloc_stream(int device) {
	CUDACHECK(cudaSetDevice(device));
	return new Stream();
}

API void free_stream(int device, Stream * stream) {
	CUDACHECK(cudaSetDevice(device));
	delete stream;
}


namespace osprey {

	template<typename T>
	__global__
	void assign_kernel(const ConfSpace<T> * conf_space, const Array<int32_t> * conf, Array<osprey::Real3<T>> * out_coords) {

		// slice up the shared memory
		extern __shared__ int8_t shared[];
		int64_t shared_offset = 0;
		auto shared_atom_pairs = reinterpret_cast<const void **>(shared + shared_offset);
		shared_offset += Assignment<T>::sizeof_atom_pairs(conf_space->num_pos);
		auto shared_conf_energies = reinterpret_cast<T *>(shared + shared_offset);

		// init the sizes for cudaMalloc'd Array instances
		out_coords->init(conf_space->max_num_conf_atoms);

		Assignment<T> assignment(*conf_space, *conf, *out_coords, shared_atom_pairs, shared_conf_energies);
	}

	template<typename T>
	static void assign(int device,
	                   Stream * stream,
	                   const ConfSpace<T> * d_conf_space,
	                   const Array<int32_t> & conf,
	                   Array<Real3<T>> & out_coords) {

		CUDACHECK(cudaSetDevice(device));

		// upload the arguments
		size_t conf_bytes = conf.get_bytes();
		Array<int32_t> * d_conf;
		CUDACHECK(cudaMalloc(&d_conf, conf_bytes));
		CUDACHECK(cudaMemcpy(d_conf, &conf, conf_bytes, cudaMemcpyHostToDevice));

		// allocate space for the coords
		size_t out_coords_bytes = out_coords.get_bytes();
		Array<Real3<T>> * d_out_coords;
		CUDACHECK(cudaMalloc(&d_out_coords, out_coords_bytes));

		// compute the shared memory size
		int64_t shared_size = 0
			+ Assignment<T>::sizeof_atom_pairs(conf.get_size())
			+ Assignment<T>::sizeof_conf_energies(conf.get_size());

		// launch the kernel
		// TODO: optimize launch bounds
		int num_threads = cuda::optimize_threads(assign_kernel<T>, shared_size, 0);
		assign_kernel<<<1, num_threads, shared_size, stream->get_stream()>>>(d_conf_space, d_conf, d_out_coords);
		cuda::check_error();

		// download the coords
		CUDACHECK(cudaMemcpy(&out_coords, d_out_coords, out_coords_bytes, cudaMemcpyDeviceToHost));

		// cleanup
		CUDACHECK(cudaFree(d_conf));
		CUDACHECK(cudaFree(d_out_coords));
	}
}

API void assign_f32(int device,
                    Stream * stream,
                    const osprey::ConfSpace<float32_t> * d_conf_space,
                    const osprey::Array<int32_t> & conf,
                    osprey::Array<osprey::Real3<float32_t>> & out_coords) {
	osprey::assign(device, stream, d_conf_space, conf, out_coords);
}
API void assign_f64(int device,
                    Stream * stream,
                    const osprey::ConfSpace<float64_t> * d_conf_space,
                    const osprey::Array<int32_t> & conf,
                    osprey::Array<osprey::Real3<float64_t>> & out_coords) {
	osprey::assign(device, stream, d_conf_space, conf, out_coords);
}


namespace osprey {

	#define CALC_THREADS 1024
	#define CALC_BLOCKS 1

	template<typename T, osprey::EnergyFunction<T> efunc>
	__global__
	__launch_bounds__(CALC_THREADS, CALC_BLOCKS)
	void calc_kernel(const osprey::ConfSpace<T> * conf_space,
	                 const osprey::Array<int32_t> * conf,
	                 const osprey::Array<osprey::PosInter<T>> * inters,
	                 osprey::Array<osprey::Real3<T>> * out_coords,
	                 T * out_energy) {

		// slice up the shared memory
		extern __shared__ int8_t _shared[];
		auto shared = _shared;

		auto shared_atom_pairs = reinterpret_cast<const void **>(shared);
		assert (cuda::is_aligned<16>(shared_atom_pairs));
		shared += cuda::pad_to_alignment<16>(Assignment<T>::sizeof_atom_pairs(conf_space->num_pos));

		auto shared_conf_energies = reinterpret_cast<T *>(shared);
		assert (cuda::is_aligned<16>(shared_conf_energies));
		shared += cuda::pad_to_alignment<16>(Assignment<T>::sizeof_conf_energies(conf_space->num_pos));

		auto thread_energy = reinterpret_cast<T *>(shared);
		assert (cuda::is_aligned<16>(shared_conf_energies));

		// init the sizes for cudaMalloc'd Array instances
		out_coords->init(conf_space->max_num_conf_atoms);

		// make the atoms
		osprey::Assignment<T> assignment(*conf_space, *conf, *out_coords, shared_atom_pairs, shared_conf_energies);

		// call the energy function
		*out_energy = efunc(assignment, *inters, thread_energy);
	}

	template<typename T, osprey::EnergyFunction<T> efunc>
	static T calc(int device,
	              Stream * stream,
	              const osprey::ConfSpace<T> * d_conf_space,
	              const osprey::Array<int32_t> & conf,
	              const osprey::Array<osprey::PosInter<T>> & inters,
	              osprey::Array<osprey::Real3<T>> * out_coords,
	              int64_t num_atoms) {

		CUDACHECK(cudaSetDevice(device));

		// upload the arguments
		size_t conf_bytes = conf.get_bytes();
		osprey::Array<int32_t> * d_conf;
		CUDACHECK(cudaMalloc(&d_conf, conf_bytes));
		CUDACHECK(cudaMemcpy(d_conf, &conf, conf_bytes, cudaMemcpyHostToDevice));

		size_t inters_bytes = inters.get_bytes();
		osprey::Array<osprey::PosInter<T>> * d_inters;
		CUDACHECK(cudaMalloc(&d_inters, inters_bytes));
		CUDACHECK(cudaMemcpy(d_inters, &inters, inters_bytes, cudaMemcpyHostToDevice));

		// TODO: combine allocations/transfers for speed?

		// allocate space for the coords
		size_t out_coords_bytes = osprey::Array<osprey::Real3<T>>::get_bytes(num_atoms);
		osprey::Array<osprey::Real3<T>> * d_out_coords;
		CUDACHECK(cudaMalloc(&d_out_coords, out_coords_bytes));

		// allocate space for the energy
		size_t out_energy_bytes = sizeof(T);
		T * d_out_energy;
		CUDACHECK(cudaMalloc(&d_out_energy, out_energy_bytes));

		// compute the shared memory size
		int64_t shared_size_static = 0
			+ cuda::pad_to_alignment<16>(Assignment<T>::sizeof_atom_pairs(conf.get_size()))
			+ cuda::pad_to_alignment<16>(Assignment<T>::sizeof_conf_energies(conf.get_size()));
		int64_t shared_size_per_thread = sizeof(T);

		// launch the kernel
		int num_threads = CALC_THREADS;
		int64_t shared_size = shared_size_static + shared_size_per_thread*num_threads;
		calc_kernel<T,efunc><<<1, num_threads, shared_size, stream->get_stream()>>>(d_conf_space, d_conf, d_inters, d_out_coords, d_out_energy);
		cuda::check_error();

		// download the energy
		T out_energy;
		CUDACHECK(cudaMemcpy(&out_energy, d_out_energy, out_energy_bytes, cudaMemcpyDeviceToHost));

		// download the coords, if needed
		if (out_coords != nullptr) {
			CUDACHECK(cudaMemcpy(out_coords, d_out_coords, out_coords_bytes, cudaMemcpyDeviceToHost));
		}

		// cleanup
		CUDACHECK(cudaFree(d_conf));
		CUDACHECK(cudaFree(d_inters));
		CUDACHECK(cudaFree(d_out_coords));
		CUDACHECK(cudaFree(d_out_energy));

		return out_energy;
	}
}

API float32_t calc_amber_eef1_f32(int device,
                                  Stream * stream,
                                  const osprey::ConfSpace<float32_t> * d_conf_space,
                                  const osprey::Array<int32_t> & conf,
                                  const osprey::Array<osprey::PosInter<float32_t>> & inters,
                                  osprey::Array<osprey::Real3<float32_t>> * out_coords,
                                  int64_t num_atoms) {
	return osprey::calc<float32_t,osprey::ambereef1::calc_energy>(device, stream, d_conf_space, conf, inters, out_coords, num_atoms);
}
API float64_t calc_amber_eef1_f64(int device,
                                  Stream * stream,
                                  const osprey::ConfSpace<float64_t> * d_conf_space,
                                  const osprey::Array<int32_t> & conf,
                                  const osprey::Array<osprey::PosInter<float64_t>> & inters,
                                  osprey::Array<osprey::Real3<float64_t>> * out_coords,
                                  int64_t num_atoms) {
	return osprey::calc<float64_t,osprey::ambereef1::calc_energy>(device, stream, d_conf_space, conf, inters, out_coords, num_atoms);
}


namespace osprey {

	// configure the launch bounds, so the compiler knows how many registers to use
	#define MINIMIZE_THREADS_F32_VOLTA 1024
	#define MINIMIZE_THREADS_F64_VOLTA 1024
	#define MINIMIZE_THREADS_F32_PASCAL 1024 // TODO: optimize for pascal?
	#define MINIMIZE_THREADS_F64_PASCAL 1024
	#define MINIMIZE_THREADS_F32_MAXWELL 1024 // TODO: optimize for maxwell?
	#define MINIMIZE_THREADS_F64_MAXWELL 1024

	// the kernels use __syncthreads a lot, so they're designed to take up an entire SM
	#define MINIMIZE_BLOCKS 1

	#if __CUDA_ARCH__ >= 700 // volta
		#define MINIMIZE_THREADS_F32 MINIMIZE_THREADS_F32_VOLTA
		#define MINIMIZE_THREADS_F64 MINIMIZE_THREADS_F64_VOLTA
	#elif __CUDA_ARCH__ >= 600 // pascal
		#define MINIMIZE_THREADS_F32 MINIMIZE_THREADS_F32_PASCAL
		#define MINIMIZE_THREADS_F64 MINIMIZE_THREADS_F64_PASCAL
	#elif __CUDA_ARCH__ >= 500 // maxwell
		#define MINIMIZE_THREADS_F32 MINIMIZE_THREADS_F32_MAXWELL
		#define MINIMIZE_THREADS_F64 MINIMIZE_THREADS_F64_MAXWELL
	#else
		// host side
		#define MINIMIZE_THREADS_F32 -1
		#define MINIMIZE_THREADS_F64 -1
	#endif

	template<typename T>
	__host__
	inline int get_minimize_threads();

	template<>
	__host__
	inline int get_minimize_threads<float32_t>() {
		int arch = cuda::get_arch();
		if (arch >= 700) {
			return MINIMIZE_THREADS_F32_VOLTA;
		} else if (arch >= 600) {
			return MINIMIZE_THREADS_F32_PASCAL;
		} else if (arch >= 500) {
			return MINIMIZE_THREADS_F32_MAXWELL;
		} else {
			throw std::runtime_error("unsupported CUDA architecture");
		}
	}

	template<>
	__host__
	inline int get_minimize_threads<float64_t>() {
		int arch = cuda::get_arch();
		if (arch >= 700) {
			return MINIMIZE_THREADS_F64_VOLTA;
		} else if (arch >= 600) {
			return MINIMIZE_THREADS_F64_PASCAL;
		} else if (arch >= 500) {
			return MINIMIZE_THREADS_F64_MAXWELL;
		} else {
			throw std::runtime_error("unsupported CUDA architecture");
		}
	}

	template<typename T, EnergyFunction<T> efunc>
	__device__
	void minimize_kernel_impl(int64_t shared_size,
	                          const ConfSpace<T> * conf_space,
	                          const Array<int32_t> * conf,
	                          const Array<PosInter<T>> * inters,
	                          Array<Real3<T>> * out_coords,
	                          DofValues<T> * out_dof_values) {

		// slice up the shared memory
		extern __shared__ int8_t shared[];
		int64_t offset = 0;

		auto shared_atom_pairs = reinterpret_cast<const void **>(shared + offset);
		assert (cuda::is_aligned<16>(shared_atom_pairs));
		offset += cuda::pad_to_alignment<16>(Assignment<T>::sizeof_atom_pairs(conf_space->num_pos));

		auto shared_conf_energies = reinterpret_cast<T *>(shared + offset);
		assert (cuda::is_aligned<16>(shared_conf_energies));
		offset += cuda::pad_to_alignment<16>(Assignment<T>::sizeof_conf_energies(conf_space->num_pos));

		auto shared_dofs_mem = reinterpret_cast<int8_t *>(shared + offset);
		assert (cuda::is_aligned<16>(shared_dofs_mem));
		offset += Dofs<T>::dofs_shared_size;

		auto shared_dof_bufs = reinterpret_cast<int8_t *>(shared + offset);
		assert (cuda::is_aligned<16>(shared_dof_bufs));
		offset += (Dof<T>::buf_size + Array<const PosInter<T> *>::get_bytes(inters->get_size()))*conf_space->max_num_dofs;

		auto shared_line_search_states = reinterpret_cast<LineSearchState<T> *>(shared + offset);
		assert (cuda::is_aligned<16>(shared_line_search_states));
		offset += cuda::pad_to_alignment<16>(sizeof(LineSearchState<T>)*conf_space->max_num_dofs);

		auto shared_here = reinterpret_cast<DofValues<T> *>(shared + offset);
		assert (cuda::is_aligned<16>(shared_here));
		offset += DofValues<T>::get_bytes(conf_space->max_num_dofs);

		auto shared_next = reinterpret_cast<DofValues<T> *>(shared + offset);
		assert (cuda::is_aligned<16>(shared_next));
		offset += DofValues<T>::get_bytes(conf_space->max_num_dofs);

		auto thread_energy = reinterpret_cast<T *>(shared + offset);
		assert (cuda::is_aligned<16>(thread_energy));
		offset += sizeof(T)*blockDim.x;

		assert (offset == shared_size);

		// init the sizes for cudaMalloc'd Array instances
		out_coords->init(conf_space->max_num_conf_atoms);
		out_dof_values->x.init(conf_space->max_num_dofs);

		// init sizes for Array instances in shared memory
		shared_here->x.init(conf_space->max_num_dofs);
		shared_next->x.init(conf_space->max_num_dofs);

		// make the atoms and the degrees of freedom
		Assignment<T> assignment(*conf_space, *conf, *out_coords, shared_atom_pairs, shared_conf_energies);
		Dofs<T> dofs(assignment, *inters, efunc, thread_energy, shared_dof_bufs, shared_dofs_mem);

		// truncate the dof values to match the actual number of dofs
		if (threadIdx.x == 0) {
			out_dof_values->x.truncate(dofs.get_size());
			shared_here->x.truncate(dofs.get_size());
			shared_next->x.truncate(dofs.get_size());
		}
		__syncthreads();

		// init the dofs to the center of the voxel
		for (int d=threadIdx.x; d<dofs.get_size(); d+=blockDim.x) {
			shared_here->x[d] = dofs[d].center();
		}
		__syncthreads();

		// use CCD to do the minimization
		minimize_ccd(dofs, *shared_here, *shared_next, shared_line_search_states);

		// copy dof values to global memory
		out_dof_values->set(*shared_here);
	}

	// set kernel launch bounds for each real type

	template<typename T, EnergyFunction<T> efunc>
	__global__
	void minimize_kernel(int64_t shared_size,
	                     const ConfSpace<T> * conf_space,
	                     const Array<int32_t> * conf,
	                     const Array<PosInter<T>> * inters,
	                     Array<Real3<T>> * out_coords,
	                     DofValues<T> * out_dof_values);

	template<>
	__global__
	__launch_bounds__(MINIMIZE_THREADS_F32, MINIMIZE_BLOCKS)
	void minimize_kernel<float32_t,osprey::ambereef1::calc_energy>(int64_t shared_size,
	                                                               const ConfSpace<float32_t> * conf_space,
	                                                               const Array<int32_t> * conf,
	                                                               const Array<PosInter<float32_t>> * inters,
	                                                               Array<Real3<float32_t>> * out_coords,
	                                                               DofValues<float32_t> * out_dof_values) {
		minimize_kernel_impl<float32_t,osprey::ambereef1::calc_energy>(shared_size, conf_space, conf, inters, out_coords, out_dof_values);
	}

	template<>
	__global__
	__launch_bounds__(MINIMIZE_THREADS_F64, MINIMIZE_BLOCKS)
	void minimize_kernel<float64_t,osprey::ambereef1::calc_energy>(int64_t shared_size,
	                                                               const ConfSpace<float64_t> * conf_space,
	                                                               const Array<int32_t> * conf,
	                                                               const Array<PosInter<float64_t>> * inters,
	                                                               Array<Real3<float64_t>> * out_coords,
	                                                               DofValues<float64_t> * out_dof_values) {
		minimize_kernel_impl<float64_t,osprey::ambereef1::calc_energy>(shared_size, conf_space, conf, inters, out_coords, out_dof_values);
	}

	template<typename T, EnergyFunction<T> efunc>
	static T minimize(int device,
	                  Stream * stream,
	                  const ConfSpace<T> * d_conf_space,
	                  const Array<int32_t> & conf,
	                  const Array<PosInter<T>> & inters,
	                  Array<Real3<T>> * out_coords,
	                  int64_t num_atoms,
	                  Array<T> * out_dofs,
	                  int64_t max_num_dofs) {

		CUDACHECK(cudaSetDevice(device));

		// split up the transfer buffer
		size_t conf_bytes = conf.get_bytes();
		size_t inters_bytes = inters.get_bytes();
		size_t out_coords_bytes = Array<Real3<T>>::get_bytes(num_atoms);
		size_t out_dof_values_bytes = DofValues<T>::get_bytes(max_num_dofs);
		if (conf_bytes + inters_bytes + out_coords_bytes + out_dof_values_bytes > Stream::d_buf_size) {
			throw std::runtime_error("exceeded transfer buffer size");
		}

		int8_t * d_xfer_buf = stream->get_d_buf();

		auto d_conf = reinterpret_cast<Array<int32_t> *>(d_xfer_buf);
		assert (cuda::is_aligned<16>(d_conf));
		d_xfer_buf += conf_bytes;

		auto d_inters = reinterpret_cast<Array<PosInter<T>> *>(d_xfer_buf);
		assert (cuda::is_aligned<16>(d_inters));
		d_xfer_buf += inters_bytes;

		auto d_out_coords = reinterpret_cast<Array<Real3<T>> *>(d_xfer_buf);
		assert (cuda::is_aligned<16>(d_out_coords));
		d_xfer_buf += out_coords_bytes;

		auto d_out_dof_values = reinterpret_cast<DofValues<T> *>(d_xfer_buf);
		assert (cuda::is_aligned<16>(d_out_dof_values));

		// upload the arguments
		CUDACHECK(cudaMemcpyAsync(d_conf, &conf, conf_bytes, cudaMemcpyHostToDevice, stream->get_stream()));
		CUDACHECK(cudaMemcpyAsync(d_inters, &inters, inters_bytes, cudaMemcpyHostToDevice, stream->get_stream()));

		// compute the shared memory size
		int64_t shared_size_static = 0
			+ cuda::pad_to_alignment<16>(Assignment<T>::sizeof_atom_pairs(conf.get_size()))
			+ cuda::pad_to_alignment<16>(Assignment<T>::sizeof_conf_energies(conf.get_size()))
			+ (Dof<T>::buf_size + Array<const PosInter<T> *>::get_bytes(inters.get_size()))*max_num_dofs
			+ Dofs<T>::dofs_shared_size
			+ cuda::pad_to_alignment<16>(sizeof(LineSearchState<T>)*max_num_dofs) // line_search_states
			+ DofValues<T>::get_bytes(max_num_dofs) // here
			+ DofValues<T>::get_bytes(max_num_dofs); // next
		int64_t shared_size_per_thread = sizeof(T);

		// launch the kernel
		int num_threads = get_minimize_threads<T>();
		int64_t shared_size = shared_size_static + shared_size_per_thread*num_threads;
		minimize_kernel<T,efunc><<<1, num_threads, shared_size, stream->get_stream()>>>(shared_size, d_conf_space, d_conf, d_inters, d_out_coords, d_out_dof_values);
		cuda::check_error();

		// download the energy
		CUDACHECK(cudaMemcpyAsync(stream->get_h_buf(), DofValues<T>::f_ptr(d_out_dof_values), sizeof(T), cudaMemcpyDeviceToHost, stream->get_stream()));

		// download the coords, if needed
		if (out_coords != nullptr) {
			CUDACHECK(cudaMemcpyAsync(out_coords, d_out_coords, out_coords_bytes, cudaMemcpyDeviceToHost, stream->get_stream()));
		}

		// download the dofs, if needed
		if (out_dofs != nullptr) {
			// TODO: does this work?
			CUDACHECK(cudaMemcpyAsync(out_dofs, DofValues<T>::x_ptr(d_out_dof_values), Array<T>::get_bytes(max_num_dofs), cudaMemcpyDeviceToHost, stream->get_stream()));
		}

		// wait for the downloads to finish
		CUDACHECK(cudaStreamSynchronize(stream->get_stream()));

		// finally, return the minimized energy
		const T * out_energy = reinterpret_cast<const T *>(stream->get_h_buf());
		return *out_energy;
	}
}

API float32_t minimize_amber_eef1_f32(int device,
                                      Stream * stream,
                                      const osprey::ConfSpace<float32_t> * d_conf_space,
                                      const osprey::Array<int32_t> & conf,
                                      const osprey::Array<osprey::PosInter<float32_t>> & inters,
                                      osprey::Array<osprey::Real3<float32_t>> * out_coords,
                                      int64_t num_atoms,
                                      osprey::Array<float32_t> * out_dofs,
                                      int64_t max_num_dofs) {
	return osprey::minimize<float32_t,osprey::ambereef1::calc_energy>(device, stream, d_conf_space, conf, inters, out_coords, num_atoms, out_dofs, max_num_dofs);
}
API float64_t minimize_amber_eef1_f64(int device,
                                      Stream * stream,
                                      const osprey::ConfSpace<float64_t> * d_conf_space,
                                      const osprey::Array<int32_t> & conf,
                                      const osprey::Array<osprey::PosInter<float64_t>> & inters,
                                      osprey::Array<osprey::Real3<float64_t>> * out_coords,
                                      int64_t num_atoms,
                                      osprey::Array<float64_t> * out_dofs,
                                      int64_t max_num_dofs) {
	return osprey::minimize<float64_t,osprey::ambereef1::calc_energy>(device, stream, d_conf_space, conf, inters, out_coords, num_atoms, out_dofs, max_num_dofs);
}
