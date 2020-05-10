
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
	// 4 bytes pad
	int64_t mem_total;
	int64_t mem_free;
};
ASSERT_JAVA_COMPATIBLE(GpuInfo, 304);

API osprey::Array<GpuInfo> * alloc_gpu_infos() {

	int num_gpus;
	cudaGetDeviceCount(&num_gpus);
	cuda::check_error();

	osprey::Array<GpuInfo> * infos = osprey::Array<GpuInfo>::make(num_gpus);
	for (int device=0; device<num_gpus; device++) {
		GpuInfo & info = (*infos)[device];

		cudaDeviceGetPCIBusId(reinterpret_cast<char *>(&info.bus_id), 16, device);
		cuda::check_error();

		cudaDeviceProp props {};
		cudaGetDeviceProperties(&props, device);
		cuda::check_error();

		memcpy(&info.name, props.name, 256);
		info.integrated = props.integrated;
		info.concurrent_kernels = props.concurrentKernels;
		info.num_processors = props.multiProcessorCount;

		size_t mem_free;
		size_t mem_total;
		cudaMemGetInfo(&mem_free, &mem_total);
		cuda::check_error();

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

		cudaSetDevice(device);
		cuda::check_error();

		// allocate the device memory
		void * p_device;
		cudaMalloc(&p_device, conf_space.size);
		cuda::check_error();

		// upload the conf space
		cudaMemcpy(p_device, &conf_space, conf_space.size, cudaMemcpyHostToDevice);
		cuda::check_error();

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

	cudaSetDevice(device);
	cuda::check_error();

	cudaFree(p);
	cuda::check_error();
}

API cudaStream_t alloc_stream(int device) {

	cudaSetDevice(device);
	cuda::check_error();

	cudaStream_t stream;
	cudaStreamCreate(&stream);
	cuda::check_error();
	return stream;
}

API void free_stream(int device, cudaStream_t stream) {

	cudaSetDevice(device);
	cuda::check_error();

	cudaStreamDestroy(stream);
	cuda::check_error();
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
	                   cudaStream_t stream,
	                   const ConfSpace<T> * d_conf_space,
	                   const Array<int32_t> & conf,
	                   Array<Real3<T>> & out_coords) {

		cudaSetDevice(device);
		cuda::check_error();

		// upload the arguments
		size_t conf_bytes = conf.get_bytes();
		Array<int32_t> * d_conf;
		cudaMalloc(&d_conf, conf_bytes);
		cuda::check_error();
		cudaMemcpy(d_conf, &conf, conf_bytes, cudaMemcpyHostToDevice);
		cuda::check_error();

		// allocate space for the coords
		size_t out_coords_bytes = out_coords.get_bytes();
		Array<Real3<T>> * d_out_coords;
		cudaMalloc(&d_out_coords, out_coords_bytes);
		cuda::check_error();

		// compute the shared memory size
		int64_t shared_size = 0
			+ Assignment<T>::sizeof_atom_pairs(conf.get_size())
			+ Assignment<T>::sizeof_conf_energies(conf.get_size());

		// launch the kernel
		// TODO: optimize launch bounds
		int num_threads = cuda::optimize_threads(assign_kernel<T>, shared_size, 0);
		assign_kernel<<<1, num_threads, shared_size, stream>>>(d_conf_space, d_conf, d_out_coords);
		cuda::check_error();

		// download the coords
		cudaMemcpy(&out_coords, d_out_coords, out_coords_bytes, cudaMemcpyDeviceToHost);
		cuda::check_error();

		// cleanup
		cudaFree(d_conf);
		cudaFree(d_out_coords);
	}
}

API void assign_f32(int device,
                    cudaStream_t stream,
                    const osprey::ConfSpace<float32_t> * d_conf_space,
                    const osprey::Array<int32_t> & conf,
                    osprey::Array<osprey::Real3<float32_t>> & out_coords) {
	osprey::assign(device, stream, d_conf_space, conf, out_coords);
}
API void assign_f64(int device,
                    cudaStream_t stream,
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
	              cudaStream_t stream,
	              const osprey::ConfSpace<T> * d_conf_space,
	              const osprey::Array<int32_t> & conf,
	              const osprey::Array<osprey::PosInter<T>> & inters,
	              osprey::Array<osprey::Real3<T>> * out_coords,
	              int64_t num_atoms) {

		cudaSetDevice(device);
		cuda::check_error();

		// upload the arguments
		size_t conf_bytes = conf.get_bytes();
		osprey::Array<int32_t> * d_conf;
		cudaMalloc(&d_conf, conf_bytes);
		cuda::check_error();
		cudaMemcpy(d_conf, &conf, conf_bytes, cudaMemcpyHostToDevice);
		cuda::check_error();

		size_t inters_bytes = inters.get_bytes();
		osprey::Array<osprey::PosInter<T>> * d_inters;
		cudaMalloc(&d_inters, inters_bytes);
		cuda::check_error();
		cudaMemcpy(d_inters, &inters, inters_bytes, cudaMemcpyHostToDevice);
		cuda::check_error();

		// TODO: combine allocations/transfers for speed?

		// allocate space for the coords
		size_t out_coords_bytes = osprey::Array<osprey::Real3<T>>::get_bytes(num_atoms);
		osprey::Array<osprey::Real3<T>> * d_out_coords;
		cudaMalloc(&d_out_coords, out_coords_bytes);
		cuda::check_error();

		// allocate space for the energy
		size_t out_energy_bytes = sizeof(T);
		T * d_out_energy;
		cudaMalloc(&d_out_energy, out_energy_bytes);
		cuda::check_error();

		// compute the shared memory size
		int64_t shared_size_static = 0
			+ cuda::pad_to_alignment<16>(Assignment<T>::sizeof_atom_pairs(conf.get_size()))
			+ cuda::pad_to_alignment<16>(Assignment<T>::sizeof_conf_energies(conf.get_size()));
		int64_t shared_size_per_thread = sizeof(T);

		// launch the kernel
		int num_threads = CALC_THREADS;
		int64_t shared_size = shared_size_static + shared_size_per_thread*num_threads;
		calc_kernel<T,efunc><<<1, num_threads, shared_size, stream>>>(d_conf_space, d_conf, d_inters, d_out_coords, d_out_energy);
		cuda::check_error();

		// download the energy
		T out_energy;
		cudaMemcpy(&out_energy, d_out_energy, out_energy_bytes, cudaMemcpyDeviceToHost);
		cuda::check_error();

		// download the coords, if needed
		if (out_coords != nullptr) {
			cudaMemcpy(out_coords, d_out_coords, out_coords_bytes, cudaMemcpyDeviceToHost);
			cuda::check_error();
		}

		// cleanup
		cudaFree(d_conf);
		cudaFree(d_inters);
		cudaFree(d_out_coords);
		cudaFree(d_out_energy);

		return out_energy;
	}
}

API float32_t calc_amber_eef1_f32(int device,
                                  cudaStream_t stream,
                                  const osprey::ConfSpace<float32_t> * d_conf_space,
                                  const osprey::Array<int32_t> & conf,
                                  const osprey::Array<osprey::PosInter<float32_t>> & inters,
                                  osprey::Array<osprey::Real3<float32_t>> * out_coords,
                                  int64_t num_atoms) {
	return osprey::calc<float32_t,osprey::ambereef1::calc_energy>(device, stream, d_conf_space, conf, inters, out_coords, num_atoms);
}
API float64_t calc_amber_eef1_f64(int device,
                                  cudaStream_t stream,
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
	void minimize_kernel_impl(const ConfSpace<T> * conf_space,
	                          const Array<int32_t> * conf,
	                          const Array<PosInter<T>> * inters,
	                          Array<Real3<T>> * out_coords,
	                          DofValues<T> * out_dof_values) {

		// slice up the shared memory
		extern __shared__ int8_t _shared[];
		auto shared = _shared;

		auto shared_atom_pairs = reinterpret_cast<const void **>(shared);
		assert (cuda::is_aligned<16>(shared_atom_pairs));
		shared += cuda::pad_to_alignment<16>(Assignment<T>::sizeof_atom_pairs(conf_space->num_pos));

		auto shared_conf_energies = reinterpret_cast<T *>(shared);
		assert (cuda::is_aligned<16>(shared_conf_energies));
		shared += cuda::pad_to_alignment<16>(Assignment<T>::sizeof_conf_energies(conf_space->num_pos));

		auto shared_dofs_mem = reinterpret_cast<int8_t *>(shared);
		assert (cuda::is_aligned<16>(shared_dofs_mem));
		shared += Dofs<T>::dofs_shared_size;

		auto shared_dofs = reinterpret_cast<Dof<T> **>(shared);
		assert (cuda::is_aligned<16>(shared_dofs));
		shared += cuda::pad_to_alignment<16>(sizeof(Dof<T> *)*conf_space->max_num_dofs);

		auto shared_line_search_states = reinterpret_cast<LineSearchState<T> *>(shared);
		assert (cuda::is_aligned<16>(shared_line_search_states));
		shared += cuda::pad_to_alignment<16>(sizeof(LineSearchState<T>)*conf_space->max_num_dofs);

		auto shared_here = reinterpret_cast<DofValues<T> *>(shared);
		assert (cuda::is_aligned<16>(shared_here));
		shared += DofValues<T>::get_bytes(conf_space->max_num_dofs);

		auto shared_next = reinterpret_cast<DofValues<T> *>(shared);
		assert (cuda::is_aligned<16>(shared_next));
		shared += DofValues<T>::get_bytes(conf_space->max_num_dofs);

		auto thread_energy = reinterpret_cast<T *>(shared);
		assert (cuda::is_aligned<16>(thread_energy));

		// init the sizes for cudaMalloc'd Array instances
		out_coords->init(conf_space->max_num_conf_atoms);
		out_dof_values->x.init(conf_space->max_num_dofs);

		// init sizes for Array instances in shared memory
		shared_here->x.init(conf_space->max_num_dofs);
		shared_next->x.init(conf_space->max_num_dofs);

		// make the atoms and the degrees of freedom
		Assignment<T> assignment(*conf_space, *conf, *out_coords, shared_atom_pairs, shared_conf_energies);
		Dofs<T> dofs(assignment, *inters, efunc, thread_energy, shared_dofs_mem, shared_dofs);

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

		// cleanup
		dofs.free();

		// copy dof values to global memory
		out_dof_values->set(*shared_here);
	}

	// set kernel launch bounds for each real type

	template<typename T, EnergyFunction<T> efunc>
	__global__
	void minimize_kernel(const ConfSpace<T> * conf_space,
	                     const Array<int32_t> * conf,
	                     const Array<PosInter<T>> * inters,
	                     Array<Real3<T>> * out_coords,
	                     DofValues<T> * out_dof_values);

	template<>
	__global__
	__launch_bounds__(MINIMIZE_THREADS_F32, MINIMIZE_BLOCKS)
	void minimize_kernel<float32_t,osprey::ambereef1::calc_energy>(const ConfSpace<float32_t> * conf_space,
	                                                               const Array<int32_t> * conf,
	                                                               const Array<PosInter<float32_t>> * inters,
	                                                               Array<Real3<float32_t>> * out_coords,
	                                                               DofValues<float32_t> * out_dof_values) {
		minimize_kernel_impl<float32_t,osprey::ambereef1::calc_energy>(conf_space, conf, inters, out_coords, out_dof_values);
	}

	template<>
	__global__
	__launch_bounds__(MINIMIZE_THREADS_F64, MINIMIZE_BLOCKS)
	void minimize_kernel<float64_t,osprey::ambereef1::calc_energy>(const ConfSpace<float64_t> * conf_space,
	                                                               const Array<int32_t> * conf,
	                                                               const Array<PosInter<float64_t>> * inters,
	                                                               Array<Real3<float64_t>> * out_coords,
	                                                               DofValues<float64_t> * out_dof_values) {
		minimize_kernel_impl<float64_t,osprey::ambereef1::calc_energy>(conf_space, conf, inters, out_coords, out_dof_values);
	}

	template<typename T, EnergyFunction<T> efunc>
	static T minimize(int device,
	                  cudaStream_t stream,
	                  const ConfSpace<T> * d_conf_space,
	                  const Array<int32_t> & conf,
	                  const Array<PosInter<T>> & inters,
	                  Array<Real3<T>> * out_coords,
	                  int64_t num_atoms,
	                  Array<T> * out_dofs,
	                  int64_t num_dofs) {

		cudaSetDevice(device);
		cuda::check_error();

		// upload the arguments
		size_t conf_bytes = conf.get_bytes();
		Array<int32_t> * d_conf;
		cudaMalloc(&d_conf, conf_bytes);
		cuda::check_error();
		cudaMemcpy(d_conf, &conf, conf_bytes, cudaMemcpyHostToDevice);
		cuda::check_error();

		size_t inters_bytes = inters.get_bytes();
		Array<PosInter<T>> * d_inters;
		cudaMalloc(&d_inters, inters_bytes);
		cuda::check_error();
		cudaMemcpy(d_inters, &inters, inters_bytes, cudaMemcpyHostToDevice);
		cuda::check_error();

		// TODO: combine allocations/transfers for speed?

		// allocate space for the coords
		size_t out_coords_bytes = Array<Real3<T>>::get_bytes(num_atoms);
		Array<Real3<T>> * d_out_coords;
		cudaMalloc(&d_out_coords, out_coords_bytes);
		cuda::check_error();

		// allocte space for the dofs
		size_t out_dof_values_bytes = DofValues<T>::get_bytes(num_dofs);
		DofValues<T> * d_out_dof_values;
		cudaMalloc(&d_out_dof_values, out_dof_values_bytes);
		cuda::check_error();

		// compute the shared memory size
		int64_t shared_size_static = 0
			+ cuda::pad_to_alignment<16>(Assignment<T>::sizeof_atom_pairs(conf.get_size()))
			+ cuda::pad_to_alignment<16>(Assignment<T>::sizeof_conf_energies(conf.get_size()))
			+ cuda::pad_to_alignment<16>(sizeof(Dof<T> *)*num_dofs)
			+ Dofs<T>::dofs_shared_size
			+ cuda::pad_to_alignment<16>(sizeof(LineSearchState<T>)*num_dofs) // line_search_states
			+ DofValues<T>::get_bytes(num_dofs) // here
			+ DofValues<T>::get_bytes(num_dofs); // next
		int64_t shared_size_per_thread = sizeof(T);

		// launch the kernel
		int num_threads = get_minimize_threads<T>();
		int64_t shared_size = shared_size_static + shared_size_per_thread*num_threads;
		minimize_kernel<T,efunc><<<1, num_threads, shared_size, stream>>>(d_conf_space, d_conf, d_inters, d_out_coords, d_out_dof_values);
		cuda::check_error();

		// download the energy
		T out_energy;
		cudaMemcpy(&out_energy, DofValues<T>::f_ptr(d_out_dof_values), sizeof(T), cudaMemcpyDeviceToHost);
		cuda::check_error();

		// download the coords, if needed
		if (out_coords != nullptr) {
			cudaMemcpy(out_coords, d_out_coords, out_coords_bytes, cudaMemcpyDeviceToHost);
			cuda::check_error();
		}

		// download the dofs, if needed
		if (out_dofs != nullptr) {
			// TODO: does this work?
			cudaMemcpy(out_dofs, DofValues<T>::x_ptr(d_out_dof_values), Array<T>::get_bytes(num_dofs), cudaMemcpyDeviceToHost);
			cuda::check_error();
		}

		// cleanup
		cudaFree(d_conf);
		cudaFree(d_inters);
		cudaFree(d_out_coords);
		cudaFree(d_out_dof_values);

		return out_energy;
	}
}

API float32_t minimize_amber_eef1_f32(int device,
                                      cudaStream_t stream,
                                      const osprey::ConfSpace<float32_t> * d_conf_space,
                                      const osprey::Array<int32_t> & conf,
                                      const osprey::Array<osprey::PosInter<float32_t>> & inters,
                                      osprey::Array<osprey::Real3<float32_t>> * out_coords,
                                      int64_t num_atoms,
                                      osprey::Array<float32_t> * out_dofs,
                                      int64_t num_dofs) {
	return osprey::minimize<float32_t,osprey::ambereef1::calc_energy>(device, stream, d_conf_space, conf, inters, out_coords, num_atoms, out_dofs, num_dofs);
}
API float64_t minimize_amber_eef1_f64(int device,
                                      cudaStream_t stream,
                                      const osprey::ConfSpace<float64_t> * d_conf_space,
                                      const osprey::Array<int32_t> & conf,
                                      const osprey::Array<osprey::PosInter<float64_t>> & inters,
                                      osprey::Array<osprey::Real3<float64_t>> * out_coords,
                                      int64_t num_atoms,
                                      osprey::Array<float64_t> * out_dofs,
                                      int64_t num_dofs) {
	return osprey::minimize<float64_t,osprey::ambereef1::calc_energy>(device, stream, d_conf_space, conf, inters, out_coords, num_atoms, out_dofs, num_dofs);
}
