
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

namespace osprey {

	template<typename T>
	static void * alloc_conf_space(const ConfSpace<T> & conf_space) {

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

API void * alloc_conf_space_f32(const osprey::ConfSpace<float32_t> & conf_space) {
	return osprey::alloc_conf_space(conf_space);
}
API void * alloc_conf_space_f64(const osprey::ConfSpace<float64_t> & conf_space) {
	return osprey::alloc_conf_space(conf_space);
}

API void free_conf_space(void * p) {
	cudaFree(p);
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
	static void assign(const ConfSpace<T> * d_conf_space,
	                   const Array<int32_t> & conf,
	                   Array<Real3<T>> & out_coords) {

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
		int num_threads = cuda::optimize_threads(assign_kernel<T>, shared_size, 0);
		assign_kernel<<<1, num_threads, shared_size>>>(d_conf_space, d_conf, d_out_coords);
		cuda::check_error();

		// download the coords
		cudaMemcpy(&out_coords, d_out_coords, out_coords_bytes, cudaMemcpyDeviceToHost);
		cuda::check_error();

		// cleanup
		cudaFree(d_conf);
		cudaFree(d_out_coords);
	}
}

API void assign_f32(const osprey::ConfSpace<float32_t> * d_conf_space,
                    const osprey::Array<int32_t> & conf,
                    osprey::Array<osprey::Real3<float32_t>> & out_coords) {
	osprey::assign(d_conf_space, conf, out_coords);
}
API void assign_f64(const osprey::ConfSpace<float64_t> * d_conf_space,
                    const osprey::Array<int32_t> & conf,
                    osprey::Array<osprey::Real3<float64_t>> & out_coords) {
	osprey::assign(d_conf_space, conf, out_coords);
}


namespace osprey {

	template<typename T, osprey::EnergyFunction<T> efunc>
	__global__
	void calc_kernel(const osprey::ConfSpace<T> * conf_space,
	                 const osprey::Array<int32_t> * conf,
	                 const osprey::Array<osprey::PosInter<T>> * inters,
	                 osprey::Array<osprey::Real3<T>> * out_coords,
	                 T * out_energy) {

		// slice up the shared memory
		extern __shared__ int8_t shared[];
		int64_t shared_offset = 0;
		auto shared_atom_pairs = reinterpret_cast<const void **>(shared + shared_offset);
		shared_offset += osprey::Assignment<T>::sizeof_atom_pairs(conf_space->num_pos);
		auto shared_conf_energies = reinterpret_cast<T *>(shared + shared_offset);
		shared_offset += osprey::Assignment<T>::sizeof_conf_energies(conf_space->num_pos);
		auto thread_energy = reinterpret_cast<T *>(shared + shared_offset);

		// init the sizes for cudaMalloc'd Array instances
		out_coords->init(conf_space->max_num_conf_atoms);

		// make the atoms
		osprey::Assignment<T> assignment(*conf_space, *conf, *out_coords, shared_atom_pairs, shared_conf_energies);

		// call the energy function
		*out_energy = efunc(assignment, *inters, thread_energy);
	}

	template<typename T, osprey::EnergyFunction<T> efunc>
	static T calc(const osprey::ConfSpace<T> * d_conf_space,
	              const osprey::Array<int32_t> & conf,
	              const osprey::Array<osprey::PosInter<T>> & inters,
	              osprey::Array<osprey::Real3<T>> * out_coords,
	              int64_t num_atoms) {

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
			+ osprey::Assignment<T>::sizeof_atom_pairs(conf.get_size())
			+ osprey::Assignment<T>::sizeof_conf_energies(conf.get_size());
		int64_t shared_size_per_thread = sizeof(T);

		// launch the kernel
		// TODO: cache thread optimization
		auto kernel = calc_kernel<T,efunc>;
		int num_threads = cuda::optimize_threads(*kernel, shared_size_static, shared_size_per_thread);
		int64_t shared_size = shared_size_static + shared_size_per_thread*num_threads;
		kernel<<<1, num_threads, shared_size>>>(d_conf_space, d_conf, d_inters, d_out_coords, d_out_energy);
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

API float32_t calc_amber_eef1_f32(const osprey::ConfSpace<float32_t> * d_conf_space,
                                  const osprey::Array<int32_t> & conf,
                                  const osprey::Array<osprey::PosInter<float32_t>> & inters,
                                  osprey::Array<osprey::Real3<float32_t>> * out_coords,
                                  int64_t num_atoms) {
	return osprey::calc<float32_t,osprey::ambereef1::calc_energy>(d_conf_space, conf, inters, out_coords, num_atoms);
}
API float64_t calc_amber_eef1_f64(const osprey::ConfSpace<float64_t> * d_conf_space,
                                  const osprey::Array<int32_t> & conf,
                                  const osprey::Array<osprey::PosInter<float64_t>> & inters,
                                  osprey::Array<osprey::Real3<float64_t>> * out_coords,
                                  int64_t num_atoms) {
	return osprey::calc<float64_t,osprey::ambereef1::calc_energy>(d_conf_space, conf, inters, out_coords, num_atoms);
}


namespace osprey {

	#if __CUDA_ARCH__ >= 700 // volta
		#define MINIMIZE_THREADS 1024 // TODO: optimize for volta
		#define MINIMIZE_BLOCKS 1
	#elif __CUDA_ARCH__ >= 600 // pascal
		#define MINIMIZE_THREADS 1024
		#define MINIMIZE_BLOCKS 1
	#elif __CUDA_ARCH__ >= 500 // maxwell
		#define MINIMIZE_THREADS 1024 // TODO: optimize for maxwell?
		#define MINIMIZE_BLOCKS 1
	#else
		// host side
		#define MINIMIZE_THREADS -1
		#define MINIMIZE_BLOCKS -1
	#endif

	__host__
	int get_minimize_threads() {
		int arch = cuda::get_arch();
		if (arch >= 700) { // volta
			return 1024;
		} else if (arch >= 600) { // pascal
			return 1024;
		} else if (arch >= 500) { // maxwell
			return 1024;
		} else {
			throw std::runtime_error("unsupported CUDA architecture");
		}
	}

	template<typename T, EnergyFunction<T> efunc>
	__global__
	__launch_bounds__(MINIMIZE_THREADS, MINIMIZE_BLOCKS)
	void minimize_kernel(const ConfSpace<T> * conf_space,
	                     const Array<int32_t> * conf,
	                     const Array<PosInter<T>> * inters,
	                     Array<Real3<T>> * out_coords,
	                     DofValues<T> * out_dof_values) {

		// slice up the shared memory
		extern __shared__ int8_t shared[];
		int64_t shared_offset = 0;
		auto shared_atom_pairs = reinterpret_cast<const void **>(shared + shared_offset);
		shared_offset += Assignment<T>::sizeof_atom_pairs(conf_space->num_pos);
		auto shared_conf_energies = reinterpret_cast<T *>(shared + shared_offset);
		shared_offset += Assignment<T>::sizeof_conf_energies(conf_space->num_pos);
		auto shared_dofs = reinterpret_cast<Dof<T> **>(shared + shared_offset);
		shared_offset += sizeof(Dof<T> *)*conf_space->max_num_dofs;
		auto shared_line_search_states = reinterpret_cast<LineSearchState<T> *>(shared + shared_offset);
		shared_offset += sizeof(LineSearchState<T>)*conf_space->max_num_dofs;
		auto shared_here = reinterpret_cast<DofValues<T> *>(shared + shared_offset);
		shared_offset += DofValues<T>::get_bytes(conf_space->max_num_dofs);
		auto shared_next = reinterpret_cast<DofValues<T> *>(shared + shared_offset);
		shared_offset += DofValues<T>::get_bytes(conf_space->max_num_dofs);
		auto thread_energy = reinterpret_cast<T *>(shared + shared_offset);

		// init the sizes for cudaMalloc'd Array instances
		out_coords->init(conf_space->max_num_conf_atoms);
		out_dof_values->x.init(conf_space->max_num_dofs);

		// init sizes for Array instances in shared memory
		shared_here->x.init(conf_space->max_num_dofs);
		shared_next->x.init(conf_space->max_num_dofs);

		// make the atoms and the degrees of freedom
		Assignment<T> assignment(*conf_space, *conf, *out_coords, shared_atom_pairs, shared_conf_energies);
		Dofs<T> dofs(assignment, *inters, efunc, thread_energy, shared_dofs);

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

	template<typename T, EnergyFunction<T> efunc>
	static T minimize(const ConfSpace<T> * d_conf_space,
	                  const Array<int32_t> & conf,
	                  const Array<PosInter<T>> & inters,
	                  Array<Real3<T>> * out_coords,
	                  int64_t num_atoms,
	                  Array<T> * out_dofs,
	                  int64_t num_dofs) {

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
			+ Assignment<T>::sizeof_atom_pairs(conf.get_size())
			+ Assignment<T>::sizeof_conf_energies(conf.get_size())
			+ sizeof(Dof<T> *)*num_dofs // thread_energy
			+ sizeof(LineSearchState<T>)*num_dofs // line_search_states
			+ DofValues<T>::get_bytes(num_dofs) // here
			+ DofValues<T>::get_bytes(num_dofs); // next
		int64_t shared_size_per_thread = sizeof(T);

		// launch the kernel
		auto kernel = minimize_kernel<T,efunc>;
		int num_threads = get_minimize_threads();
		int64_t shared_size = shared_size_static + shared_size_per_thread*num_threads;
		kernel<<<1, num_threads, shared_size>>>(d_conf_space, d_conf, d_inters, d_out_coords, d_out_dof_values);
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

API float32_t minimize_amber_eef1_f32(const osprey::ConfSpace<float32_t> * d_conf_space,
                                      const osprey::Array<int32_t> & conf,
                                      const osprey::Array<osprey::PosInter<float32_t>> & inters,
                                      osprey::Array<osprey::Real3<float32_t>> * out_coords,
                                      int64_t num_atoms,
                                      osprey::Array<float32_t> * out_dofs,
                                      int64_t num_dofs) {
	return osprey::minimize<float32_t,osprey::ambereef1::calc_energy>(d_conf_space, conf, inters, out_coords, num_atoms, out_dofs, num_dofs);
}
API float64_t minimize_amber_eef1_f64(const osprey::ConfSpace<float64_t> * d_conf_space,
                                      const osprey::Array<int32_t> & conf,
                                      const osprey::Array<osprey::PosInter<float64_t>> & inters,
                                      osprey::Array<osprey::Real3<float64_t>> * out_coords,
                                      int64_t num_atoms,
                                      osprey::Array<float64_t> * out_dofs,
                                      int64_t num_dofs) {
	return osprey::minimize<float64_t,osprey::ambereef1::calc_energy>(d_conf_space, conf, inters, out_coords, num_atoms, out_dofs, num_dofs);
}
