
#include "global.h"
#include "cuda.h"
#include "array.h"
#include "formats.h"
#include "rotation.h"
#include "atoms.h"
#include "confspace.h"
#include "assignment.h"
/* TEMP
#include "energy.h"
#include "energy_ambereef1.h"
#include "motions.h"
#include "minimization.h"
*/

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

template<typename T>
static void * alloc_conf_space(const osprey::ConfSpace<T> & conf_space) {

	// allocate the device memory
	void * p_device;
	cudaMalloc(&p_device, conf_space.size);
	cuda::check_error();

	// upload the conf space
	cudaMemcpy(p_device, &conf_space, conf_space.size, cudaMemcpyHostToDevice);
	cuda::check_error();

	return p_device;
}

API void * alloc_conf_space_f32(const osprey::ConfSpace<float32_t> & conf_space) {
	return alloc_conf_space(conf_space);
}
API void * alloc_conf_space_f64(const osprey::ConfSpace<float64_t> & conf_space) {
	return alloc_conf_space(conf_space);
}

API void free_conf_space(void * p) {
	cudaFree(p);
	cuda::check_error();
}


template<typename T>
__global__
void assign_kernel(const osprey::ConfSpace<T> * conf_space, const osprey::Array<int32_t> * conf, osprey::Array<osprey::Real3<T>> * out_coords) {

	cg::thread_block threads = cg::this_thread_block();

	// slice up the shared memory
	extern __shared__ int8_t shared[];
	int64_t shared_offset = 0;
	auto shared_atom_pairs = reinterpret_cast<const void **>(shared + shared_offset);
	shared_offset += osprey::Assignment<T>::sizeof_atom_pairs(conf_space->num_pos);
	auto shared_conf_energies = reinterpret_cast<T *>(shared + shared_offset);

	// init the sizes for device-allocated Array instances
	out_coords->init(conf_space->max_num_conf_atoms, threads);

	osprey::Assignment<T> assignment(*conf_space, *conf, *out_coords, shared_atom_pairs, shared_conf_energies, threads);
}

template<typename T>
static void assign(const osprey::ConfSpace<T> * d_conf_space, const osprey::Array<int32_t> & conf, osprey::Array<osprey::Real3<T>> & out_coords) {

	// upload the arguments
	size_t conf_bytes = conf.get_bytes();
	osprey::Array<int32_t> * d_conf;
	cudaMalloc(&d_conf, conf_bytes);
	cuda::check_error();
	cudaMemcpy(d_conf, &conf, conf_bytes, cudaMemcpyHostToDevice);
	cuda::check_error();

	// allocate space for the coords
	size_t out_coords_bytes = out_coords.get_bytes();
	osprey::Array<osprey::Real3<T>> * d_out_coords;
	cudaMalloc(&d_out_coords, out_coords_bytes);
	cuda::check_error();

	// compute the shared memory size
	int64_t shared_size = 0
		+ osprey::Assignment<T>::sizeof_atom_pairs(conf.get_size())
		+ osprey::Assignment<T>::sizeof_conf_energies(conf.get_size());

	// launch the kernel
	int num_threads = cuda::optimize_threads(assign_kernel<T>, shared_size);
	assign_kernel<<<1,num_threads, shared_size>>>(d_conf_space, d_conf, d_out_coords);
	cuda::check_error();

	// download the coords
	cudaMemcpy(&out_coords, d_out_coords, out_coords_bytes, cudaMemcpyDeviceToHost);
	cuda::check_error();

	// cleanup
	cudaFree(d_conf);
	cudaFree(d_out_coords);
}

API void assign_f32(const osprey::ConfSpace<float32_t> * d_conf_space, const osprey::Array<int32_t> & conf, osprey::Array<osprey::Real3<float32_t>> & out_coords) {
	assign(d_conf_space, conf, out_coords);
}
API void assign_f64(const osprey::ConfSpace<float64_t> * d_conf_space, const osprey::Array<int32_t> & conf, osprey::Array<osprey::Real3<float64_t>> & out_coords) {
	assign(d_conf_space, conf, out_coords);
}


/* TODO

template<typename T>
static T calc(const osprey::ConfSpace<T> & conf_space, const int32_t conf[], const osprey::Array<osprey::PosInter<T>> & inters,
              osprey::EnergyFunction<T> efunc, osprey::Array<osprey::Real3<T>> * out_coords) {

	osprey::Assignment<T> assignment(conf_space, conf);
	T energy = efunc(assignment, inters);
	if (out_coords != nullptr) {
		out_coords->copy_from(assignment.atoms);
	}
	return energy;
}

API float32_t calc_amber_eef1_f32(const osprey::ConfSpace<float32_t> & conf_space, const int32_t conf[],
                                  const osprey::Array<osprey::PosInter<float32_t>> & inters,
                                  osprey::Array<osprey::Real3<float32_t>> * out_coords) {
	return calc<float32_t>(conf_space, conf, inters, osprey::ambereef1::calc_energy, out_coords);
}
API float64_t calc_amber_eef1_f64(const osprey::ConfSpace<float64_t> & conf_space, const int32_t conf[],
                                  const osprey::Array<osprey::PosInter<float64_t>> & inters,
                                  osprey::Array<osprey::Real3<float64_t>> * out_coords) {
	return calc<float64_t>(conf_space, conf, inters, osprey::ambereef1::calc_energy, out_coords);
}

template<typename T>
static T minimize(const osprey::ConfSpace<T> & conf_space, const int32_t conf[],
                  const osprey::Array<osprey::PosInter<T>> & inters,
                  osprey::EnergyFunction<T> efunc,
                  osprey::Array<osprey::Real3<T>> * out_coords, osprey::Array<T> * out_dofs) {

	// make the coords and the degrees of freedom
	osprey::Assignment<T> assignment(conf_space, conf);
	osprey::Dofs<T> dofs(assignment, inters, efunc);

	// init the dofs to the center of the voxel
	osprey::DofValues<T> vals(dofs.get_size());
	for (int d=0; d<dofs.get_size(); d++) {
		vals.x[d] = dofs[d].center();
	}

	// use CCD to do the minimization
	osprey::minimize_ccd(dofs, vals);

	// set out the results
	if (out_coords != nullptr) {
		out_coords->copy_from(assignment.atoms);
	}
	if (out_dofs != nullptr) {
		out_dofs->copy_from(vals.x, dofs.get_size());
		out_dofs->truncate(dofs.get_size());
	}
	return vals.f;
}

API float32_t minimize_amber_eef1_f32(const osprey::ConfSpace<float32_t> & conf_space, const int32_t conf[],
                                      const osprey::Array<osprey::PosInter<float32_t>> & inters,
                                      osprey::Array<osprey::Real3<float32_t>> * out_coords, osprey::Array<float32_t> * out_dofs) {
	return minimize<float32_t>(conf_space, conf, inters, osprey::ambereef1::calc_energy, out_coords, out_dofs);
}
API float64_t minimize_amber_eef1_f64(const osprey::ConfSpace<float64_t> & conf_space, const int32_t conf[],
                                      const osprey::Array<osprey::PosInter<float64_t>> & inters,
                                      osprey::Array<osprey::Real3<float64_t>> * out_coords, osprey::Array<float64_t> * out_dofs) {
	return minimize<float64_t>(conf_space, conf, inters, osprey::ambereef1::calc_energy, out_coords, out_dofs);
}

*/
