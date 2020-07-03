
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
	int err = cudaGetLastError();
	switch (err) {
		case cudaSuccess: break;
		case cudaErrorInsufficientDriver: return -1;
		case cudaErrorNoDevice: return -2;
		case cudaErrorInitializationError: return -3;
		case cudaErrorCompatNotSupportedOnDevice: return -4;
		default:
			std::cout << "Unrecognized CUDA error: " << err << std::endl;
			return INT_MIN;
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

API void * alloc_conf_space_f32(int device, const osprey::ConfSpace<float32_t> & conf_space) {
	return osprey::alloc_conf_space(device, conf_space);
}
API void * alloc_conf_space_f64(int device, const osprey::ConfSpace<float64_t> & conf_space) {
	return osprey::alloc_conf_space(device, conf_space);
}

API void free_conf_space(int device, void * p) {
	osprey::free_conf_space(device, p);
}


API Stream * alloc_stream(int device, int64_t host_bytes, int64_t device_bytes) {
	return osprey::alloc_stream(device, host_bytes, device_bytes);
}
API void free_stream(int device, Stream * stream) {
	osprey::free_stream(device, stream);
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


API void minimize_batch_amber_eef1_f32(int device,
                                       Stream * stream,
                                       const osprey::ConfSpace<float32_t> * d_conf_space,
                                       const osprey::ConfSpaceSizes * conf_space_sizes,
                                       const osprey::MinimizationJobs<float32_t> * jobs,
                                       osprey::Array<float32_t> * out_energies) {
	osprey::minimize_batch<float32_t,osprey::ambereef1::calc_energy>(device, stream, d_conf_space, conf_space_sizes, jobs, out_energies);
}
API void minimize_batch_amber_eef1_f64(int device,
                                       Stream * stream,
                                       const osprey::ConfSpace<float64_t> * d_conf_space,
                                       const osprey::ConfSpaceSizes * conf_space_sizes,
                                       const osprey::MinimizationJobs<float64_t> * jobs,
                                       osprey::Array<float64_t> * out_energies) {
	osprey::minimize_batch<float64_t,osprey::ambereef1::calc_energy>(device, stream, d_conf_space, conf_space_sizes, jobs, out_energies);
}

API int64_t minimize_batch_bufsize_device_f32(const osprey::ConfSpaceSizes * sizes, int64_t batch_size) {
	return osprey::MinimizationBuffers<float32_t>(sizes, batch_size).get_bytes_device();
}
API int64_t minimize_batch_bufsize_device_f64(const osprey::ConfSpaceSizes * sizes, int64_t batch_size) {
	return osprey::MinimizationBuffers<float64_t>(sizes, batch_size).get_bytes_device();
}

API int64_t minimize_batch_bufsize_host_f32(const osprey::ConfSpaceSizes * sizes, int64_t batch_size) {
	return osprey::MinimizationBuffers<float32_t>(sizes, batch_size).get_bytes_host();
}
API int64_t minimize_batch_bufsize_host_f64(const osprey::ConfSpaceSizes * sizes, int64_t batch_size) {
	return osprey::MinimizationBuffers<float64_t>(sizes, batch_size).get_bytes_host();
}


API float32_t minimize_amber_eef1_f32(int device,
                                      Stream * stream,
                                      const osprey::ConfSpace<float32_t> * d_conf_space,
                                      const osprey::ConfSpaceSizes * conf_space_sizes,
                                      const osprey::MinimizationJobs<float32_t> * jobs,
                                      osprey::Array<osprey::Real3<float32_t>> * out_coords,
                                      osprey::Array<float32_t> * out_dofs) {
	return osprey::minimize<float32_t,osprey::ambereef1::calc_energy>(device, stream, d_conf_space, conf_space_sizes, jobs, out_coords, out_dofs);
}
API float64_t minimize_amber_eef1_f64(int device,
                                      Stream * stream,
                                      const osprey::ConfSpace<float64_t> * d_conf_space,
                                      const osprey::ConfSpaceSizes * conf_space_sizes,
                                      const osprey::MinimizationJobs<float64_t> * jobs,
                                      osprey::Array<osprey::Real3<float64_t>> * out_coords,
                                      osprey::Array<float64_t> * out_dofs) {
	return osprey::minimize<float64_t,osprey::ambereef1::calc_energy>(device, stream, d_conf_space, conf_space_sizes, jobs, out_coords, out_dofs);
}
