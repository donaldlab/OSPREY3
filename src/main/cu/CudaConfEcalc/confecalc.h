
#ifndef CUDACONFECALC_CONFECALC_H
#define CUDACONFECALC_CONFECALC_H


class Stream {
	public:

		__host__
		inline Stream(int64_t host_bytes, int64_t device_bytes) {
			CUDACHECK(cudaStreamCreate(&stream));
			CUDACHECK(cudaMallocHost(&h_buf, host_bytes));
			CUDACHECK(cudaMalloc(&d_buf, device_bytes));
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

namespace osprey {

	void exception_handler();

	template<typename T>
	void * alloc_conf_space(int device, const ConfSpace<T> & conf_space) {

		CUDACHECK(cudaSetDevice(device));

		// allocate the device memory
		void * p_device;
		CUDACHECK(cudaMalloc(&p_device, conf_space.size));

		// upload the conf space
		CUDACHECK(cudaMemcpy(p_device, &conf_space, conf_space.size, cudaMemcpyHostToDevice));

		return p_device;
	}

	void free_conf_space(int device, void * p);


	Stream * alloc_stream(int device, int64_t host_bytes, int64_t device_bytes);
	void free_stream(int device, Stream * stream);

	template<typename T>
	__global__
	void assign_kernel(const ConfSpace<T> * conf_space, const Array<int32_t> * conf, Array<osprey::Real3<T>> * out_coords) {

		// slice up the shared memory
		extern __shared__ int8_t _shared[];
		auto shared = _shared;

		auto shared_index_offsets = reinterpret_cast<int32_t *>(shared);
		assert (cuda::is_aligned<16>(shared_index_offsets));
		shared += cuda::pad_to_alignment<16>(Assignment<T>::sizeof_index_offsets(conf_space->num_pos));

		auto shared_atom_pairs = reinterpret_cast<const void **>(shared);
		assert (cuda::is_aligned<16>(shared_atom_pairs));
		shared += cuda::pad_to_alignment<16>(Assignment<T>::sizeof_atom_pairs(conf_space->num_pos));

		auto shared_conf_energies = reinterpret_cast<T *>(shared);
		assert (cuda::is_aligned<16>(shared_conf_energies));
		shared += cuda::pad_to_alignment<16>(Assignment<T>::sizeof_conf_energies(conf_space->num_pos));

		// init the sizes for cudaMalloc'd Array instances
		out_coords->init(conf_space->max_num_conf_atoms);

		Assignment<T> assignment(*conf_space, *conf, *out_coords, shared_index_offsets, shared_atom_pairs, shared_conf_energies);
	}

	template<typename T>
	void assign(int device,
	            Stream * stream,
	            const ConfSpace<T> * d_conf_space,
	            const Array<int32_t> & conf,
	            Array<Real3<T>> & out_coords) {

		std::set_terminate(exception_handler);

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
			+ cuda::pad_to_alignment<16>(Assignment<T>::sizeof_index_offsets(conf.get_size()))
			+ cuda::pad_to_alignment<16>(Assignment<T>::sizeof_atom_pairs(conf.get_size()))
			+ cuda::pad_to_alignment<16>(Assignment<T>::sizeof_conf_energies(conf.get_size()));
		CUDACHECK_SHAREDSIZE(device, shared_size);

		// launch the kernel
		// TODO: optimize launch bounds
		int num_threads = cuda::optimize_threads(assign_kernel<T>, shared_size, 0);
		assign_kernel<<<1, num_threads, shared_size, stream->get_stream()>>>(d_conf_space, d_conf, d_out_coords);
		CUDACHECK_ERR();

		// download the coords
		CUDACHECK(cudaMemcpy(&out_coords, d_out_coords, out_coords_bytes, cudaMemcpyDeviceToHost));

		// cleanup
		CUDACHECK(cudaFree(d_conf));
		CUDACHECK(cudaFree(d_out_coords));
	}

	#define CALC_THREADS 512
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

		auto shared_index_offsets = reinterpret_cast<int32_t *>(shared);
		assert (cuda::is_aligned<16>(shared_index_offsets));
		shared += cuda::pad_to_alignment<16>(Assignment<T>::sizeof_index_offsets(conf_space->num_pos));

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
		osprey::Assignment<T> assignment(*conf_space, *conf, *out_coords, shared_index_offsets, shared_atom_pairs, shared_conf_energies);

		// call the energy function
		*out_energy = efunc(assignment, *inters, thread_energy);
	}

	template<typename T, osprey::EnergyFunction<T> efunc>
	T calc(int device,
	       Stream * stream,
	       const osprey::ConfSpace<T> * d_conf_space,
	       const osprey::Array<int32_t> & conf,
	       const osprey::Array<osprey::PosInter<T>> & inters,
	       osprey::Array<osprey::Real3<T>> * out_coords,
	       int64_t num_atoms) {

		std::set_terminate(exception_handler);

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
			+ cuda::pad_to_alignment<16>(Assignment<T>::sizeof_index_offsets(conf.get_size()))
			+ cuda::pad_to_alignment<16>(Assignment<T>::sizeof_atom_pairs(conf.get_size()))
			+ cuda::pad_to_alignment<16>(Assignment<T>::sizeof_conf_energies(conf.get_size()));
		int64_t shared_size_per_thread = sizeof(T);

		// launch the kernel
		int num_threads = CALC_THREADS;
		int64_t shared_size = shared_size_static + shared_size_per_thread*num_threads;
		CUDACHECK_SHAREDSIZE(device, shared_size);
		calc_kernel<T,efunc><<<1, num_threads, shared_size, stream->get_stream()>>>(d_conf_space, d_conf, d_inters, d_out_coords, d_out_energy);
		CUDACHECK_ERR();

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


	template<typename T>
	class MinimizationBuffers {
		public:
			int64_t conf_bytes;
			int64_t inters_bytes;
			int64_t coords_bytes;
			int64_t dof_values_bytes;
			int64_t mol_motions_bytes;
			int64_t batch_size;

			__host__
			inline MinimizationBuffers(const ConfSpaceSizes * sizes, int64_t _batch_size) {
				conf_bytes = Array<int32_t>::get_bytes(sizes->num_pos);
				inters_bytes = Array<PosInter<T>>::get_bytes(sizes->max_num_inters);
				coords_bytes = Array<Real3<T>>::get_bytes(sizes->num_atoms);
				dof_values_bytes = DofValues<T>::get_bytes(sizes->max_num_dofs);
				mol_motions_bytes = Dofs<T>::get_mol_motion_bytes(sizes->num_pos, sizes->num_mol_motions, sizes->num_atoms);
				batch_size = _batch_size;
				check_alignments();
			}

			__device__
			inline MinimizationBuffers(const ConfSpace<T> * conf_space, int64_t max_num_inters, int64_t _batch_size) {
				conf_bytes = Array<int32_t>::get_bytes(conf_space->num_pos);
				inters_bytes = Array<PosInter<T>>::get_bytes(max_num_inters);
				coords_bytes = Array<Real3<T>>::get_bytes(conf_space->max_num_conf_atoms);
				dof_values_bytes = DofValues<T>::get_bytes(conf_space->max_num_dofs);
				mol_motions_bytes = Dofs<T>::get_mol_motion_bytes(conf_space->num_pos, conf_space->num_mol_motions, conf_space->max_num_conf_atoms);
				batch_size = _batch_size;
				check_alignments();
			}

			MinimizationBuffers(const MinimizationBuffers & other) = default;
			~MinimizationBuffers() = default;

			__host__
			inline int64_t get_bytes_device() const {
				return (conf_bytes + inters_bytes + coords_bytes + dof_values_bytes + mol_motions_bytes)*batch_size;
			}

			__host__
			inline int64_t get_bytes_host() const {
				return cuda::pad_to_alignment<16>(sizeof(T)*batch_size);
			}

			__host__ __device__
			inline Array<int32_t> * get_conf_buf(int job, int8_t * buf) const {
				return reinterpret_cast<Array<int32_t> *>(get_job_buf(job, buf));
			}

			__host__ __device__
			inline Array<PosInter<T>> * get_inters_buf(int job, int8_t * buf) const {
				return reinterpret_cast<Array<PosInter<T>> *>(get_job_buf(job, buf) + conf_bytes);
			}

			__host__ __device__
			inline Array<Real3<T>> * get_coords_buf(int job, int8_t * buf) const {
				return reinterpret_cast<Array<Real3<T>> *>(get_job_buf(job, buf) + conf_bytes + inters_bytes);
			}

			__host__
			inline T * get_energy_buf(int job, int8_t * buf) const {
				return reinterpret_cast<T *>(buf + sizeof(T)*job);
			}

			__host__ __device__
			inline DofValues<T> * get_dof_values_buf(int job, int8_t * buf) const {
				return reinterpret_cast<DofValues<T> *>(get_job_buf(job, buf) + conf_bytes + inters_bytes + coords_bytes);
			}

			__host__ __device__
			inline int8_t * get_mol_motion_bufs(int job, int8_t * buf) const {
				return get_job_buf(job, buf) + conf_bytes + inters_bytes + coords_bytes + dof_values_bytes;
			}

		private:

			__host__ __device__
			inline void check_alignments() const {
				assert (cuda::is_aligned<16>(conf_bytes));
				assert (cuda::is_aligned<16>(inters_bytes));
				assert (cuda::is_aligned<16>(coords_bytes));
				assert (cuda::is_aligned<16>(dof_values_bytes));
				assert (cuda::is_aligned<16>(mol_motions_bytes));
			}

			__host__ __device__
			inline int8_t * get_job_buf(int job, int8_t * buf) const {
				return buf + (conf_bytes + inters_bytes + coords_bytes + dof_values_bytes + mol_motions_bytes)*job;
			}
	};

	template<typename T>
	class alignas(16) MinimizationJobs {
		public:

			MinimizationJobs() = delete; // created only on the Java side

			__host__
			inline int64_t get_size() const {
				return num_jobs;
			}

			__host__
			inline const Array<int32_t> * get_conf(int j) const {
				auto p = job_offset(j);
				return reinterpret_cast<const Array<int32_t> *>(p);
			}

			__host__
			inline const Array<PosInter<T>> * get_inters(int j) const {
				auto conf = get_conf(j);
				auto p = reinterpret_cast<const int8_t *>(conf) + conf->get_bytes();
				return reinterpret_cast<const Array<PosInter<T>> *>(p);
			}

		private:
			int64_t num_jobs;

			__host__
			inline const int8_t * job_offset(int j) const {
				auto offsets = reinterpret_cast<const Array<int64_t> *>(this);
				return reinterpret_cast<const int8_t *>(this) + (*offsets)[j];
			}
	};
	ASSERT_JAVA_COMPATIBLE_REALS(MinimizationJobs, 16, 16);


	// configure the launch bounds, so the compiler knows how many registers to use
	#define MINIMIZE_THREADS_F32_VOLTA 512
	#define MINIMIZE_THREADS_F64_VOLTA 512
	#define MINIMIZE_THREADS_F32_PASCAL 512 // TODO: optimize for pascal?
	#define MINIMIZE_THREADS_F64_PASCAL 512
	#define MINIMIZE_THREADS_F32_MAXWELL 512 // TODO: optimize for maxwell?
	#define MINIMIZE_THREADS_F64_MAXWELL 512

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
	                          int64_t max_num_inters,
	                          int8_t * xfer_buf) {

		// slice up the shared memory
		extern __shared__ int8_t shared[];
		int64_t offset = 0;

		auto shared_index_offsets = reinterpret_cast<int32_t *>(shared + offset);
		assert (cuda::is_aligned<16>(shared_index_offsets));
		offset += cuda::pad_to_alignment<16>(Assignment<T>::sizeof_index_offsets(conf_space->num_pos));

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
		offset += (Dof<T>::buf_size + Array<const PosInter<T> *>::get_bytes(max_num_inters))*conf_space->max_num_dofs;

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

		// slice up the transfer buf: one minimization job per block
		MinimizationBuffers<T> bufs(conf_space, max_num_inters, gridDim.x);

		// init the sizes for Array instances in the transfer buffer
		bufs.get_coords_buf(blockIdx.x, xfer_buf)->init(conf_space->max_num_conf_atoms);
		bufs.get_dof_values_buf(blockIdx.x, xfer_buf)->x.init(conf_space->max_num_dofs);

		// init sizes for Array instances in shared memory
		shared_here->x.init(conf_space->max_num_dofs);
		shared_next->x.init(conf_space->max_num_dofs);

		// make the atoms and the degrees of freedom
		Assignment<T> assignment(
			*conf_space,
			*bufs.get_conf_buf(blockIdx.x, xfer_buf),
			*bufs.get_coords_buf(blockIdx.x, xfer_buf),
			shared_index_offsets,
			shared_atom_pairs,
			shared_conf_energies
		);
		Dofs<T> dofs(
			assignment,
			*bufs.get_inters_buf(blockIdx.x, xfer_buf),
			bufs.get_mol_motion_bufs(blockIdx.x, xfer_buf),
			efunc,
			thread_energy,
			shared_dof_bufs,
			shared_dofs_mem
		);

		// truncate the dof values to match the actual number of dofs
		if (threadIdx.x == 0) {
			bufs.get_dof_values_buf(blockIdx.x, xfer_buf)->x.truncate(dofs.get_size());
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
		bufs.get_dof_values_buf(blockIdx.x, xfer_buf)->set(*shared_here);
	}

	template<typename T, EnergyFunction<T> efunc>
	__global__
	void minimize_kernel(int64_t shared_size,
	                     const ConfSpace<T> * conf_space,
	                     int64_t max_num_inters,
	                     int8_t * xfer_buf);


	template<typename T>
	__host__
	inline int64_t minimize_kernel_shared_size(const ConfSpaceSizes * conf_space_sizes) {

		// TODO: move quadratically-sized things out of shared memory?
		//  this doesn't scale well to 10s of design positions

		return 0
			+ cuda::pad_to_alignment<16>(Assignment<T>::sizeof_index_offsets(conf_space_sizes->num_pos))
			+ cuda::pad_to_alignment<16>(Assignment<T>::sizeof_atom_pairs(conf_space_sizes->num_pos)) // TODO: quadratically-sized
			+ cuda::pad_to_alignment<16>(Assignment<T>::sizeof_conf_energies(conf_space_sizes->num_pos))
			+ (Dof<T>::buf_size + Array<const PosInter<T> *>::get_bytes(conf_space_sizes->max_num_inters))*conf_space_sizes->max_num_dofs // TODO: quadratically-sized
			+ Dofs<T>::dofs_shared_size
			+ cuda::pad_to_alignment<16>(sizeof(LineSearchState<T>)*conf_space_sizes->max_num_dofs) // line_search_states
			+ DofValues<T>::get_bytes(conf_space_sizes->max_num_dofs) // here
			+ DofValues<T>::get_bytes(conf_space_sizes->max_num_dofs); // next
	}

	template<typename T>
	__host__
	inline int64_t minimize_kernel_shared_size_thread() {
		return sizeof(T);
	}

	template<typename T, EnergyFunction<T> efunc>
	void minimize_batch(int device,
	                    Stream * stream,
	                    const ConfSpace<T> * d_conf_space,
	                    const ConfSpaceSizes * conf_space_sizes,
	                    const MinimizationJobs<T> * jobs,
	                    Array<T> * out_energies) {

		std::set_terminate(exception_handler);

		// make sure we get some jobs, otherwise CUDA will complain with cryptic errors
		if (jobs->get_size() <= 0) {
			throw std::runtime_error("no jobs submitted to batch minimizer");
		}

		CUDACHECK(cudaSetDevice(device));

		// upload the arguments
		MinimizationBuffers<T> buf(conf_space_sizes, jobs->get_size());
		for (int j=0; j<jobs->get_size(); j++) {

			const Array<int32_t> * h_conf = jobs->get_conf(j);
			Array<int32_t> * d_conf = buf.get_conf_buf(j, stream->get_d_buf());
			CUDACHECK(cudaMemcpyAsync(d_conf, h_conf, buf.conf_bytes, cudaMemcpyHostToDevice, stream->get_stream()));

			const Array<PosInter<T>> * h_inters = jobs->get_inters(j);
			Array<PosInter<T>> * d_inters = buf.get_inters_buf(j, stream->get_d_buf());
			CUDACHECK(cudaMemcpyAsync(d_inters, h_inters, buf.inters_bytes, cudaMemcpyHostToDevice, stream->get_stream()));
		}

		// launch the kernel
		int num_threads = get_minimize_threads<T>();
		int64_t shared_size = minimize_kernel_shared_size<T>(conf_space_sizes) + minimize_kernel_shared_size_thread<T>()*num_threads;
		CUDACHECK_SHAREDSIZE(device, shared_size);
		minimize_kernel<T,efunc><<<jobs->get_size(), num_threads, shared_size, stream->get_stream()>>>(
			shared_size,
			d_conf_space,
			conf_space_sizes->max_num_inters,
			stream->get_d_buf()
		);
		CUDACHECK_ERR();

		// download the energies
		for (int j=0; j<jobs->get_size(); j++) {
			auto h_out_energy = buf.get_energy_buf(j, stream->get_h_buf());
			auto d_out_energy = DofValues<T>::f_ptr(buf.get_dof_values_buf(j, stream->get_d_buf()));
			CUDACHECK(cudaMemcpyAsync(h_out_energy, d_out_energy, sizeof(T), cudaMemcpyDeviceToHost, stream->get_stream()));
		}

		// wait for the downloads to finish
		CUDACHECK(cudaStreamSynchronize(stream->get_stream()));

		// finally, copy out the minimized energies
		for (int j=0; j<jobs->get_size(); j++) {
			auto h_out_energy = buf.get_energy_buf(j, stream->get_h_buf());
			(*out_energies)[j] = *reinterpret_cast<const T *>(h_out_energy);
		}
	}


	template<typename T, EnergyFunction<T> efunc>
	T minimize(int device,
	           Stream * stream,
	           const ConfSpace<T> * d_conf_space,
	           const ConfSpaceSizes * conf_space_sizes,
	           const MinimizationJobs<T> * jobs,
	           Array<Real3<T>> * out_coords,
	           Array<T> * out_dofs) {

		std::set_terminate(exception_handler);

		CUDACHECK(cudaSetDevice(device));

		// upload the arguments
		MinimizationBuffers<T> buf(conf_space_sizes, 1);

		const Array<int32_t> * h_conf = jobs->get_conf(0);
		Array<int32_t> * d_conf = buf.get_conf_buf(0, stream->get_d_buf());
		CUDACHECK(cudaMemcpyAsync(d_conf, h_conf, buf.conf_bytes, cudaMemcpyHostToDevice, stream->get_stream()));

		const Array<PosInter<T>> * h_inters = jobs->get_inters(0);
		Array<PosInter<T>> * d_inters = buf.get_inters_buf(0, stream->get_d_buf());
		CUDACHECK(cudaMemcpyAsync(d_inters, h_inters, buf.inters_bytes, cudaMemcpyHostToDevice, stream->get_stream()));

		// launch the kernel
		int num_threads = get_minimize_threads<T>();
		int64_t shared_size = minimize_kernel_shared_size<T>(conf_space_sizes) + minimize_kernel_shared_size_thread<T>()*num_threads;
		CUDACHECK_SHAREDSIZE(device, shared_size);
		minimize_kernel<T,efunc><<<jobs->get_size(), num_threads, shared_size, stream->get_stream()>>>(
			shared_size,
			d_conf_space,
			conf_space_sizes->max_num_inters,
			stream->get_d_buf()
		);
		CUDACHECK_ERR();

		// download the outputs
		auto h_out_energy = buf.get_energy_buf(0, stream->get_h_buf());
		auto d_out_energy = DofValues<T>::f_ptr(buf.get_dof_values_buf(0, stream->get_d_buf()));
		CUDACHECK(cudaMemcpyAsync(h_out_energy, d_out_energy, sizeof(T), cudaMemcpyDeviceToHost, stream->get_stream()));

		auto d_out_dofs = DofValues<T>::x_ptr(buf.get_dof_values_buf(0, stream->get_d_buf()));
		CUDACHECK(cudaMemcpyAsync(out_dofs, d_out_dofs, Array<T>::get_bytes(conf_space_sizes->max_num_dofs), cudaMemcpyDeviceToHost, stream->get_stream()));

		auto d_out_coords = buf.get_coords_buf(0, stream->get_d_buf());
		CUDACHECK(cudaMemcpyAsync(out_coords, d_out_coords, buf.coords_bytes, cudaMemcpyDeviceToHost, stream->get_stream()));

		// wait for the downloads to finish
		CUDACHECK(cudaStreamSynchronize(stream->get_stream()));

		// finally, return the energy
		return *h_out_energy;
	}
}


#endif //CUDACONFECALC_CONFECALC_H
