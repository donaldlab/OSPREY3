
#ifndef CONFECALC_MOTIONS_DIHEDRAL_H
#define CONFECALC_MOTIONS_DIHEDRAL_H


namespace osprey {
	namespace motions {

		template<typename T>
		struct alignas(16) Dihedral {

			const static int64_t Id = 0;

			int64_t id;
			T min_radians;
			T max_radians;
			int32_t a_index;
			int32_t b_index;
			int32_t c_index;
			int32_t d_index;
			int32_t num_rotated;
			int32_t modified_posi;
			// 8 bytes padding if T = float32_t

			// only created by java side
			Dihedral() = delete;

			__host__ __device__
			static int64_t get_motion_bytes() {
				// no extra storage space needed here
				return 0;
			}

			__device__
			inline int32_t get_rotated_index(int i) const {

				// just in case ...
				assert (i >= 0);
				assert (i < num_rotated);

				return get_rotated_indices()[i];
			}

			__device__
			inline const int32_t * get_rotated_indices() const {
				return reinterpret_cast<const int32_t *>(this + 1);
			}
		};
		ASSERT_JAVA_COMPATIBLE_REALS(Dihedral, 48, 48);


		template<typename T>
		class alignas(16) DihedralDof {

			public:

				// only allocated in shared memory
				DihedralDof() = delete;
				DihedralDof(const DihedralDof & other) = delete;
				~DihedralDof() = delete;

				Dof<T> base;
				const Dihedral<T> * dihedral;

				__device__
				void init(const void * desc, const Array<PosInter<T>> & inters) {

					dihedral = reinterpret_cast<const Dihedral<T> *>(desc);

					base.min = dihedral->min_radians;
					base.max = dihedral->max_radians;
					base.initial_step_size = 0.004363323; // 0.25 degrees, in radians
					base.setter = set;
					base.inters = reinterpret_cast<Array<PosInter<T>> *>(this + 1);
					base.filter_inters(inters, dihedral->modified_posi);
				}

				__device__
				static void set(Assignment<T> & assignment, void * dof, T radians, int8_t * shared_mem);
		};
		static_assert(sizeof(DihedralDof<float32_t>) <= Dof<float32_t>::buf_size);
		static_assert(sizeof(DihedralDof<float64_t>) <= Dof<float64_t>::buf_size);
		static_assert(offsetof(DihedralDof<float32_t>, base) == 0);
		static_assert(offsetof(DihedralDof<float64_t>, base) == 0);

		template<typename T>
		__device__
		void make_dofs_dihedral(const void * desc,
		                        int8_t * desc_buf,
		                        int8_t * dof_bufs,
		                        int & size,
		                        Assignment<T> & assignment,
		                        const Array<PosInter<T>> & inters) {
			if (threadIdx.x == 0) {
				auto dof = reinterpret_cast<DihedralDof<T> *>(Dof<T>::get_buf(dof_bufs, size, inters));
				dof->init(desc, inters);
			}
			__syncthreads();
			size += 1;
		}
	}
}


#endif // CONFECALC_MOTIONS_DIHEDRAL_H
