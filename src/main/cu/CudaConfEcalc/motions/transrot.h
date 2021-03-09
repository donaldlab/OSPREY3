
#ifndef CONFECALC_MOTIONS_TRANSROT_H
#define CONFECALC_MOTIONS_TRANSROT_H


namespace osprey {
	namespace motions {

		template<typename T>
		struct alignas(16) TranslationRotation {

			const static int64_t Id = 1;

			int64_t id;
			T max_distance;
			T max_radians;
			Real3<T> centroid;
			int32_t moli;
			// 12 bytes padding

			// only created by java side
			TranslationRotation() = delete;

			__host__ __device__
			static int64_t get_motion_bytes(int32_t num_pos, int32_t max_num_conf_atoms) {
				// need space to store the indices for the modified positions and the translated/rotated atoms
				return cuda::pad_to_alignment<16>(sizeof(bool)*num_pos)
					+ cuda::pad_to_alignment<16>(sizeof(bool)*max_num_conf_atoms);
			}
		};
		ASSERT_JAVA_COMPATIBLE_REALS(TranslationRotation, 48, 64);

		// need a prototype here since TransRotDofs and TransRotDof form circular references
		template<typename T>
		class alignas(16) TransRotDofs;

		template<typename T>
		class alignas(16) TransRotDof {
			public:

				// only allocated in shared memory
				TransRotDof() = delete;
				TransRotDof(const TransRotDof & other) = delete;
				~TransRotDof() = delete;

				Dof<T> base;
				TransRotDofs<T> * dofs;
				T value;

				__device__
				void init(TransRotDofs<T> * _dofs, const Array<PosInter<T>> & inters, T radius, T step, bool * pos_mask) {

					dofs = _dofs;
					value = (T)0.0;

					base.min = -radius;
					base.max = radius;
					base.initial_step_size = step;
					base.setter = set;
					base.inters = reinterpret_cast<Array<PosInter<T>> *>(this + 1);
					base.filter_inters(inters, pos_mask);
				}

				__device__
				static void set(Assignment<T> & assignment, void * dof, T value, int8_t * shared_mem);
		};
		static_assert(sizeof(TransRotDof<float32_t>) <= Dof<float32_t>::buf_size);
		static_assert(sizeof(TransRotDof<float64_t>) <= Dof<float64_t>::buf_size);
		static_assert(offsetof(TransRotDof<float32_t>, base) == 0);
		static_assert(offsetof(TransRotDof<float64_t>, base) == 0);

		template<typename T>
		struct Transform {

			Real3<T> translation;
			Rotation<T> rotation;

			__device__
			Transform(T rx, T ry, T rz, T dx, T dy, T dz) {
				init(rx, ry, rz, dx, dy, dz);
			}

			__device__
			void init(T rx, T ry, T rz, T dx, T dy, T dz) {
				translation = real3(dx, dy, dz);
				rotation.set_xyz(rx, ry, rz);
			}

			__device__
			void set_identity() {
				translation = real3<T>(0, 0, 0);
				rotation.set_identity();
			}
		};

		template<typename T>
		class alignas(16) TransRotDofs {
			public:

				// only allocated with malloc-style calls
				TransRotDofs() = delete;
				TransRotDofs(const TransRotDofs & other) = delete;
				~TransRotDofs() = delete;

				TransRotDof<T> * psi;
				TransRotDof<T> * theta;
				TransRotDof<T> * phi;
				TransRotDof<T> * x;
				TransRotDof<T> * y;
				TransRotDof<T> * z;

				const TranslationRotation<T> * transrot;
				bool * atom_mask;
				Transform<T> transform_current;

				__device__
				void init(const TranslationRotation<T> * _transrot, bool * _atom_mask) {
					transrot = _transrot;
					atom_mask = _atom_mask;
					transform_current.set_identity();
				}
		};

		template<typename T>
		__device__
		void make_dofs_transrot(const void * desc,
		                        int8_t * desc_buf,
		                        int8_t * dof_bufs,
		                        int & size,
		                        Assignment<T> & assignment,
		                        const Array<PosInter<T>> & inters) {

			auto transrot = reinterpret_cast<const TranslationRotation<T> *>(desc);

			// slice up the descriptor buffer
			size_t offset = 0;
			auto pos_mask = reinterpret_cast<bool *>(desc_buf + offset);
			offset += cuda::pad_to_alignment<16>(sizeof(bool)*assignment.conf_space.num_pos);
			auto atom_mask = reinterpret_cast<bool *>(desc_buf + offset);
			offset += cuda::pad_to_alignment<16>(sizeof(bool)*assignment.conf_space.max_num_conf_atoms);
			auto transrot_dofs = reinterpret_cast<TransRotDofs<T> *>(desc_buf + offset);

			// initialize the masks to false
			for (int i=threadIdx.x; i<assignment.conf_space.num_pos; i+=blockDim.x) {
				pos_mask[i] = false;
			}
			for (int i=threadIdx.x; i<assignment.conf_space.max_num_conf_atoms; i+=blockDim.x) {
				atom_mask[i] = false;
			}

			// collect the static atoms
			for (int atomi=threadIdx.x; atomi<assignment.conf_space.get_static_atom_molis().get_size(); atomi+=blockDim.x) {
				int32_t atom_moli = assignment.conf_space.get_static_atom_molis()[atomi];
				if (atom_moli == transrot->moli) {
					pos_mask[StaticPos] = true;
					atom_mask[assignment.get_static_index(atomi)] = true;
				}
			}

			// collect the conf atoms
			for (int posi=0; posi<assignment.conf_space.num_pos; posi++) {
				const Pos & pos = assignment.conf_space.get_pos(posi);

				// get the assigned conf at this pos, if any
				int32_t confi = assignment.conf[posi];
				if (confi < 0) {
					continue;
				}
				const Conf<T> & conf = assignment.conf_space.get_conf(pos, confi);

				const Array<int32_t> & atoms = assignment.conf_space.get_conf_atom_molis(conf);
				for (int atomi=threadIdx.x; atomi<atoms.get_size(); atomi+=blockDim.x) {
					int32_t atom_moli = assignment.conf_space.get_conf_atom_molis(conf)[atomi];
					if (atom_moli == transrot->moli) {
						pos_mask[posi] = true;
						atom_mask[assignment.get_index(posi, atomi)] = true;
					}
				}
			}

			__syncthreads();

			// build the DoFs
			if (threadIdx.x == 0) {

				transrot_dofs->psi   = reinterpret_cast<TransRotDof<T> *>(Dof<T>::get_buf(dof_bufs, size++, inters));
				transrot_dofs->theta = reinterpret_cast<TransRotDof<T> *>(Dof<T>::get_buf(dof_bufs, size++, inters));
				transrot_dofs->phi   = reinterpret_cast<TransRotDof<T> *>(Dof<T>::get_buf(dof_bufs, size++, inters));
				transrot_dofs->x     = reinterpret_cast<TransRotDof<T> *>(Dof<T>::get_buf(dof_bufs, size++, inters));
				transrot_dofs->y     = reinterpret_cast<TransRotDof<T> *>(Dof<T>::get_buf(dof_bufs, size++, inters));
				transrot_dofs->z     = reinterpret_cast<TransRotDof<T> *>(Dof<T>::get_buf(dof_bufs, size++, inters));

				static constexpr T translation_step_size = 0.01; // angstroms
				static constexpr T rotation_step_size = 0.004363323; // 0.25 degrees

				transrot_dofs->psi  ->init(transrot_dofs, inters, transrot->max_radians, rotation_step_size, pos_mask);
				transrot_dofs->theta->init(transrot_dofs, inters, transrot->max_radians, rotation_step_size, pos_mask);
				transrot_dofs->phi  ->init(transrot_dofs, inters, transrot->max_radians, rotation_step_size, pos_mask);
				transrot_dofs->x    ->init(transrot_dofs, inters, transrot->max_distance, translation_step_size, pos_mask);
				transrot_dofs->y    ->init(transrot_dofs, inters, transrot->max_distance, translation_step_size, pos_mask);
				transrot_dofs->z    ->init(transrot_dofs, inters, transrot->max_distance, translation_step_size, pos_mask);

				transrot_dofs->init(transrot, atom_mask);

			} else {
				size += 6;
			}
			__syncthreads();
		}
	}
}


#endif // CONFECALC_MOTIONS_TRANSROT_H
