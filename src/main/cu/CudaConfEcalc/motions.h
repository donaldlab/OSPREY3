
#ifndef CONFECALC_MOTIONS_H
#define CONFECALC_MOTIONS_H


namespace osprey {

	template<typename T>
	using DofSetter = void (*)(Assignment<T> & assignment, void * dof, T val, int8_t * shared_mem);

	// degree of freedom
	// uses C style polymorphism (rather than C++ style)
	// since we need to pre-allocate all the shared memory
	template<typename T>
	class alignas(16) Dof {
		public:

			// size of shared memory for each dof
			static const int64_t buf_size = 128;

			// only allocated in shared memory
			Dof() = delete;
			Dof(const Dof & other) = delete;
			~Dof() = delete;

			T min;
			T max;
			T initial_step_size;
			DofSetter<T> setter;
			Array<PosInter<T>> * inters;

			__device__
			inline void set(Assignment<T> & assignment, T val, int8_t * shared_mem) {
				assert (setter != nullptr);
				setter(assignment, this, val, shared_mem);
			}

			__device__
			inline T center() const {
				return (min + max)/2;
			}

			__device__
			const Array<PosInter<T>> & get_inters() const {
				assert (inters != nullptr);
				return *inters;
			}

			__device__
			void make_inters(const Array<PosInter<T>> & all_inters, const int32_t modified_posi[], int num_modified) {

				// how many of the interactions are affected by this degree of freedom?
				int64_t inters_size = 0;
				for (int i=0; i<all_inters.get_size(); i++) {
					if (is_inter_affected(all_inters[i], modified_posi, num_modified)) {
						inters_size += 1;
					}
				}

				assert (inters != nullptr);
				inters->init0(inters_size);

				// copy the interactions
				inters_size = 0;
				for (int i=0; i<all_inters.get_size(); i++) {
					if (is_inter_affected(all_inters[i], modified_posi, num_modified)) {
						(*inters)[inters_size++] = all_inters[i];
					}
				}
			}

			__device__
			static bool is_inter_affected(const PosInter<T> & inter, const int32_t modified_posi[], int num_modified) {

				// is one of the modified positions in this interaction?
				for (int i=0; i<num_modified; i++) {
					int posi = modified_posi[i];
					if (inter.posi1 == posi || inter.posi2 == posi) {
						return true;
					}
				}

				return false;
			};
	};


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
				void init(const Dihedral<T> * _dihedral, const Array<PosInter<T>> & inters) {

					dihedral = _dihedral;

					base.min = dihedral->min_radians;
					base.max = dihedral->max_radians;
					base.initial_step_size = 0.004363323; // 0.25 degrees, in radians
					base.setter = set;
					base.inters = reinterpret_cast<Array<PosInter<T>> *>(this + 1);

					// filter the inters
					base.make_inters(inters, &dihedral->modified_posi, 1);
				}

				__device__
				static void set(Assignment<T> & assignment, void * dof, T radians, int8_t * shared_mem);
		};
		static_assert(sizeof(DihedralDof<float32_t>) <= Dof<float32_t>::buf_size);
		static_assert(sizeof(DihedralDof<float64_t>) <= Dof<float64_t>::buf_size);
		static_assert(offsetof(DihedralDof<float32_t>, base) == 0);
		static_assert(offsetof(DihedralDof<float64_t>, base) == 0);
	}
}


#endif //CONFECALC_MOTIONS_H
