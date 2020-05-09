
#ifndef CONFECALC_MOTIONS_H
#define CONFECALC_MOTIONS_H


namespace osprey {

	// degree of freedom
	// make sure to only construct/destruct on one thread!
	template<typename T>
	class Dof {
		public:

			const T min;
			const T max;
			const T initial_step_size;

			__device__
			Dof(T min, T max, T initial_step_size):
					min(min), max(max), initial_step_size(initial_step_size) {}

			__device__
			Dof(const Dof & other) = delete;
			~Dof() = default;

			__device__
			virtual void set(T val, int8_t * shared_mem) = 0;

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
			inline void free() {
				if (inters != nullptr) {
					std::free(inters);
					inters = nullptr;
				}
			}

		protected:

			__device__
			void make_inters(const Array<PosInter<T>> & all_inters, const int32_t modified_posi[], int num_modified) {

				// how many of the interactions are affected by this degree of freedom?
				int inters_size = 0;
				for (int i=0; i<all_inters.get_size(); i++) {
					if (is_inter_affected(all_inters[i], modified_posi, num_modified)) {
						inters_size += 1;
					}
				}

				// allocate space for the filtered interactions
				inters = Array<PosInter<T>>::make(inters_size);

				// copy the interactions
				inters_size = 0;
				for (int i=0; i<all_inters.get_size(); i++) {
					if (is_inter_affected(all_inters[i], modified_posi, num_modified)) {
						(*inters)[inters_size++] = all_inters[i];
					}
				}
			}

		private:
			Array<PosInter<T>> * inters;

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

			class Dof: public osprey::Dof<T> {
				public:

					static constexpr T step_size = 0.004363323; // 0.25 degrees

					__device__
					Dof(const Dihedral<T> & dihedral, Assignment<T> & assignment, const Array<PosInter<T>> & inters):
						osprey::Dof<T>(dihedral.min_radians, dihedral.max_radians, step_size),
						dihedral(dihedral), assignment(assignment) {

						// filter the inters
						make_inters(inters, &dihedral.modified_posi, 1);
					}

					const Dihedral<T> & dihedral;
					Assignment<T> & assignment;

					__device__
					virtual void set(T radians, int8_t * shared_mem);
			};

			__device__
			inline Dof * make_dof(Assignment<T> & assignment, const Array<PosInter<T>> & inters) const {
				return new Dof(*this, assignment, inters);
			}
		};
		ASSERT_JAVA_COMPATIBLE_REALS(Dihedral, 48, 48);

		template<>
		__device__
		void Dihedral<float32_t>::Dof::set(float32_t radians, int8_t * shared_mem) {

			Real3<float32_t> * atoms = assignment.atoms.items();

			Real3<float32_t> & a = atoms[dihedral.a_index];
			Real3<float32_t> & b = atoms[dihedral.b_index];
			Real3<float32_t> & c = atoms[dihedral.c_index];
			Real3<float32_t> & d = atoms[dihedral.d_index];

			// rotate into a coordinate system where:
			//   b->c is along the -z axis
			//   b->a is in the xz plane
			Rotation<float32_t> r;
			r.set_look(c - b, a - b);

			// rotate about z to set the desired dihedral angle
			float32_t dx = dot<float32_t>(r.xaxis, d - b);
			float32_t dy = dot<float32_t>(r.yaxis, d - b);
			RotationZ<float32_t> r_z;
			r_z.set(dx, dy, -radians);

			// transform all the rotated atoms, in parallel
			for (uint i=threadIdx.x; i<dihedral.num_rotated; i+=blockDim.x) {
				Real3<float32_t> & p = atoms[dihedral.get_rotated_index(i)];
				assert (!isnan3<float32_t>(p));
				p -= b;
				r.mul(p);
				r_z.mul(p);
				r.mul_inv(p); // p = r^T p
				p += b;
				assert (!isnan3<float32_t>(p));
			}
			__syncthreads();
		}

		template<>
		__device__
		void Dihedral<float64_t>::Dof::set(float64_t radians, int8_t * shared_mem) {

			// NOTE: doubles use twice as many registers as floats
			// so this function causes a *lot* of register pressure!
			// this function has been heavily optimized to reduce register usage.
			// that's why it looks so weird!

			// slice up the shared memory
			auto b = reinterpret_cast<Real3<float64_t> *>(shared_mem);
			shared_mem += sizeof(Real3<float64_t>);
			auto r = reinterpret_cast<Rotation<float64_t> *>(shared_mem);
			shared_mem += sizeof(Rotation<float64_t>);
			auto r_z = reinterpret_cast<RotationZ<float64_t> *>(shared_mem);
			shared_mem += sizeof(RotationZ<float64_t>);

			// copy b to shared mem
			Real3<float64_t> * atoms = assignment.atoms.items();
			*b = atoms[dihedral.b_index];

			// rotate into a coordinate system where:
			//   b->c is along the -z axis
			//   b->a is in the xz plane
			if (threadIdx.x == 0) {
				Real3<float64_t> & a = atoms[dihedral.a_index];
				Real3<float64_t> & c = atoms[dihedral.c_index];
				r->set_look(c - *b, a - *b);
			}
			__syncthreads();

			// rotate about z to set the desired dihedral angle
			if (threadIdx.x == 0) {
				Real3<float64_t> & d = atoms[dihedral.d_index];
				float64_t dx = dot<float64_t>(r->xaxis, d - *b);
				float64_t dy = dot<float64_t>(r->yaxis, d - *b);
				r_z->set(dx, dy, -radians);
			}
			__syncthreads();

			// transform all the rotated atoms, in parallel
			for (uint i=threadIdx.x; i<dihedral.num_rotated; i+=blockDim.x) {
				Real3<float64_t> & p = atoms[dihedral.get_rotated_index(i)];
				assert (!isnan3<float64_t>(p));
				p -= *b;
				r->mul(p);
				r_z->mul(p);
				r->mul_inv(p); // p = r^T p
				p += *b;
				assert (!isnan3<float64_t>(p));
			}
			__syncthreads();
		}
	}
}


#endif //CONFECALC_MOTIONS_H
