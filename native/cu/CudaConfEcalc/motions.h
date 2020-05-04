
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
			virtual T get() const = 0;

			__device__
			virtual void set(T val, cg::thread_group threads) = 0;

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
		struct alignas(8) Dihedral {

			const static int32_t id = 0;

			T min_radians;
			T max_radians;
			int32_t a_index;
			int32_t b_index;
			int32_t c_index;
			int32_t d_index;
			int32_t num_rotated;
			int32_t modified_posi;

			// only created by java side
			Dihedral() = delete;

			__device__
			inline int32_t get_rotated_index(int i) const {

				// just in case ...
				assert (i >= 0);
				assert (i < num_rotated);

				return reinterpret_cast<const int32_t *>(this + 1)[i];
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
					virtual T get() const {

						Real3<T> a = assignment.atoms[dihedral.a_index];
						Real3<T> b = assignment.atoms[dihedral.b_index];
						Real3<T> c = assignment.atoms[dihedral.c_index];
						Real3<T> d = assignment.atoms[dihedral.d_index];

						// translate so b is at the origin
						a -= b;
						c -= b;
						d -= b;

						// rotate into a coordinate system where:
						//   b->c is along the -z axis
						//   b->a is in the yz plane
						Rotation<T> r;
						r.set_look(c, a);
						d = r*d;

						return HalfPi<T> - atan2(d.y, d.x);
					}

					__device__
					virtual void set(T radians, cg::thread_group threads) {

						Real3<T> a = assignment.atoms[dihedral.a_index];
						Real3<T> b = assignment.atoms[dihedral.b_index];
						Real3<T> c = assignment.atoms[dihedral.c_index];
						Real3<T> d = assignment.atoms[dihedral.d_index];

						// translate so b is at the origin
						a -= b;
						c -= b;
						d -= b;

						// rotate into a coordinate system where:
						//   b->c is along the -z axis
						//   b->a is in the yz plane
						Rotation<T> r_in;
						r_in.set_look(c, a);
						d = r_in*d;

						// rotate about z to set the desired dihedral angle
						Rotation<T> r_z;
						r_z.set_z(HalfPi<T> - radians - atan2(d.y, d.x));

						// rotate back into the world frame
						Rotation<T> r_out(r_in);
						r_out.invert();

						// transform all the rotated atoms, in parallel
						for (int i=threads.thread_rank(); i<dihedral.num_rotated; i+=threads.size()) {
							Real3<T> & p = assignment.atoms[dihedral.get_rotated_index(i)];
							assert (!isnan3<T>(p));
							p = r_out*(r_z*(r_in*(p - b))) + b;
							assert (!isnan3<T>(p));
						}
						threads.sync();
					}
			};

			__device__
			inline Dof * make_dof(Assignment<T> & assignment, const Array<PosInter<T>> & inters) const {
				return new Dof(*this, assignment, inters);
			}
		};
		ASSERT_JAVA_COMPATIBLE_REALS(Dihedral, 32, 40);
	}
}


#endif //CONFECALC_MOTIONS_H
