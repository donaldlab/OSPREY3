
#ifndef CONFECALC_MOTIONS_H
#define CONFECALC_MOTIONS_H


namespace osprey {

	// degree of freedom
	template<typename T>
	class Dof {
		public:

			const T min;
			const T max;
			const T initial_step_size;
			Array<int32_t> modified_posi;

			Dof(T min, T max, T initial_step_size, int num_modified_posi):
					min(min), max(max), initial_step_size(initial_step_size),
					modified_posi(num_modified_posi) {}

			Dof(const Dof & other) = delete;

			~Dof() {
				if (inters != nullptr) {
					delete inters;
				}
			}

			virtual T get() const = 0;
			virtual void set(T val) = 0;

			inline T center() const {
				return (min + max)/2;
			}

			void set_inters(const Array<PosInter<T>> & all_inters) {

				// how many of the interactions are affected by this degree of freedom?
				int inters_size = 0;
				for (int i=0; i<all_inters.get_size(); i++) {
					if (is_inter_affected(all_inters[i])) {
						inters_size += 1;
					}
				}

				// allocate space for the filtered interactions
				inters = new Array<PosInter<T>>(inters_size);

				// copy the interactions
				inters_size = 0;
				for (int i=0; i<all_inters.get_size(); i++) {
					if (is_inter_affected(all_inters[i])) {
						(*inters)[inters_size++] = all_inters[i];
					}
				}
			}

			const Array<PosInter<T>> & get_inters() const {
				assert (inters != nullptr);
				return *inters;
			}

		private:
			Array<PosInter<T>> * inters;

			bool is_inter_affected(const PosInter<T> & inter) const {

				// is one of the modified positions in this interaction?
				for (int i=0; i < modified_posi.get_size(); i++) {
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

			inline int32_t get_rotated_index(int i) const {

				// just in case ...
				assert (i >= 0);
				assert (i < num_rotated);

				return reinterpret_cast<const int32_t *>(this + 1)[i];
			}

			class Dof: public osprey::Dof<T> {
				public:

					static constexpr T step_size = 0.004363323; // 0.25 degrees

					Dof(const Dihedral<T> & dihedral, Assignment<T> & assignment):
						osprey::Dof<T>(dihedral.min_radians, dihedral.max_radians, step_size, 1),
						dihedral(dihedral), assignment(assignment) {

						// set the modified position indices
						Dof::modified_posi[0] = dihedral.modified_posi;
					}

					const Dihedral<T> & dihedral;
					Assignment<T> & assignment;

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

					virtual void set(T radians) {

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

						// transform all the rotated atoms
						for (int i=0; i<dihedral.num_rotated; i++) {
							Real3<T> & p = assignment.atoms[dihedral.get_rotated_index(i)];
							assert (!p.isnan());
							p = r_out*(r_z*(r_in*(p - b))) + b;
							assert (!p.isnan());
						}
					}
			};

			inline Dof * make_dof(Assignment<T> & assignment) const {
				return new Dof(*this, assignment);
			}
		};
		ASSERT_JAVA_COMPATIBLE_REALS(Dihedral, 32, 40);

		template<typename T>
		struct alignas(8) TranslationRotation {

			const static int32_t id = 1;

			T max_distance;
			T max_radians;
			Real3<T> centroid;
			int32_t num_atoms; // TODO: switch to array
			int32_t num_modified_pos;

			// only created by java side
			TranslationRotation() = delete;

			inline int32_t get_atomi(int i) const {

				// just in case ...
				assert (i >= 0);
				assert (i < num_atoms);

				// atom indices are right after the struct
				return reinterpret_cast<const int32_t *>(this + 1)[i];
			}

			inline int32_t get_modified_posi(int i) const {

				// just in case ...
				assert (i >= 0);
				assert (i > num_modified_pos);

				// modified position indices are right after the atom indices
				return reinterpret_cast<const int32_t *>(this + 1)[num_atoms + i];
			}

			// TODO: implement DoFs
		};
		ASSERT_JAVA_COMPATIBLE_REALS(TranslationRotation, 32, 48);
	}
}


#endif //CONFECALC_MOTIONS_H
