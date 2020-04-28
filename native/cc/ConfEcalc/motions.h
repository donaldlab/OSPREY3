
#ifndef CONFECALC_MOTIONS_H
#define CONFECALC_MOTIONS_H


namespace osprey::motions {

	// degree of freedom
	template<typename T>
	class Dof {
		public:

			const T min;
			const T max;
			const T initial_step_size;

			Dof(T min, T max, T initial_step_size, int num_modified_pos):
					min(min), max(max), initial_step_size(initial_step_size),
					num_modified_pos(num_modified_pos), modified_pos(new int32_t[num_modified_pos]) {}

			Dof(const Dof & other) = delete;

			~Dof() {
				delete[] modified_pos;
			}

			virtual T get() const = 0;
			virtual void set(T val) = 0;

			inline T center() const {
				return (min + max)/2;
			}

		private:
			int num_modified_pos;
			int32_t * modified_pos;
	};

	template<typename T>
	struct alignas(8) Dihedral {

		const static int32_t id = 0;

		T min_radians;
		T max_radians;
		int32_t a; // atom indices
		int32_t b;
		int32_t c;
		int32_t d;
		int32_t num_rotated;
		int32_t modified_posi;

		// only created by java side
		Dihedral() = delete;

		inline const int32_t get_rotated(int i) const {

			// just in case ...
			assert (i >= 0);
			assert (i < num_rotated);

			return reinterpret_cast<const int32_t *>(this + 1)[i];
		}

		class Dof: public motions::Dof<T> {
			public:

				static constexpr T step_size = 0.004363323; // 0.25 degrees

				Dof(const Dihedral<T> & dihedral, Assignment<T> & assignment):
					motions::Dof<T>(dihedral.min_radians, dihedral.max_radians, step_size, 1),
					dihedral(dihedral), assignment(assignment) {}

				const Dihedral<T> & dihedral;
				Assignment<T> & assignment;

				virtual T get() const {

					Real3<T> a = assignment.atoms[dihedral.a];
					Real3<T> b = assignment.atoms[dihedral.b];
					Real3<T> c = assignment.atoms[dihedral.c];
					Real3<T> d = assignment.atoms[dihedral.d];

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

					Real3<T> a = assignment.atoms[dihedral.a];
					Real3<T> b = assignment.atoms[dihedral.b];
					Real3<T> c = assignment.atoms[dihedral.c];
					Real3<T> d = assignment.atoms[dihedral.d];

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
						Real3<T> & p = assignment.atoms[dihedral.get_rotated(i)];
						p = r_out*(r_z*(r_in*(p - b))) + b;
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
		int32_t num_atoms;
		int32_t num_modified_pos;

		// only created by java side
		TranslationRotation() = delete;

		inline const int32_t get_atomi(int i) const {

			// just in case ...
			assert (i >= 0);
			assert (i < num_atoms);

			// atom indices are right after the struct
			return reinterpret_cast<const int32_t *>(this + 1)[i];
		}

		inline const int32_t get_modified_posi(int i) const {

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


#endif //CONFECALC_MOTIONS_H
