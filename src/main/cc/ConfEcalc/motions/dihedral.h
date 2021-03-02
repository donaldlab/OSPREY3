
#ifndef CONFECALC_MOTIONS_DIHEDRAL_H
#define CONFECALC_MOTIONS_DIHEDRAL_H


namespace osprey::motions {

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

		inline const int32_t get_rotated_index(int i) const {

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
					Dof::modified_posi.add(dihedral.modified_posi);
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

		inline void make_dofs(AutoArray<osprey::Dof<T> *> * dofs, Assignment<T> & assignment) const {
			dofs->add(new Dof(*this, assignment));
		}
	};
	ASSERT_JAVA_COMPATIBLE_REALS(Dihedral, 32, 40);
}


#endif //CONFECALC_MOTIONS_DIHEDRAL_H
