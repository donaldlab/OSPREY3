
#ifndef CONFECALC_MOTIONS_TRANSROT_H
#define CONFECALC_MOTIONS_TRANSROT_H


namespace osprey::motions {

	template<typename T>
	struct alignas(8) TranslationRotation {

		const static int32_t id = 1;

		T max_distance;
		T max_radians;
		Real3<T> centroid;
		// 4 bytes pad, if T = float32_t
		int32_t moli;
		// 4 bytes pad

		// only created by java side
		TranslationRotation() = delete;

		class TransRotDofs;

		class Dof: public osprey::Dof<T> {
			public:

				Dof(TransRotDofs * transrot, T dist, T step_size):
					osprey::Dof<T>(-dist, dist, step_size, transrot->assignment.conf_space.num_pos),
					transrot(transrot), value(0.0) {}

				~Dof() {
					// cleanup the TransRotDofs if this is the x dof
					if (this == &transrot->dof_x) {
						delete transrot;
					}
				}

				virtual T get() const {
					return value;
				}

				virtual void set(T val) {
					value = val;
					transrot->apply();
				}

			private:
				TransRotDofs * transrot;
				T value;
		};

		// initializes to identity transform
		struct Transform {
			Real3<T> translation;
			Rotation<T> rotation;
		};

		class TransRotDofs {
			public:

				static constexpr T translation_step_size = 0.01; // angstroms
				static constexpr T rotation_step_size = 0.004363323; // 0.25 degrees

				const TranslationRotation<T> & desc;
				Assignment<T> & assignment;
				Dof dof_psi;
				Dof dof_theta;
				Dof dof_phi;
				Dof dof_x;
				Dof dof_y;
				Dof dof_z;

				TransRotDofs(const TranslationRotation<T> & desc, Assignment<T> & assignment):
					desc(desc), assignment(assignment),
					dof_psi(Dof(this, desc.max_radians, rotation_step_size)),
					dof_theta(Dof(this, desc.max_radians, rotation_step_size)),
					dof_phi(Dof(this, desc.max_radians, rotation_step_size)),
					dof_x(Dof(this, desc.max_distance, translation_step_size)),
					dof_y(Dof(this, desc.max_distance, translation_step_size)),
					dof_z(Dof(this, desc.max_distance, translation_step_size)),
					modified_atomi(assignment.conf_space.get_static_atom_coords().get_size() + assignment.conf_space.max_num_conf_atoms)
				{
					// collect the modified positions and atoms

					// collect the static atoms
					for (int atomi=0; atomi<assignment.conf_space.get_static_atom_molis().get_size(); atomi++) {
						int32_t atom_moli = assignment.conf_space.get_static_atom_molis()[atomi];
						if (atom_moli == desc.moli) {
							modified_atomi.add(assignment.get_static_index(atomi));
							add_modified_posi(StaticPos);
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
						for (int atomi=0; atomi<atoms.get_size(); atomi++) {
							int32_t atom_moli = assignment.conf_space.get_conf_atom_molis(conf)[atomi];
							if (atom_moli == desc.moli) {
								modified_atomi.add(assignment.get_index(posi, atomi));
								add_modified_posi(posi);
							}
						}
					}
				}

				void apply() {

					// invert the current transform
					transform_current.translation.negate();
					transform_current.rotation.invert();

					// build the next transform
					Transform transform_next;
					transform_next.translation = { dof_x.get(), dof_y.get(), dof_z.get() };
					transform_next.rotation.set_xyz(dof_psi.get(), dof_theta.get(), dof_phi.get());

					for (int i=0; i<modified_atomi.get_size(); i++) {
						int64_t atomi = modified_atomi[i];
						Real3<T> & p = assignment.atoms[atomi];

						p -= desc.centroid;

						// undo the current transformation
						p = transform_current.rotation*(p + transform_current.translation);

						// apply the next transformation
						p = transform_next.rotation*p + transform_next.translation;

						p += desc.centroid;
					}

					transform_current = transform_next;
				}

			private:
				AutoArray<int32_t> modified_atomi;
				Transform transform_current;

				void add_modified_posi(int32_t posi) {

					// but don't add the same one twice
					// (hopefully linear search here isn't very slow... we only have like 10s of positions at most)
					for (int i=0; i<dof_x.modified_posi.get_size(); i++) {
						if (dof_x.modified_posi[i] == posi) {
							return;
						}
					}

					// add to all the dofs simultaneously
					// NOTE: this order needs to match the java code for the unit tests to pass
					dof_psi.modified_posi.add(posi);
					dof_theta.modified_posi.add(posi);
					dof_phi.modified_posi.add(posi);
					dof_x.modified_posi.add(posi);
					dof_y.modified_posi.add(posi);
					dof_z.modified_posi.add(posi);
				}
		};

		void make_dofs(AutoArray<osprey::Dof<T> *> * dofs, Assignment<T> & assignment) const {
			auto transrot = new TransRotDofs(*this, assignment);
			dofs->add(&transrot->dof_psi);
			dofs->add(&transrot->dof_theta);
			dofs->add(&transrot->dof_phi);
			dofs->add(&transrot->dof_x);
			dofs->add(&transrot->dof_y);
			dofs->add(&transrot->dof_z);
		}
	};
	ASSERT_JAVA_COMPATIBLE_REALS(TranslationRotation, 32, 48);
}


#endif //CONFECALC_MOTIONS_TRANSROT_H
