
#ifndef CONFECALC_MINIMIZATION_H
#define CONFECALC_MINIMIZATION_H


namespace osprey {

	template<typename T>
	class Minimization {
		public:

			Minimization(Assignment<T> & assignment, const PosInter<T> inters[], int64_t inters_size):
				assignment(assignment), inters(inters), inters_size(inters_size) {

				energy = 0.0;
				num_dofs = 0;

				// allocate space for the dofs
				auto max_num_dofs = assignment.conf_space.max_num_dofs;
				dofs = new motions::Dof<T> *[max_num_dofs];
				dof_values = new T[max_num_dofs];

				// make molecule dofs
				for (int motioni=0; motioni < assignment.conf_space.num_molecule_motions; motioni++) {
					switch (assignment.conf_space.get_molecule_motion_id(motioni)) {

						case motions::Dihedral<T>::id: {
							const motions::Dihedral<T> & dihedral = *reinterpret_cast<const motions::Dihedral<T> *>(assignment.conf_space.get_molecule_motion(motioni));
							dofs[num_dofs] = dihedral.make_dof(assignment);
							num_dofs += 1;
						} break;

						default: throw std::invalid_argument("unrecognized motion id");
					}
				}

				// make the conf dofs
				for (int posi=0; posi<assignment.conf_space.num_pos; posi++) {
					const Pos & pos = assignment.conf_space.get_pos(posi);
					const Conf<T> & conf = assignment.conf_space.get_conf(pos, assignment.conf[posi]);
					for (int motioni=0; motioni<conf.num_motions; motioni++) {
						switch (assignment.conf_space.get_conf_motion_id(conf, motioni)) {

							case motions::Dihedral<T>::id: {
								const motions::Dihedral<T> & dihedral = *reinterpret_cast<const motions::Dihedral<T> *>(assignment.conf_space.get_conf_motion(conf, motioni));
								dofs[num_dofs] = dihedral.make_dof(assignment);
								num_dofs += 1;
							} break;

							default: throw std::invalid_argument("unrecognized motion id");
						}
					}
				}
			}

			Minimization(const Minimization & other) = delete;

			~Minimization() {
				delete[] dofs;
			}

			Assignment<T> & assignment;
			const PosInter<T> * const inters;
			const int64_t inters_size;

			T energy;

			inline int get_num_dofs() const {
				return num_dofs;
			}

			inline T get_dof(int i) const {

				// just in case ...
				assert (i >= 0);
				assert (i < num_dofs);

				return dof_values[i];
			}

			void minimize(EnergyFunction<T> efunc) {

				// init the dofs to the center of the voxel
				for (int i=0; i<num_dofs; i++) {
					dof_values[i] = dofs[i]->center();
				}

				// TODO: implement CCD

				// TEMP: just call the energy function
				energy = efunc(assignment, inters, inters_size);
			}

		private:
			int num_dofs;
			motions::Dof<T> ** dofs;
			T * dof_values;
	};
}


#endif //CONFECALC_MINIMIZATION_H
