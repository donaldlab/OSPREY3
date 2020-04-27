
#ifndef CONFECALC_MINIMIZATION_H
#define CONFECALC_MINIMIZATION_H


namespace osprey {

	// degree of freedom
	template<typename T>
	struct alignas(8) Dof {
		T min;
		T max;
		T initial_step_size;
		uint32_t num_modified_pos;
		// 4 bytes pad, if T = float32_t
	};
	ASSERT_JAVA_COMPATIBLE_REALS(Dof, 16, 32);

	template<typename T>
	class Minimization {
		public:

			Minimization(Assignment<T> & assignment, const PosInter<T> inters[], int64_t inters_size):
				assignment(assignment), inters(inters), inters_size(inters_size) {

				energy = 0.0;

				// TODO: make all the DoFs

				num_dofs = 5;
				dofs = new T[num_dofs];
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

				return dofs[i];
			}

			void minimize(EnergyFunction<T> efunc) {

				// TEMP
				for (int i=0; i<num_dofs; i++) {
					dofs[i] = 5.6;
				}

				// TODO: init dof values to centers

				energy = efunc(assignment, inters, inters_size);
			}

		private:
			int num_dofs;
			T * dofs;
	};
}


#endif //CONFECALC_MINIMIZATION_H
