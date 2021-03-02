
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
			AutoArray<int32_t> modified_posi;

			Dof(T min, T max, T initial_step_size, int max_modified_posi):
					min(min), max(max), initial_step_size(initial_step_size),
					modified_posi(max_modified_posi) {}

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
}


#include "motions/dihedral.h"
#include "motions/transrot.h"

#endif //CONFECALC_MOTIONS_H
