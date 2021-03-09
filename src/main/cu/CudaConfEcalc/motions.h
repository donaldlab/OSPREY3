
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

			__device__
			static Dof<T> * get_buf(int8_t * dof_bufs, int d, const Array<PosInter<T>> & inters) {
				int64_t size = buf_size + Array<const PosInter<T> *>::get_bytes(inters.get_size());
				return reinterpret_cast<Dof<T> *>(dof_bufs + size*d);
			}

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
			void filter_inters(const Array<PosInter<T>> & all_inters, const int32_t modified_posi) {

				assert (inters != nullptr);
				inters->init0(all_inters.get_size());

				// copy the interactions
				int64_t size = 0;
				for (int i=0; i<all_inters.get_size(); i++) {
					const PosInter<T> & inter = all_inters[i];
					if (inter.posi1 == modified_posi || inter.posi2 == modified_posi) {
						(*inters)[size++] = all_inters[i];
					}
				}

				inters->truncate(size);
			}

			__device__
			void filter_inters(const Array<PosInter<T>> & all_inters, const bool pos_mask[]) {

				assert (inters != nullptr);
				inters->init0(all_inters.get_size());

				// copy the interactions
				int64_t size = 0;
				for (int i=0; i<all_inters.get_size(); i++) {
					const PosInter<T> & inter = all_inters[i];
					if (pos_mask[inter.posi1] || pos_mask[inter.posi2]) {
						(*inters)[size++] = all_inters[i];
					}
				}

				inters->truncate(size);
			}
	};
}


#include "motions/dihedral.h"
#include "motions/transrot.h"


#endif //CONFECALC_MOTIONS_H
