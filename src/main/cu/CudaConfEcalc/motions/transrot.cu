
#include "../global.h"
#include "../cuda.h"
#include "../array.h"
#include "../formats.h"
#include "../rotation.h"
#include "../atoms.h"
#include "../confspace.h"
#include "../assignment.h"
#include "../energy.h"
#include "../motions.h"
#include "../minimization.h"


namespace osprey {
	namespace motions {

		template<>
		__device__
		void TransRotDof<float32_t>::set(Assignment<float32_t> & assignment, void * pdof, float32_t value, int8_t * shared_mem) {

			auto dof = reinterpret_cast<TransRotDof<float32_t> *>(pdof);
			auto dofs = dof->dofs;

			// save the dof value
			dof->value = value;

			// slice up the shared memory
			// (this fuction uses too many registers, so offload some storage to shared memory)
			static_assert(osprey::Dofs<float32_t>::dofs_shared_size >= sizeof(Transform<float32_t>), "using too much shared memory");
			auto transform_next = reinterpret_cast<Transform<float32_t> *>(shared_mem);

			// build the next transform
			if (threadIdx.x == 0) {
				// TODO: parallelize this?
				transform_next->init(
					dofs->psi->value, dofs->theta->value, dofs->phi->value,
					dofs->x->value, dofs->y->value, dofs->z->value
				);
			}
			__syncthreads();

			// for each masked atom ...
			for (int atomi=threadIdx.x; atomi < assignment.conf_space.max_num_conf_atoms; atomi+=blockDim.x) {
				if (!dofs->atom_mask[atomi]) {
					continue;
				}

				Real3<float32_t> & p = assignment.atoms[atomi];

				p -= dofs->transrot->centroid;

				// undo the current transformation
				p -= dofs->transform_current.translation;
				dofs->transform_current.rotation.mul_inv(p);

				// apply the next transformation
				transform_next->rotation.mul(p);
				p += transform_next->translation;

				p += dofs->transrot->centroid;
			}
			__syncthreads();

			if (threadIdx.x == 0) {
				dofs->transform_current = *transform_next;
			}
			__syncthreads();
		}

		template<>
		__device__
		void TransRotDof<float64_t>::set(Assignment<float64_t> & assignment, void * pdof, float64_t value, int8_t * shared_mem) {

			auto dof = reinterpret_cast<TransRotDof<float64_t> *>(pdof);
			auto dofs = dof->dofs;

			// save the dof value
			dof->value = value;

			// slice up the shared memory
			// (this fuction uses too many registers, so offload some storage to shared memory)
			static_assert(osprey::Dofs<float64_t>::dofs_shared_size >= sizeof(Transform<float64_t>), "using too much shared memory");
			auto transform_next = reinterpret_cast<Transform<float64_t> *>(shared_mem);

			// build the next transform
			if (threadIdx.x == 0) {
				// TODO: parallelize this?
				transform_next->init(
					dofs->psi->value, dofs->theta->value, dofs->phi->value,
					dofs->x->value, dofs->y->value, dofs->z->value
				);
			}
			__syncthreads();

			// for each masked atom ...
			for (int atomi=threadIdx.x; atomi < assignment.conf_space.max_num_conf_atoms; atomi+=blockDim.x) {
				if (!dofs->atom_mask[atomi]) {
					continue;
				}

				Real3<float64_t> & p = assignment.atoms[atomi];

				p -= dofs->transrot->centroid;

				// undo the current transformation
				p -= dofs->transform_current.translation;
				dofs->transform_current.rotation.mul_inv(p);

				// apply the next transformation
				transform_next->rotation.mul(p);
				p += transform_next->translation;

				p += dofs->transrot->centroid;
			}
			__syncthreads();

			// save this transformation as current
			if (threadIdx.x == 0) {
				dofs->transform_current = *transform_next;
			}
			__syncthreads();
		}
	}
}
