
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
		void DihedralDof<float32_t>::set(Assignment<float32_t> & assignment, void * pdof, float32_t radians, int8_t * shared_mem) {

			auto dof = reinterpret_cast<DihedralDof<float32_t> *>(pdof);
			Real3<float32_t> * atoms = assignment.atoms.items();

			Real3<float32_t> & a = atoms[dof->dihedral->a_index];
			Real3<float32_t> & b = atoms[dof->dihedral->b_index];
			Real3<float32_t> & c = atoms[dof->dihedral->c_index];
			Real3<float32_t> & d = atoms[dof->dihedral->d_index];

			// rotate into a coordinate system where:
			//   b->c is along the -z axis
			//   b->a is in the xz plane
			Rotation<float32_t> r;
			r.set_look(c - b, a - b);

			// rotate about z to set the desired dihedral angle
			float32_t dx = dot<float32_t>(r.xaxis, d - b);
			float32_t dy = dot<float32_t>(r.yaxis, d - b);
			RotationZ<float32_t> r_z;
			r_z.set(dx, dy, -radians);

			// transform all the rotated atoms, in parallel
			for (uint i=threadIdx.x; i<dof->dihedral->num_rotated; i+=blockDim.x) {
				Real3<float32_t> & p = atoms[dof->dihedral->get_rotated_index(i)];
				assert (!isnan3<float32_t>(p));
				p -= b;
				r.mul(p);
				r_z.mul(p);
				r.mul_inv(p); // p = r^T p
				p += b;
				assert (!isnan3<float32_t>(p));
			}
			__syncthreads();
		}

		template<>
		__device__
		void DihedralDof<float64_t>::set(Assignment<float64_t> & assignment, void * pdof, float64_t radians, int8_t * shared_mem) {

			// TODO: go back and clean this up?

			// NOTE: doubles use twice as many registers as floats
			// so this function causes a *lot* of register pressure!
			// this function has been heavily optimized to reduce register usage.
			// that's why it looks so weird!

			auto dof = reinterpret_cast<DihedralDof<float64_t> *>(pdof);

			// slice up the shared memory
			static_assert(osprey::Dofs<float64_t>::dofs_shared_size >=
				sizeof(Real3<float64_t>)
				+ sizeof(Rotation<float64_t>)
				+ sizeof(RotationZ<float64_t>)
			, "using too much shared memory");
			auto b = reinterpret_cast<Real3<float64_t> *>(shared_mem);
			shared_mem += sizeof(Real3<float64_t>);
			auto r = reinterpret_cast<Rotation<float64_t> *>(shared_mem);
			shared_mem += sizeof(Rotation<float64_t>);
			auto r_z = reinterpret_cast<RotationZ<float64_t> *>(shared_mem);

			// copy b to shared mem
			Real3<float64_t> * atoms = assignment.atoms.items();
			if (threadIdx.x == 0) {
				*b = atoms[dof->dihedral->b_index];
			}
			__syncthreads();

			// rotate into a coordinate system where:
			//   b->c is along the -z axis
			//   b->a is in the xz plane
			if (threadIdx.x == 0) {
				Real3<float64_t> & a = atoms[dof->dihedral->a_index];
				Real3<float64_t> & c = atoms[dof->dihedral->c_index];
				r->set_look(c - *b, a - *b);
			}
			__syncthreads();

			// rotate about z to set the desired dihedral angle
			if (threadIdx.x == 0) {
				Real3<float64_t> & d = atoms[dof->dihedral->d_index];
				float64_t dx = dot<float64_t>(r->xaxis, d - *b);
				float64_t dy = dot<float64_t>(r->yaxis, d - *b);
				r_z->set(dx, dy, -radians);
			}
			__syncthreads();

			// transform all the rotated atoms, in parallel
			for (uint i=threadIdx.x; i<dof->dihedral->num_rotated; i+=blockDim.x) {
				Real3<float64_t> & p = atoms[dof->dihedral->get_rotated_index(i)];
				assert (!isnan3<float64_t>(p));
				p -= *b;
				r->mul(p);
				r_z->mul(p);
				r->mul_inv(p); // p = r^T p
				p += *b;
				assert (!isnan3<float64_t>(p));
			}
			__syncthreads();
		}
	}
}
