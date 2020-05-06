
#ifndef CONFECALC_MINIMIZATION_H
#define CONFECALC_MINIMIZATION_H


namespace osprey {

	template<typename T>
	class DofValues {
		public:

			// created only by malloc-style allocations
			__device__
			DofValues() = delete;
			__device__
			DofValues(const DofValues<T> & other) = delete;
			~DofValues() = default;

			__device__
			inline void set(const DofValues<T> & other) {

				// TODO: do all the copies simultaneously, but on different threads?
				if (threadIdx.x == 0) {
					f = other.f;
				}
				__syncthreads();

				assert (x.get_size() == other.x.get_size());
				x.copy_from_device(other.x);
			}

			T f;
			Array<T> x;

			__device__
			static inline DofValues<T> make(int size, void * shared_ptr) {
				if (threadIdx.x == 0) {
					shared_ptr = std::malloc(get_bytes(size));
				}
				__syncthreads();
				return reinterpret_cast<DofValues<T> *>(shared_ptr);
			}

			__host__ __device__
			static inline int64_t get_bytes(int size) {
				return 8 + Array<T>::get_bytes(size);
			}

			__host__
			static inline const T * f_ptr(const DofValues<T> * p) {
				return reinterpret_cast<const T *>(p);
			}

			__host__
			static inline const T * x_ptr(const DofValues<T> * p) {
				return reinterpret_cast<const T *>(reinterpret_cast<const int8_t *>(p) + 8);
			}
	};
	ASSERT_MALLOCABLE_REALS(DofValues, 8 + sizeof(Array<float32_t>), 8 + sizeof(Array<float64_t>));

	template<typename T>
	class Dofs {
		public:

			Assignment<T> & assignment;
			const Array<PosInter<T>> & inters;
			const EnergyFunction<T> efunc;
			T * thread_energy;

			__device__
			Dofs(Assignment<T> & assignment,
			     const Array<PosInter<T>> & inters,
			     EnergyFunction<T> efunc,
			     T thread_energy[],
			     Dof<T> * shared_dofs[]):
					assignment(assignment), efunc(efunc), thread_energy(thread_energy), inters(inters), shared_dofs(shared_dofs) {

				size = 0;

				// TODO: parallelize, so thread 0 isn't doing all the work

				// make molecule dofs
				for (int motioni=0; motioni < assignment.conf_space.num_molecule_motions; motioni++) {
					switch (assignment.conf_space.get_molecule_motion_id(motioni)) {

						case motions::Dihedral<T>::id: {
							if (threadIdx.x == 0) {
								const motions::Dihedral<T> & dihedral = *reinterpret_cast<const motions::Dihedral<T> *>(assignment.conf_space.get_molecule_motion(motioni));
								shared_dofs[size] = dihedral.make_dof(assignment, inters);
							}
							size += 1;
						} break;

						// TODO: translation/rotation

						default: assert(false);
					}
				}

				// make the conf dofs
				for (int posi=0; posi<assignment.conf_space.num_pos; posi++) {
					const Pos & pos = assignment.conf_space.get_pos(posi);
					const Conf<T> & conf = assignment.conf_space.get_conf(pos, assignment.conf[posi]);
					for (int motioni=0; motioni<conf.num_motions; motioni++) {
						switch (assignment.conf_space.get_conf_motion_id(conf, motioni)) {

							case motions::Dihedral<T>::id: {
								if (threadIdx.x == 0) {
									const motions::Dihedral<T> & dihedral = *reinterpret_cast<const motions::Dihedral<T> *>(assignment.conf_space.get_conf_motion(conf, motioni));
									shared_dofs[size] = dihedral.make_dof(assignment, inters);
								}
								size += 1;
							} break;

							default: assert(false);
						}
					}
				}

				// sync the writes to shared_dofs
				__syncthreads();
			}

			__device__
			Dofs(const Dofs<T> & other) = delete;
			~Dofs() = default;

			__device__
			inline void free() const {
				for (int d=threadIdx.x; d<size; d+=blockDim.x) {
					shared_dofs[d]->free();
				}
				__syncthreads();
			}

			__device__
			inline int get_size() const {
				return size;
			}

			__device__
			Dof<T> & operator [] (int i) {
				return *shared_dofs[i];
			}

			__device__
			inline void set(const Array<T> & x) {

				// WARNING: it may be tempting to set dofs in parallel here
				// but don't do it: different dofs can move the same atoms! it'll race!

				assert (x.get_size() == size);
				for (int d=0; d<x.get_size(); d++) {
					assert (shared_dofs[d] != nullptr);
					shared_dofs[d]->set(x[d]);
				}
			}

			__device__
			inline T eval_efunc(Array<T> & x) {
				set(x);
				return efunc(assignment, inters, thread_energy);
			}

			__device__
			inline T eval_efunc(int d, T x) {
				Dof<T> & dof = *shared_dofs[d];
				dof.set(x);
				return efunc(assignment, dof.get_inters(), thread_energy);
			}

		private:
			int size;
			Dof<T> ** shared_dofs;
	};


	template<typename T>
	static const T tolerance;
	template<>
	const float32_t tolerance<float32_t> = 1e-3;
	template<>
	const float64_t tolerance<float64_t> = 1e-6;

	// scale abs(f) by tolerance, unless f is very small
	template<typename T>
	__device__
	static T scaled_tolerance(T f) {
		return tolerance<T>*max(static_cast<T>(1.0), std::abs(f));
	}


	// search the line by fitting a local quadratic model, taking a step, and then surfing the slope
	template<typename T>
	__device__
	static T line_search_surf(Dofs<T> & dofs, int d, T x, T & step) {

		auto f = [&dofs, d](T x) -> T {
			return dofs.eval_efunc(d, x);
		};

		T fx = f(x);

		T fxmin = NAN;
		T fxmax = NAN;

		// make sure the step isn't so big that the quadratic approximation is worthless
		T xmin = dofs[d].min;
		T xmax = dofs[d].max;
		while (x - step < xmin && x + step > xmax) {
			step /= 2;
		}

		// get the positive (p) and negative (n) neighbors for our current pos
		T xp = x + step;
		T xm = x - step;
		T fxp = INFINITY;
		if (xp <= xmax) {
			fxp = f(xp);
		}
		T fxm = INFINITY;
		if (xm >= xmin) {
			fxm = f(xm);
		}

		// fit a quadratic to the objective function, locally:
		// q(x) = fx + a*(x - xd)^2 + b*(x - xd)
		// a*step^2 + b*step = fxp - fx
		// a*step^2 - b*step = fxm - fx

		// solve for the shape of the parabola
		T shape = fxp + fxm - 2*fx;
		const T shape_epsilon = 1e-12;
		T xstar;
		if ((shape < -shape_epsilon) || isnan<T>(shape) || isinf<T>(shape)) {

			// negative shape means quadratic is concave down
			// infinite or nan a means we're hitting a constraint or impossible conformation
			// so just minimize over the endpoints of the interval
			if (fxm < fxp) {
				xstar = xm;
			} else {
				xstar = xp;
			}

		} else if (shape <= shape_epsilon) {

			// shape near zero means it's basically flat here
			// so don't step anywhere
			xstar = x;

		} else {

			// positive shape means quadratic is concave up
			// step to the optimum
			xstar = x + (fxm - fxp)*step/2/shape;
		}

		xstar = clamp(xstar, xmin, xmax);
		T fxstar = f(xstar);

		// did we go downhill?
		if (fxstar < fx) {

			// surf along f locally to try to find better minimum
			T xsurfHere = xstar;
			T fxsurfHere = fxstar;
			while (true) {

				// take a step twice as far as we did last time
				T xsurfNext = x + 2*(xsurfHere - x);

				// did we step off the min?
				if (xsurfNext < xmin) {

					// if the min is better, go there instead
					if (isnan<T>(fxmin)) {
						fxmin = f(xmin);
					}
					if (fxmin < fxsurfHere) {
						xsurfHere = xmin;
						fxsurfHere = fxmin;
					}

					break;

				// did we step off the max?
				} else if (xsurfNext > xmax) {

					// if the max is better, go there instead
					if (isnan<T>(fxmax)) {
						fxmax = f(xmax);
					}
					if (fxmax < fxsurfHere) {
						xsurfHere = xmax;
						fxsurfHere = fxmax;
					}

					break;
				}

				T fxsurfNext = f(xsurfNext);

				// did we improve the min enough to keep surfing?
				if (fxsurfNext < fxsurfHere - scaled_tolerance(fxsurfHere)) {

					// yeah, keep going
					xsurfHere = xsurfNext;
					fxsurfHere = fxsurfNext;

				} else {

					// nope, stop surfing
					break;
				}
			}

			// update the minimum estimate so far
			xstar = xsurfHere;
			fxstar = fxsurfHere;

		// did we go significantly uphill?
		} else if (fxstar > fx + tolerance<T>) {

			// try to surf back downhill
			T xsurfHere = xstar;
			T fxsurfHere = fxstar;
			while (true) {

				// cut the step in half
				T xsurfNext = x + (xsurfHere - x)/2;
				T fxsurfNext = f(xsurfNext);

				// did we improve the min enough to keep surfing?
				if (fxsurfNext < fxsurfHere - scaled_tolerance(fxsurfHere)) {

					// yeah, keep going
					xsurfHere = xsurfNext;
					fxsurfHere = fxsurfNext;

				} else {

					// nope, stop surfing
					break;
				}
			}

			// did the quadratic step help at all?
			if (fxstar < fx) {

				// yeah, keep it!

			} else {

				// nope, the original spot was lower
				xstar = x;
				fxstar = fx;
			}

			// did surfing help at all?
			if (fxsurfHere < fxstar) {

				// yeah, use the surf spot
				xstar = xsurfHere;
				fxstar = fxsurfHere;
			}
		}

		// update step before wall jumping
		step = xstar - x;

		// try to jump over walls arbitrarily
		// look in a 1-degree step for a better minimum

		// NOTE: skipping this can make minimization a bit faster,
		// but skipping this causes a noticeable rise in final energies too
		// it's best to keep doing it I think

		xm = xstar - 1;
		xp = xstar + 1;

		if (xm >= xmin) {
			fxm = f(xm);
			if (fxm < fxstar) {
				xstar = xm;
				fxstar = fxm;
			}
		}

		if (xp <= xmax) {
			fxp = f(xp);
			if (fxp < fxstar) {
				xstar = xp;
				fxstar = fxp;
			}
		}

		// set the coords back to the best option
		dofs[d].set(xstar);

		return xstar;
	}

	template<typename T>
	struct LineSearchState {
		T first_step = 1.0;
		T last_step = 1.0;
	};

	template<typename T>
	using LineSearchFunction = T (*)(Dofs<T> &, int, T, T&);

	// TODO: use launch bounds to limit register usage?
	//__launch_bounds__(maxThreadsPerBlock, minBlocksPerMultiprocessor)
	// see: https://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html#launch-bounds

	template<typename T>
	__device__
	static void minimize_ccd(Dofs<T> & dofs,
	                         DofValues<T> & here,
	                         DofValues<T> & next,
	                         LineSearchState<T> shared_line_search_states[]) {

		// get the current objective function value
		here.f = dofs.eval_efunc(here.x);

		// ccd is pretty simple actually
		// just do a line search along each dimension until we stop improving
		// we deal with cycles by just capping the number of iterations

		const int max_iterations = 30;
		const T convergence_threshold = 0.001;

		const LineSearchFunction<T> line_search = line_search_surf;

		for (int iter=0; iter<max_iterations; iter++) {

			// update all the dofs using line search
			next.set(here);
			for (int d=0; d<dofs.get_size(); d++) {

				Dof<T> & dof = dofs[d];
				LineSearchState<T> & state = shared_line_search_states[d];

				// get the step size, try to make it adaptive (based on historical steps if possible; else on step #)
				T step;
				if (std::abs(state.last_step) > tolerance<T> && std::abs(state.first_step) > tolerance<T>) {
					step = dof.initial_step_size*std::abs(state.last_step/state.first_step);
				} else {
					step = dof.initial_step_size/std::pow(iter + 1, 3);
				}

				// get the next x value for this dof
				next.x[d] = line_search(dofs, d, next.x[d], step);

				// update step stats
				if (threadIdx.x == 0) {
					if (iter == 0) {
						state.first_step = step;
					}
					state.last_step = step;
				}
				__syncthreads();
			}

			// how much did we improve?
			next.f = dofs.eval_efunc(next.x);
			T improvement = here.f - next.f;
			if (improvement > 0) {

				// take the step
				here.set(next);

				if (improvement < convergence_threshold) {
					break;
				}

			} else {

				// don't take the step, revert the coords
				dofs.set(here.x);

				break;
			}
		}
	}
}


#endif //CONFECALC_MINIMIZATION_H
