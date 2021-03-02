
#ifndef CONFECALC_MINIMIZATION_H
#define CONFECALC_MINIMIZATION_H


namespace osprey {

	template<typename T>
	class DofValues {
		public:

			DofValues(int size): f(0.0), x(size) {}
			DofValues(const DofValues<T> & other) : DofValues(other.size) {
				set(other);
			}
			~DofValues() = default;

			inline void set(const DofValues<T> & other) {
				f = other.f;
				x.copy_from(other.x);
			}

			T f;
			Array<T> x;
	};

	template<typename T>
	class Dofs {
		public:

			Assignment<T> & assignment;
			const Array<PosInter<T>> & inters;
			const EnergyFunction<T> efunc;

			Dofs(Assignment<T> & assignment, const Array<PosInter<T>> & inters, EnergyFunction<T> efunc):
				assignment(assignment), efunc(efunc), inters(inters) {

				// allocate space for the dofs
				dofs = new AutoArray<Dof<T> *>(assignment.conf_space.max_num_dofs);

				// make molecule dofs
				for (int motioni=0; motioni < assignment.conf_space.num_molecule_motions; motioni++) {
					switch (auto motionid = assignment.conf_space.get_molecule_motion_id(motioni)) {

						case motions::Dihedral<T>::id: {
							const motions::Dihedral<T> & dihedral = *reinterpret_cast<const motions::Dihedral<T> *>(assignment.conf_space.get_molecule_motion(motioni));
							dihedral.make_dofs(dofs, assignment);
						} break;

						case motions::TranslationRotation<T>::id: {
							const motions::TranslationRotation<T> & transrot = *reinterpret_cast<const motions::TranslationRotation<T> *>(assignment.conf_space.get_molecule_motion(motioni));
							transrot.make_dofs(dofs, assignment);
						} break;

						default: {
							std::cout << "motion id: " << motionid << std::endl;
							throw std::invalid_argument("unrecognized motion id for molecule");
						}
					}
				}

				// make the conf dofs
				for (int posi=0; posi<assignment.conf_space.num_pos; posi++) {
					const Pos & pos = assignment.conf_space.get_pos(posi);

					// is this pos assigned?
					int32_t confi = assignment.conf[posi];
					if (confi >= 0) {
						const Conf<T> & conf = assignment.conf_space.get_conf(pos, confi);

						// yup, make the dofs
						for (int motioni=0; motioni<conf.num_motions; motioni++) {
							switch (auto motionid = assignment.conf_space.get_conf_motion_id(conf, motioni)) {

								case motions::Dihedral<T>::id: {
									const motions::Dihedral<T> & dihedral = *reinterpret_cast<const motions::Dihedral<T> *>(assignment.conf_space.get_conf_motion(conf, motioni));
									dihedral.make_dofs(dofs, assignment);
								} break;

								default: {
									std::cout << "motion id: " << motionid << std::endl;
									throw std::invalid_argument("unrecognized motion id for conformation");
								}
							}
						}
					}
				}

				// make the inters for each dof
				for (int d=0; d<dofs->get_size(); d++) {
					(*dofs)[d]->set_inters(inters);
				}
			}

			Dofs(const Dofs<T> & other) = delete;

			~Dofs() {
				delete dofs;
			}

			inline int get_size() const {
				return dofs->get_size();
			}

			Dof<T> & operator [] (int i) {
				return *(*dofs)[i];
			}

			inline void set(const Array<T> & x) {
				for (int d=0; d<x.get_size(); d++) {
					(*dofs)[d]->set(x[d]);
				}
			}

			inline T eval_efunc(Array<T> & x) {
				set(x);
				return efunc(assignment, inters);
			}

			inline T eval_efunc(int d, T x) {
				Dof<T> & dof = *(*dofs)[d];
				dof.set(x);
				return efunc(assignment, dof.get_inters());
			}

		private:
			AutoArray<Dof<T> *> * dofs;
	};


	template<typename T>
	static const T tolerance;
	template<>
	const float32_t tolerance<float32_t> = 1e-3;
	template<>
	const float64_t tolerance<float64_t> = 1e-6;

	// scale abs(f) by tolerance, unless f is very small
	template<typename T>
	static T scaled_tolerance(T f) {
		return tolerance<T>*std::max(static_cast<T>(1.0), std::abs(f));
	}


	// search the line by fitting a local quadratic model, taking a step, and then surfing the slope
	template<typename T>
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
		T fxp = std::numeric_limits<T>::infinity();
		if (xp <= xmax) {
			fxp = f(xp);
		}
		T fxm = std::numeric_limits<T>::infinity();
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
		T xstar = 0;
		if ((shape < -shape_epsilon) || std::isnan(shape) || std::isinf(shape)) {

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

		xstar = std::clamp(xstar, xmin, xmax);
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
					if (std::isnan(fxmin)) {
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
					if (std::isnan(fxmax)) {
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
	using LineSearchFunction = T (*)(Dofs<T> &, int, T, T &);

	template<typename T>
	static void minimize_ccd(Dofs<T> & dofs, DofValues<T> & here) {

		// get the current objective function value
		here.f = dofs.eval_efunc(here.x);

		// ccd is pretty simple actually
		// just do a line search along each dimension until we stop improving
		// we deal with cycles by just capping the number of iterations

		const int max_iterations = 30;
		const T convergence_threshold = 0.001;

		// init line search
		auto line_search_states = Array<LineSearchState<T>>(dofs.get_size());
		LineSearchFunction<T> line_search = line_search_surf;

		DofValues<T> next(dofs.get_size());
		for (int iter=0; iter<max_iterations; iter++) {

			// update all the dofs using line search
			next.set(here);
			for (int d=0; d<dofs.get_size(); d++) {

				Dof<T> & dof = dofs[d];
				LineSearchState<T> & state = line_search_states[d];

				// get the step size, try to make it adaptive (based on historical steps if possible; else on step #)
				T step;
				if (std::abs(state.last_step) > tolerance<T> && std::abs(state.first_step) > tolerance<T>) {
					step = dof.initial_step_size*std::abs(state.last_step/state.first_step);
				} else {
					step = dof.initial_step_size/std::pow(iter + 1, 3);
				}

				// get the next x value for this dof
				next.x[d] = line_search(dofs, d, next.x[d], step);

				if (iter == 0) {
					state.first_step = step;
				}
				state.last_step = step;
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
