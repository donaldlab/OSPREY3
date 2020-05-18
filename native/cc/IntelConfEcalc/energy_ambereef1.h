
#ifndef CONFECALC_ENERGY_AMBEREEF1_H
#define CONFECALC_ENERGY_AMBEREEF1_H


namespace osprey::ambereef1 {

	struct alignas(8) Params {
		bool distance_dependent_dielectric;
		// 7 bytes pad
	};
	ASSERT_JAVA_COMPATIBLE(Params, 8);

	struct AtomPairs {
		int32_t num_amber;
		int32_t num_eef1;
	};
	ASSERT_JAVA_COMPATIBLE(AtomPairs, 8);

	template<typename T>
	struct alignas(8) AtomPairAmber {
		int32_t atomi1;
		int32_t atomi2;
		T esQ;
		T vdwA;
		T vdwB;
		// 4 bytes pad, if T = float32_t

		inline T calc(T r, T r2, bool distance_dependent_dielectric) const {

			// just in case ...
			assert(r >= 0.0);
			assert(r2 >= 0.0);

			// TODO: optimize out division?

			// calculate the electrostatics energy
			T es;
			if (distance_dependent_dielectric) {
				es = esQ/r2;
			} else {
				es = esQ/r;
			}

			// calculate the van der Waals energy
			T r6 = r2*r2*r2;
			T r12 = r6*r6;
			T vdw = vdwA/r12 - vdwB/r6;

			return es + vdw;
		}
	};
	ASSERT_JAVA_COMPATIBLE_REALS(AtomPairAmber, 24, 32);

	template<typename T>
	struct alignas(8) AtomPairEef1 {
		int32_t atomi1;
		int32_t atomi2;
		T vdwRadius1;
		T lambda1;
		T vdwRadius2;
		T lambda2;
		T alpha1;
		T alpha2;

		inline T calc(T r, T r2) const {

			// just in case ...
			assert(r >= 0.0);
			assert(r2 >= 0.0);
			assert(lambda1 >= 0.0);
			assert(lambda2 >= 0.0);

			const T cutoff = 9.0;
			if (r > cutoff) {
				return 0.0;
			}

			// TODO: optimize out division?

			T Xij = (r - vdwRadius1)/lambda1;
			T Xji = (r - vdwRadius2)/lambda2;
			return -(alpha1*exp(-Xij*Xij) + alpha2*exp(-Xji*Xji))/r2;
		}
	};
	ASSERT_JAVA_COMPATIBLE_REALS(AtomPairEef1, 32, 56);

	template<typename T>
	static T calc(const Array<Real3<T>> & atoms, const Params & params, const AtomPairs & pairs) {

		T energy = 0.0;

		// add the amber interactions
		auto pair_amber = reinterpret_cast<const AtomPairAmber<T> *>(&pairs + 1);
		for (int i=0; i<pairs.num_amber; i++) {
			Real3<T> atom1 = atoms[pair_amber->atomi1];
			Real3<T> atom2 = atoms[pair_amber->atomi2];
			T r2 = distance_sq(atom1, atom2);
			assert (!std::isnan(r2));
			T r = std::sqrt(r2);
			energy += pair_amber->calc(r, r2, params.distance_dependent_dielectric);
			pair_amber += 1;
		}

		// add the eef1 interactions
		auto pair_eef1 = reinterpret_cast<const AtomPairEef1<T> *>(pair_amber);
		for (int i=0; i<pairs.num_eef1; i++) {
			Real3<T> atom1 = atoms[pair_eef1->atomi1];
			Real3<T> atom2 = atoms[pair_eef1->atomi2];
			T r2 = distance_sq(atom1, atom2);
			assert (!std::isnan(r2));
			T r = std::sqrt(r2);
			energy += pair_eef1->calc(r, r2);
			pair_eef1 += 1;
		}

		return energy;
	}

	template<typename T>
	T calc_energy(Assignment<T> & assignment, const Array<PosInter<T>> & inters) {

		const Params & params = *reinterpret_cast<const Params *>(assignment.conf_space.get_params());

		// sum all the interaction energies
		T energy = 0.0;
		for (int i=0; i<inters.get_size(); i++) {
			PosInter<T> inter = inters[i];

			T inter_energy = 0.0;

			// add extra static/internal energy if needed
			if (inter.posi1 == inter.posi2) {
				if (inter.posi1 == StaticPos) {
					// TODO: add a switch to turn this off?
					inter_energy += assignment.conf_space.static_energy;
				} else {
					inter_energy += assignment.get_conf_energy(inter.posi1);
				}
			}

			// add the energy for the atom pairs
			const AtomPairs & atom_pairs = *reinterpret_cast<const AtomPairs *>(assignment.get_atom_pairs(inter.posi1, inter.posi2));
			inter_energy += calc(assignment.atoms, params, atom_pairs);

			// apply weight and offset
			energy += inter.weight*(inter_energy + inter.offset);
		}

		return energy;
	}
}


#endif //CONFECALC_ENERGY_AMBEREEF1_H
