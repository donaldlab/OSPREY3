
#ifndef CONFECALC_ENERGY_AMBEREEF1_H
#define CONFECALC_ENERGY_AMBEREEF1_H


namespace osprey { namespace ambereef1 {

	struct alignas(16) Params {
		bool distance_dependent_dielectric;
		// 15 bytes pad
	};
	ASSERT_JAVA_COMPATIBLE(Params, 16);

	struct alignas(16) AtomPairs {
		int32_t num_amber;
		int32_t num_eef1;
	};
	ASSERT_JAVA_COMPATIBLE(AtomPairs, 16);

	// NOTE: try to aim for atom pair structs to 16 bytes,
	// (less is ok, but not more!)
	// so we can efficiently coalesce LDs in a warp

	struct alignas(8) AtomPairAmberF32a {
		int32_t atomi1;
		int32_t atomi2;
	};
	ASSERT_JAVA_COMPATIBLE(AtomPairAmberF32a, 8);

	struct alignas(16) AtomPairAmberF32b {
		float32_t esQ;
		float32_t vdwA;
		float32_t vdwB;
		// 4 bytes pad
	};
	ASSERT_JAVA_COMPATIBLE(AtomPairAmberF32b, 16);

	struct alignas(16) AtomPairAmberF64a {
		int32_t atomi1;
		int32_t atomi2;
		float64_t esQ;
	};
	ASSERT_JAVA_COMPATIBLE(AtomPairAmberF64a, 16);

	struct alignas(16) AtomPairAmberF64b {
		float64_t vdwA;
		float64_t vdwB;
	};
	ASSERT_JAVA_COMPATIBLE(AtomPairAmberF64b, 16);

	// TODO: optimize EEF1 memory layout
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

		__device__
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
	__device__
	inline static T calc_amber(const AtomPairs & atom_pairs, int pairi, bool distance_dependent_dielectric, const Array<Real3<T>> & atoms) {}

	template<>
	__device__
	inline float32_t calc_amber<float32_t>(const AtomPairs & atom_pairs, int pairi, bool distance_dependent_dielectric, const Array<Real3<float32_t>> & atoms) {

		// TODO: optimize/cache pointer offset calculations?

		// read the atom indices for this atom pair
		auto p = reinterpret_cast<const int8_t *>(&atom_pairs + 1);
		const AtomPairAmberF32a paira = *reinterpret_cast<const AtomPairAmberF32a *>(p
			+ sizeof(AtomPairAmberF32a)*pairi);

		// calculate the radius
		const Real3<float32_t> atom1 = atoms[paira.atomi1];
		const Real3<float32_t> atom2 = atoms[paira.atomi2];
		float32_t oor2 = 1.0f/distance_sq<float32_t>(atom1, atom2);
		assert (!isnan<float32_t>(oor2));
		assert (!isinf<float32_t>(oor2));

		// read the amber params for this atom pair
		const AtomPairAmberF32b pairb = *reinterpret_cast<const AtomPairAmberF32b *>(p
			+ cuda::pad_to_alignment<alignof(AtomPairAmberF32b)>(sizeof(AtomPairAmberF32a)*atom_pairs.num_amber)
			+ sizeof(AtomPairAmberF32b)*pairi);

		// calculate the electrostatics energy
		float32_t energy = pairb.esQ;
		if (distance_dependent_dielectric) {
			energy *= oor2;
		} else {
			energy *= std::sqrt(oor2);
		}

		// calculate the van der Waals energy
		float32_t oor6 = oor2*oor2*oor2;
		energy -= pairb.vdwB*oor6;
		float32_t oor12 = oor6*oor6;
		energy += pairb.vdwA*oor12;

		return energy;
	}

	template<>
	__device__
	inline float64_t calc_amber<float64_t>(const AtomPairs & atom_pairs, int pairi, bool distance_dependent_dielectric, const Array<Real3<float64_t>> & atoms) {

		// read the atom indices for this atom pair
		auto p = reinterpret_cast<const int8_t *>(&atom_pairs + 1);
		const AtomPairAmberF64a paira = *reinterpret_cast<const AtomPairAmberF64a *>(p
			+ sizeof(AtomPairAmberF64a)*pairi);

		// calculate the radius
		const Real3<float64_t> atom1 = atoms[paira.atomi1];
		const Real3<float64_t> atom2 = atoms[paira.atomi2];
		float64_t oor2 = 1.0/distance_sq<float64_t>(atom1, atom2);
		assert (!isnan<float32_t>(oor2));
		assert (!isinf<float32_t>(oor2));

		// calculate the electrostatics energy
		float32_t energy = paira.esQ;
		if (distance_dependent_dielectric) {
			energy *= oor2;
		} else {
			energy *= std::sqrt(oor2);
		}

		// read the amber params for this atom pair
		const AtomPairAmberF64b pairb = *reinterpret_cast<const AtomPairAmberF64b *>(p
			+ sizeof(AtomPairAmberF64a)*atom_pairs.num_amber
			+ sizeof(AtomPairAmberF64b)*pairi);

		// calculate the van der Waals energy
		float64_t oor6 = oor2*oor2*oor2;
		energy -= pairb.vdwB*oor6;
		float64_t oor12 = oor6*oor6;
		energy += pairb.vdwA*oor12;

		return energy;
	}

	template<typename T>
	__device__
	inline static T calc_eef1(const AtomPairs & atom_pairs, int pairi, const Array<Real3<T>> & atoms) {}

	template<>
	__device__
	inline float32_t calc_eef1<float32_t>(const AtomPairs & atom_pairs, int pairi, const Array<Real3<float32_t>> & atoms) {

		// TODO: optimize memory accesses

		// eef1 atom pairs are placed right after the AtomPairs struct and the amber atom pairs
		auto p = reinterpret_cast<const int8_t *>(&atom_pairs + 1);
		auto pairs_eef1 = reinterpret_cast<const AtomPairEef1<float32_t> *>(p
			+ cuda::pad_to_alignment<alignof(AtomPairAmberF32b)>(sizeof(AtomPairAmberF32a)*atom_pairs.num_amber)
			+ sizeof(AtomPairAmberF32b)*atom_pairs.num_amber
		);

		// calculate energy for the next atom pair
		auto pair = pairs_eef1[pairi];
		const Real3<float32_t> & atom1 = atoms[pair.atomi1];
		const Real3<float32_t> & atom2 = atoms[pair.atomi2];
		float32_t r2 = distance_sq<float32_t>(atom1, atom2);
		assert (!isnan<float32_t>(r2));
		assert (r2 > 0.0);
		float32_t r = std::sqrt(r2);

		return pair.calc(r, r2);
	}

	template<>
	__device__
	inline float64_t calc_eef1<float64_t>(const AtomPairs & atom_pairs, int pairi, const Array<Real3<float64_t>> & atoms) {

		// TODO: optimize memory accesses

		// eef1 atom pairs are placed right after the AtomPairs struct and the amber atom pairs
		auto p = reinterpret_cast<const int8_t *>(&atom_pairs + 1);
		const AtomPairEef1<float64_t> * pairs_eef1 = reinterpret_cast<const AtomPairEef1<float64_t> *>(p
			+ cuda::pad_to_alignment<alignof(AtomPairAmberF64b)>(sizeof(AtomPairAmberF64a)*atom_pairs.num_amber)
			+ sizeof(AtomPairAmberF64b)*atom_pairs.num_amber
		);

		// calculate energy for the next atom pair
		auto pair = pairs_eef1[pairi];
		const Real3<float64_t> & atom1 = atoms[pair.atomi1];
		const Real3<float64_t> & atom2 = atoms[pair.atomi2];
		float64_t r2 = distance_sq<float64_t>(atom1, atom2);
		assert (!isnan<float64_t>(r2));
		assert (r2 > 0.0);
		float64_t r = std::sqrt(r2);
		return pair.calc(r, r2);
	}

	template<typename T>
	__device__
	T calc_energy(Assignment<T> & assignment, const Array<PosInter<T>> & inters, T thread_energy[]) {

		T energy = 0.0;

		// sum all the interaction energies first
		for (int i=threadIdx.x; i<inters.get_size(); i+=blockDim.x) {
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

			// apply the weight and offset
			energy += inter.weight*(inter_energy + inter.offset);
		}

		const AtomPairs * atom_pairs = nullptr;
		int64_t start_global_pairi = 0;
		int64_t stop_global_pairi = 0;
		int interi = -1;
		T inter_weight = 0.0;

		// sum all the EEF1 energies
		for (int64_t global_pairi = threadIdx.x; ; global_pairi += blockDim.x) {

			// advance to the next group of atom pairs, if needed
			while (global_pairi >= stop_global_pairi) {

				// advance to the next interaction, if any
				interi++;
				if (interi >= inters.get_size()) {
					goto finish_eef1; // break out of both loops
				}
				const PosInter<T> & inter = inters[interi];
				inter_weight = inter.weight;

				// update the pair index range
				atom_pairs = reinterpret_cast<const AtomPairs *>(assignment.get_atom_pairs(inter.posi1, inter.posi2));
				start_global_pairi = stop_global_pairi;
				stop_global_pairi += atom_pairs->num_eef1;
			}

			int64_t pairi = global_pairi - start_global_pairi;
			energy += calc_eef1<T>(*atom_pairs, pairi, assignment.atoms)*inter_weight;
		}

		finish_eef1:

		// reset the counters
		atom_pairs = nullptr;
		start_global_pairi = 0;
		stop_global_pairi = 0;
		interi = -1;
		inter_weight = 0.0;

		const Params & params = *reinterpret_cast<const Params *>(assignment.conf_space.get_params());
		bool distance_dependent_dielectric = params.distance_dependent_dielectric;

		// sum all the Amber energies
		for (int64_t global_pairi = threadIdx.x; ; global_pairi += blockDim.x) {

			// advance to the next group of atom pairs, if needed
			while (global_pairi >= stop_global_pairi) {

				// advance to the next interaction, if any
				interi++;
				if (interi >= inters.get_size()) {
					goto finish_amber; // break out of both loops
				}
				PosInter<T> inter = inters[interi];
				inter_weight = inter.weight;

				// update the pair index range
				atom_pairs = reinterpret_cast<const AtomPairs *>(assignment.get_atom_pairs(inter.posi1, inter.posi2));
				start_global_pairi = stop_global_pairi;
				stop_global_pairi += atom_pairs->num_amber;
			}

			int64_t pairi = global_pairi - start_global_pairi;
			energy += calc_amber<T>(*atom_pairs, pairi, distance_dependent_dielectric, assignment.atoms)*inter_weight;
		}

		finish_amber:

		// do a sum reduction across all threads
		thread_energy[threadIdx.x] = energy;
		__syncthreads();

		for (uint offset=1; offset<blockDim.x; offset<<=1) {

			// sum this level of the reduction tree
			int mask = (offset << 1) - 1;
			if ((threadIdx.x & mask) == 0) {
				int pos = threadIdx.x + offset;
				if (pos < blockDim.x) {
					thread_energy[threadIdx.x] += thread_energy[pos];
				}
			}
			__syncthreads();
		}

		return thread_energy[0];
	}
}}


#endif //CONFECALC_ENERGY_AMBEREEF1_H
