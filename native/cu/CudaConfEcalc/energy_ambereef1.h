
#ifndef CONFECALC_ENERGY_AMBEREEF1_H
#define CONFECALC_ENERGY_AMBEREEF1_H


namespace osprey { namespace ambereef1 {

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

		__device__
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

		int64_t global_pairi = threadIdx.x;
		int64_t start_global_pairi = 0;
		int64_t stop_global_pairi = 0;
		int interi = -1;
		T inter_weight = 0.0;

		// sum all the EEF1 energies
		const AtomPairEef1<T> * pairs_eef1 = nullptr;
		while (true) {

			// advance to the next group of atom pairs, if needed
			while (global_pairi >= stop_global_pairi) {

				// advance to the next interaction, if any
				interi++;
				if (interi >= inters.get_size()) {
					goto finish_eef1; // break out of both loops
				}
				PosInter<T> inter = inters[interi];
				inter_weight = inter.weight;

				// update the pair index range
				const AtomPairs & atom_pairs = *reinterpret_cast<const AtomPairs *>(assignment.get_atom_pairs(inter.posi1, inter.posi2));
				start_global_pairi = stop_global_pairi;
				stop_global_pairi += atom_pairs.num_eef1;

				// get the EEF1 atom pairs
				auto p = reinterpret_cast<const int8_t *>(&atom_pairs + 1);
				pairs_eef1 = reinterpret_cast<const AtomPairEef1<T> *>(p + sizeof(AtomPairAmber<T>)*atom_pairs.num_amber);
			}

			// calculate energy for the next atom pair
			auto pair = pairs_eef1[global_pairi - start_global_pairi];
			T r2 = distance_sq<T>(assignment.atoms[pair.atomi1], assignment.atoms[pair.atomi2]);
			assert (!isnan<T>(r2));
			T r = std::sqrt(r2);
			energy += pair.calc(r, r2)*inter_weight;

			// advance to the next atom pair
			global_pairi += blockDim.x;
		}

		finish_eef1:

		// reset the counters
		global_pairi = threadIdx.x;
		start_global_pairi = 0;
		stop_global_pairi = 0;
		interi = -1;
		inter_weight = 0.0;

		const Params & params = *reinterpret_cast<const Params *>(assignment.conf_space.get_params());

		// sum all the Amber energies
		const AtomPairAmber<T> * pairs_amber = nullptr;
		while (true) {

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
				const AtomPairs & atom_pairs = *reinterpret_cast<const AtomPairs *>(assignment.get_atom_pairs(inter.posi1, inter.posi2));
				start_global_pairi = stop_global_pairi;
				stop_global_pairi += atom_pairs.num_amber;

				// get the Amber atom pairs
				pairs_amber = reinterpret_cast<const AtomPairAmber<T> *>(&atom_pairs + 1);
			}

			// TODO: OPTIMIZATION: optimize memory access patterns here

			// calculate energy for the next atom pair
			auto pair = pairs_amber[global_pairi - start_global_pairi];
			T r2 = distance_sq<T>(assignment.atoms[pair.atomi1], assignment.atoms[pair.atomi2]);
			assert (!isnan<T>(r2));
			T r = std::sqrt(r2);
			energy += pair.calc(r, r2, params.distance_dependent_dielectric)*inter_weight;

			// advance to the next atom pair
			global_pairi += blockDim.x;
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
