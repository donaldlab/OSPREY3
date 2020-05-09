
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

	// Amber atom pair structs

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

	class AmberAtomPairPtrs {
		public:
			const void * a;
			const void * b;

			__device__
			AmberAtomPairPtrs(): a(nullptr), b(nullptr) {}
			AmberAtomPairPtrs(const AmberAtomPairPtrs & other) = default;
			~AmberAtomPairPtrs() = default;

			template<typename T>
			__device__
			inline void set(const AtomPairs * atom_pairs);
	};

	template<>
	__device__
	inline void AmberAtomPairPtrs::set<float32_t>(const AtomPairs * atom_pairs) {
		auto p = reinterpret_cast<const int8_t *>(atom_pairs + 1);
		a = reinterpret_cast<const void *>(p);
		p += cuda::pad_to_alignment<alignof(AtomPairAmberF32b)>(sizeof(AtomPairAmberF32a)*atom_pairs->num_amber);
		b = reinterpret_cast<const void *>(p);
	}

	template<>
	__device__
	inline void AmberAtomPairPtrs::set<float64_t>(const AtomPairs * atom_pairs) {
		auto p = reinterpret_cast<const int8_t *>(atom_pairs + 1);
		a = reinterpret_cast<const void *>(p);
		p += cuda::pad_to_alignment<alignof(AtomPairAmberF64b)>(sizeof(AtomPairAmberF64a)*atom_pairs->num_amber);
		b = reinterpret_cast<const void *>(p);
	}

	template<typename T>
	__device__
	inline static T calc_amber(const AmberAtomPairPtrs & pair_ptrs, int pairi, bool distance_dependent_dielectric, const Real3<T> * atoms) {}

	template<>
	__device__
	inline float32_t calc_amber<float32_t>(const AmberAtomPairPtrs & pair_ptrs, int pairi, bool distance_dependent_dielectric, const Real3<float32_t> * atoms) {

		// TODO: optimize/cache pointer offset calculations?

		// read the atom indices for this atom pair
		const AtomPairAmberF32a paira = reinterpret_cast<const AtomPairAmberF32a *>(pair_ptrs.a)[pairi];

		// calculate the radius
		const Real3<float32_t> atom1 = atoms[paira.atomi1];
		const Real3<float32_t> atom2 = atoms[paira.atomi2];
		float32_t oor2 = rcp_intr<float32_t>(distance_sq<float32_t>(atom1, atom2));
		assert (!isnan<float32_t>(oor2));
		assert (!isinf<float32_t>(oor2));

		// read the amber params for this atom pair
		const AtomPairAmberF32b pairb = reinterpret_cast<const AtomPairAmberF32b *>(pair_ptrs.b)[pairi];

		// calculate the electrostatics energy
		float32_t energy = pairb.esQ;
		if (distance_dependent_dielectric) {
			energy *= oor2;
		} else {
			energy *= sqrt_intr<float32_t>(oor2);
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
	inline float64_t calc_amber<float64_t>(const AmberAtomPairPtrs & pair_ptrs, int pairi, bool distance_dependent_dielectric, const Real3<float64_t> * atoms) {

		// read the atom indices for this atom pair
		const AtomPairAmberF64a paira = reinterpret_cast<const AtomPairAmberF64a *>(pair_ptrs.a)[pairi];

		// calculate the radius
		const Real3<float64_t> atom1 = atoms[paira.atomi1];
		const Real3<float64_t> atom2 = atoms[paira.atomi2];
		float64_t oor2 = rcp_intr<float64_t>(distance_sq<float64_t>(atom1, atom2));
		assert (!isnan<float32_t>(oor2));
		assert (!isinf<float32_t>(oor2));

		// calculate the electrostatics energy
		float32_t energy = paira.esQ;
		if (distance_dependent_dielectric) {
			energy *= oor2;
		} else {
			energy *= sqrt_intr<float64_t>(oor2);
		}

		// read the amber params for this atom pair
		const AtomPairAmberF64b pairb = reinterpret_cast<const AtomPairAmberF64b *>(pair_ptrs.b)[pairi];

		// calculate the van der Waals energy
		float64_t oor6 = oor2*oor2*oor2;
		energy -= pairb.vdwB*oor6;
		float64_t oor12 = oor6*oor6;
		energy += pairb.vdwA*oor12;

		return energy;
	}


	// EEF1 atom pair structs

	struct alignas(8) AtomPairEef1F32a {
		int32_t atomi1;
		int32_t atomi2;
	};
	ASSERT_JAVA_COMPATIBLE(AtomPairEef1F32a, 8);

	struct alignas(16) AtomPairEef1F32b {
		float32_t vdwRadius;
		float32_t oolambda; // 1/lambda
		float32_t alpha;
		// 4 byte pad
	};
	ASSERT_JAVA_COMPATIBLE(AtomPairEef1F32b, 16);

	struct alignas(8) AtomPairEef1F64a {
		int32_t atomi1;
		int32_t atomi2;
	};
	ASSERT_JAVA_COMPATIBLE(AtomPairEef1F64a, 8);

	struct alignas(16) AtomPairEef1F64b {
		float64_t vdwRadius;
		float64_t oolambda; // 1/lambda
	};
	ASSERT_JAVA_COMPATIBLE(AtomPairEef1F64b, 16);

	struct alignas(16) AtomPairEef1F64c {
		float64_t alpha1;
		float64_t alpha2;
	};
	ASSERT_JAVA_COMPATIBLE(AtomPairEef1F64c, 16);

	class Eef1AtomPairPtrs {
		public:
			const void * a;
			const void * b;
			const void * c;
			const void * d;

			__device__
			Eef1AtomPairPtrs(): a(nullptr), b(nullptr), c(nullptr), d(nullptr) {}
			Eef1AtomPairPtrs(const Eef1AtomPairPtrs & other) = default;
			~Eef1AtomPairPtrs() = default;

			template<typename T>
			__device__
			inline void set(const AtomPairs * atom_pairs);
	};

	template<>
	__device__
	inline void Eef1AtomPairPtrs::set<float32_t>(const AtomPairs * atom_pairs) {
		auto p = reinterpret_cast<const int8_t *>(atom_pairs + 1);

		// skip all the amber atom pairs
		p += cuda::pad_to_alignment<alignof(AtomPairAmberF32b)>(sizeof(AtomPairAmberF32a)*atom_pairs->num_amber);
		p += sizeof(AtomPairAmberF32b)*atom_pairs->num_amber;

		a = reinterpret_cast<const void *>(p);
		p += cuda::pad_to_alignment<alignof(AtomPairEef1F32b)>(sizeof(AtomPairEef1F32a)*atom_pairs->num_eef1);
		b = reinterpret_cast<const void *>(p);
		p += sizeof(AtomPairEef1F32b)*atom_pairs->num_eef1;
		c = reinterpret_cast<const void *>(p);
		d = nullptr;
	}

	template<>
	__device__
	inline void Eef1AtomPairPtrs::set<float64_t>(const AtomPairs * atom_pairs) {
		auto p = reinterpret_cast<const int8_t *>(atom_pairs + 1);

		// skip all the amber atom pairs
		p += cuda::pad_to_alignment<alignof(AtomPairAmberF64b)>(sizeof(AtomPairAmberF64a)*atom_pairs->num_amber);
		p += sizeof(AtomPairAmberF64b)*atom_pairs->num_amber;

 		a = reinterpret_cast<const void *>(p);
		p += cuda::pad_to_alignment<alignof(AtomPairEef1F64b)>(sizeof(AtomPairEef1F64a)*atom_pairs->num_eef1);
		b = reinterpret_cast<const void *>(p);
		p += cuda::pad_to_alignment<alignof(AtomPairEef1F64b)>(sizeof(AtomPairEef1F64b)*atom_pairs->num_eef1);
		c = reinterpret_cast<const void *>(p);
		p += cuda::pad_to_alignment<alignof(AtomPairEef1F64c)>(sizeof(AtomPairEef1F64b)*atom_pairs->num_eef1);
		d = reinterpret_cast<const void *>(p);
	}

	template<typename T>
	__device__
	inline static T calc_eef1(const Eef1AtomPairPtrs & pair_ptrs, int pairi, const Real3<T> * atoms) {}

	template<>
	__device__
	inline float32_t calc_eef1<float32_t>(const Eef1AtomPairPtrs & pair_ptrs, int pairi, const Real3<float32_t> * atoms) {

		// read the atom indices for this atom pair
		const AtomPairEef1F32a paira = reinterpret_cast<const AtomPairEef1F32a *>(pair_ptrs.a)[pairi];

		// calculate the radius
		const Real3<float32_t> atom1 = atoms[paira.atomi1];
		const Real3<float32_t> atom2 = atoms[paira.atomi2];
		float32_t r2 = distance_sq<float32_t>(atom1, atom2);
		assert (!isnan<float32_t>(r2));
		assert (r2 > 0.0);

		// cutoff when r > 9
		if (r2 > 81.0) {
			return 0.0;
		}

		float32_t r = sqrt_intr<float32_t>(r2);
		float32_t oor2 = rcp_intr(r2);

		const AtomPairEef1F32b pairb = reinterpret_cast<const AtomPairEef1F32b *>(pair_ptrs.b)[pairi];
		float32_t Xij = (r - pairb.vdwRadius)*pairb.oolambda;
		float32_t energy = -pairb.alpha*exp_intr<float32_t>(-Xij*Xij)*oor2;

		const AtomPairEef1F32b pairc = reinterpret_cast<const AtomPairEef1F32b *>(pair_ptrs.c)[pairi];
		float32_t Xji = (r - pairc.vdwRadius)*pairc.oolambda;
		energy -= pairc.alpha*exp_intr<float32_t>(-Xji*Xji)*oor2;

		return energy;
	}

	template<>
	__device__
	inline float64_t calc_eef1<float64_t>(const Eef1AtomPairPtrs & pair_ptrs, int pairi, const Real3<float64_t> * atoms) {

		// read the atom indices for this atom pair
		const AtomPairEef1F64a paira = reinterpret_cast<const AtomPairEef1F64a *>(pair_ptrs.a)[pairi];

		// calculate the radius
		const Real3<float64_t> atom1 = atoms[paira.atomi1];
		const Real3<float64_t> atom2 = atoms[paira.atomi2];
		float64_t r2 = distance_sq<float64_t>(atom1, atom2);
		assert (!isnan<float64_t>(r2));
		assert (r2 > 0.0);

		// cutoff when r > 9
		if (r2 > 81.0) {
			return 0.0;
		}

		float64_t r = sqrt_intr<float64_t>(r2);

		const AtomPairEef1F64b pairb = reinterpret_cast<const AtomPairEef1F64b *>(pair_ptrs.b)[pairi];
		float64_t Xij = (r - pairb.vdwRadius)*pairb.oolambda;

		const AtomPairEef1F64b pairc = reinterpret_cast<const AtomPairEef1F64b *>(pair_ptrs.c)[pairi];
		float64_t Xji = (r - pairc.vdwRadius)*pairc.oolambda;

		const AtomPairEef1F64c paird = reinterpret_cast<const AtomPairEef1F64c *>(pair_ptrs.d)[pairi];
		return -(paird.alpha1*exp_intr<float64_t>(-Xij*Xij) + paird.alpha2*exp_intr<float64_t>(-Xji*Xji))*rcp_intr(r2);
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

		const Real3<T> * atoms = assignment.atoms.items();

		int64_t start_global_pairi = 0;
		int64_t stop_global_pairi = 0;
		int interi = -1;
		T inter_weight = 0.0;

		Eef1AtomPairPtrs eef1_pairs;

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

				// update the pair index range
				auto atom_pairs = reinterpret_cast<const AtomPairs *>(assignment.get_atom_pairs(inter.posi1, inter.posi2));
				start_global_pairi = stop_global_pairi;
				stop_global_pairi += atom_pairs->num_eef1;

				// continue early, if possible
				if (global_pairi >= stop_global_pairi) {
					continue;
				}

				eef1_pairs.set<T>(atom_pairs);
				inter_weight = inter.weight;
			}

			int64_t pairi = global_pairi - start_global_pairi;
			energy += calc_eef1<T>(eef1_pairs, pairi, atoms)*inter_weight;
		}

		finish_eef1:

		// reset the counters
		start_global_pairi = 0;
		stop_global_pairi = 0;
		interi = -1;
		inter_weight = 0.0;

		const Params & params = *reinterpret_cast<const Params *>(assignment.conf_space.get_params());
		bool distance_dependent_dielectric = params.distance_dependent_dielectric;
		AmberAtomPairPtrs amber_pairs;

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

				// update the pair index range
				auto atom_pairs = reinterpret_cast<const AtomPairs *>(assignment.get_atom_pairs(inter.posi1, inter.posi2));
				start_global_pairi = stop_global_pairi;
				stop_global_pairi += atom_pairs->num_amber;

				// continue early, if possible
				if (global_pairi >= stop_global_pairi) {
					continue;
				}

				amber_pairs.set<T>(atom_pairs);
				inter_weight = inter.weight;
			}

			int64_t pairi = global_pairi - start_global_pairi;
			energy += calc_amber<T>(amber_pairs, pairi, distance_dependent_dielectric, atoms)*inter_weight;
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
