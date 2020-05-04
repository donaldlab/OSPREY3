
#ifndef CONFECALC_ENERGY_H
#define CONFECALC_ENERGY_H


namespace osprey {

	template<typename T>
	struct PosInter {
		int32_t posi1;
		int32_t posi2;
		T weight;
		T offset;
	};
	ASSERT_JAVA_COMPATIBLE_REALS(PosInter, 16, 24);

	template<typename T>
	using EnergyFunction = T (*)(Assignment<T> &, const Array<PosInter<T>> &, cg::thread_group, T[]);
}


#endif //CONFECALC_ENERGY_H
