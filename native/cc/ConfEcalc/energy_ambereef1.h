
#ifndef CONFECALC_ENERGY_AMBEREEF1_H
#define CONFECALC_ENERGY_AMBEREEF1_H


namespace osprey {

	struct alignas(8) ParamsAmberEef1 {
		bool distance_dependent_dielectric;
		// 7 bytes pad
	};
	ASSERT_JAVA_COMPATIBLE(ParamsAmberEef1, 8);

	template<typename T>
	T calc_amber_eef1(const osprey::ConfSpace<T> & conf_space, const int32_t conf[], const osprey::PosInter<T> inters[], int64_t inters_size) {

		osprey::Assignment<T> assignment(conf_space, conf);

		// TODO

		return 0.3;
	}
}


#endif //CONFECALC_ENERGY_AMBEREEF1_H
