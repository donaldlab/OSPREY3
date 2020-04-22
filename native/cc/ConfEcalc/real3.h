
#ifndef CONFECALC_REAL3_H
#define CONFECALC_REAL3_H


namespace osprey {

	template<typename T>
	class alignas(8) Real3 {

		public:
			T x;
			T y;
			T z;
			// 4 bytes pad, if T = float32_t

			Real3() : x(0.0), y(0.0), z(0.0) {}

			Real3(T _x, T _y, T _z) {
				x = _x;
				y = _y;
				z = _z;
			}
	};
	ASSERT_JAVA_COMPATIBLE_REALS(Real3, 16, 24);
	ASSERT_COPYABLE_REALS(Real3);

	template<typename T>
	std::ostream & operator << (std::ostream & out, const Real3<T> & v) {
		auto w = fmt::real_width(out);
		out << "("
			<< std::setw(w) << v.x << ", "
			<< std::setw(w) << v.y << ", "
			<< std::setw(w) << v.z << ")";
		return out;
	}
}

#endif //CONFECALC_REAL3_H
