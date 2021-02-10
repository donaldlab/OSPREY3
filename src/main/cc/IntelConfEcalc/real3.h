
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

			inline Real3() : x(0.0), y(0.0), z(0.0) {}

			inline Real3(T _x, T _y, T _z) {
				x = _x;
				y = _y;
				z = _z;
			}

			inline void operator += (const Real3<T> & v) {
				x += v.x;
				y += v.y;
				z += v.z;
			}

			inline void operator -= (const Real3<T> & v) {
				x -= v.x;
				y -= v.y;
				z -= v.z;
			}

			inline T len_sq() const {
				return x*x + y*y + z*z;
			}
			inline T len() const {
				return std::sqrt(len_sq());
			}

			inline void normalize() {
				T invlen = 1.0/len();
				x *= invlen;
				y *= invlen;
				z *= invlen;
			}

			inline void negate() {
				x = -x;
				y = -y;
				z = -z;
			}

			inline T dot(const Real3<T> v) const {
				return x*v.x + y*v.y + z*v.z;
			}

			inline bool isnan() const {
				return std::isnan(x) || std::isnan(y) || std::isnan(z);
			}
	};
	ASSERT_JAVA_COMPATIBLE_REALS(Real3, 16, 24);
	ASSERT_COPYABLE_REALS(Real3);

	template<typename T>
	inline Real3<T> operator - (const Real3<T> & v) {
		return {
			-v.x,
			-v.y,
			-v.z
		};
	}

	template<typename T>
	inline Real3<T> operator + (const Real3<T> & a, const Real3<T> & b) {
		return {
			a.x + b.x,
			a.y + b.y,
			a.z + b.z
		};
	}

	template<typename T>
	inline Real3<T> operator - (const Real3<T> & a, const Real3<T> & b) {
		return {
			a.x - b.x,
			a.y - b.y,
			a.z - b.z
		};
	}

	template<typename T>
	inline Real3<T> cross(const Real3<T> & a, const Real3<T> & b) {
		return {
			a.y*b.z - a.z*b.y,
			a.z*b.x - a.x*b.z,
			a.x*b.y - a.y*b.x
		};
	}

	template<typename T>
	std::ostream & operator << (std::ostream & out, const Real3<T> & v) {
		auto w = fmt::real_width(out);
		out << "("
			<< std::setw(w) << v.x << ", "
			<< std::setw(w) << v.y << ", "
			<< std::setw(w) << v.z << ")";
		return out;
	}

	template<typename T>
	static T distance_sq(const Real3<T> & a, const Real3<T> & b) {
		T dx = a.x - b.x;
		T dy = a.y - b.y;
		T dz = a.z - b.z;
		return dx*dx + dy*dy + dz*dz;
	}

	template<typename T>
	static T distance(const Real3<T> & a, const Real3<T> & b) {
		return std::sqrt(distance_sq(a, b));
	}

}

#endif //CONFECALC_REAL3_H
