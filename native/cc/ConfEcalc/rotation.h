
#ifndef CONFECALC_ROTATION_H
#define CONFECALC_ROTATION_H


namespace osprey {

	template<typename T>
	const T Pi = 3.14159265358979323846;

	template<typename T>
	const T TwoPi = Pi<T>*2;

	template<typename T>
	const T HalfPi = Pi<T>/2;

	template<typename T>
	T normalize_mpi_pi(T radians) {

		assert (std::isfinite(radians));

		while (radians <= -Pi<T>) {
			radians += TwoPi<T>;
		}
		while (radians > Pi<T>) {
			radians -= TwoPi<T>;
		}

		return radians;
	}

	// a 3x3 matrix representation of a rotation
	template<typename T>
	class Rotation {
		public:

			Real3<T> xaxis;
			Real3<T> yaxis;
			Real3<T> zaxis;

			inline Rotation(const Real3<T> & xaxis, const Real3<T> & yaxis, const Real3<T> & zaxis):
				xaxis(xaxis), yaxis(yaxis), zaxis(zaxis) {}

			inline Rotation(): xaxis(1, 0, 0), yaxis(0, 1, 0), zaxis(0, 0, 1) {}

			inline void set_z(T radians) {
				T sin = std::sin(radians);
				T cos = std::cos(radians);
				xaxis = { cos, sin, 0 };
				yaxis = { -sin, cos, 0 };
				zaxis = { 0, 0, 1 };
			}

			inline void set_look(const Real3<T> & zaxis_unnorm, const Real3<T> & yaxis_unnorm) {

				// normalize z and negate it
				zaxis = -zaxis_unnorm;
				zaxis.normalize();

				// x = y x z
				xaxis = cross(yaxis_unnorm, zaxis);
				xaxis.normalize();

				// y = z x x
				yaxis = cross(zaxis, xaxis);
			}

			inline void invert() {
				// transpose the matrix
				std::swap(xaxis.y, yaxis.x);
				std::swap(xaxis.z, zaxis.x);
				std::swap(yaxis.z, zaxis.y);
			}
	};

	template<typename T>
	inline Real3<T> operator * (const Rotation<T> & r, const Real3<T> & v) {
		return {
			r.xaxis.dot(v),
			r.yaxis.dot(v),
			r.zaxis.dot(v)
		};
	}
}


#endif //CONFECALC_ROTATION_H
