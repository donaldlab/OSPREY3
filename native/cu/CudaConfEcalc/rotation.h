
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
	__device__
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

			__device__
			inline Rotation(const Real3<T> & xaxis, const Real3<T> & yaxis, const Real3<T> & zaxis):
				xaxis(xaxis), yaxis(yaxis), zaxis(zaxis) {}
			__device__
			inline Rotation(const Rotation<T> & other):
				xaxis(other.xaxis), yaxis(other.yaxis), zaxis(other.zaxis) {}
			~Rotation() = default;

			__device__
			inline Rotation() {
				set_identity();
			}

			__device__
			inline void set_identity() {
				xaxis = real3<T>(1, 0, 0);
				yaxis = real3<T>(0, 1, 0);
				zaxis = real3<T>(0, 0, 1);
			}

			__device__
			inline void set_z(T radians) {
				T sin = std::sin(radians);
				T cos = std::cos(radians);
				xaxis = real3<T>(cos, -sin, 0.0);
				yaxis = real3<T>(sin, cos, 0.0);
				zaxis = real3<T>(0, 0, 1);
			}

			__device__
			inline void set_look(const Real3<T> & zaxis_unnorm, const Real3<T> & yaxis_unnorm) {

				// normalize z and negate it
				zaxis = -zaxis_unnorm;
				normalize<T>(zaxis);

				// x = y x z
				xaxis = cross<T>(yaxis_unnorm, zaxis);
				normalize<T>(xaxis);

				// y = z x x
				yaxis = cross<T>(zaxis, xaxis);
			}

			__device__
			inline void invert() {
				// transpose the matrix
				swap(xaxis.y, yaxis.x);
				swap(xaxis.z, zaxis.x);
				swap(yaxis.z, zaxis.y);
			}
	};

	template<typename T>
	__device__
	inline Real3<T> operator * (const Rotation<T> & r, const Real3<T> & v) {
		return real3<T>(
			dot<T>(r.xaxis, v),
			dot<T>(r.yaxis, v),
			dot<T>(r.zaxis, v)
		);
	}
}


#endif //CONFECALC_ROTATION_H
