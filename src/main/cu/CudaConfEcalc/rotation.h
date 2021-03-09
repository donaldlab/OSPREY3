
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

	// a 3x3 column-major matrix representation of a rotation
	template<typename T>
	class Rotation {
		public:

			// these are row vectors
			Real3<T> xaxis;
			Real3<T> yaxis;
			Real3<T> zaxis;

			inline Rotation() = default;
			inline Rotation(const Rotation<T> & other) = default;
			~Rotation() = default;

			__device__
			inline Rotation(const Real3<T> & xaxis, const Real3<T> & yaxis, const Real3<T> & zaxis):
				xaxis(xaxis), yaxis(yaxis), zaxis(zaxis) {}

			__device__
			inline void set_identity() {
				xaxis = real3<T>(1, 0, 0);
				yaxis = real3<T>(0, 1, 0);
				zaxis = real3<T>(0, 0, 1);
			}

			__device__
			inline void set_xyz(T radians_x, T radians_y, T radians_z) {

				T sin_y, cos_y;
				sincos_intr(radians_y, sin_y, cos_y);
				T sin_z, cos_z;
				sincos_intr(radians_z, sin_z, cos_z);

				xaxis = real3(
					cos_y*cos_z,
					-cos_y*sin_z,
					sin_y
				);

				T sin_x, cos_x;
				sincos_intr(radians_x, sin_x, cos_x);

				yaxis = real3(
					sin_x*sin_y*cos_z + cos_x*sin_z,
					-sin_x*sin_y*sin_z + cos_x*cos_z,
					-sin_x*cos_y
				);

				zaxis = real3(
					-cos_x*sin_y*cos_z + sin_x*sin_z,
					cos_x*sin_y*sin_z + sin_x*cos_z,
					cos_x*cos_y
				);
			}

			__device__
			inline void set_look(const Real3<T> & zaxis_unnorm, const Real3<T> & xaxis_unnorm) {

				// normalize Z and negate it
				zaxis = -zaxis_unnorm;
				normalize<T>(zaxis);

				// Y = ZxX
				yaxis = cross<T>(zaxis, xaxis_unnorm);
				normalize<T>(yaxis);

				// X = YxZ
				xaxis = cross<T>(yaxis, zaxis);
			}

			__device__
			inline void invert() {
				// transpose the matrix
				swap(xaxis.y, yaxis.x);
				swap(xaxis.z, zaxis.x);
				swap(yaxis.z, zaxis.y);
			}

			__device__
			inline void mul(Real3<T> & v) {
				T x = dot<T>(xaxis, v);
				T y = dot<T>(yaxis, v);
				v.z = dot<T>(zaxis, v);
				v.x = x;
				v.y = y;
			}

			__device__
			inline void mul_inv(Real3<T> & v) {
				T x = xaxis.x*v.x + yaxis.x*v.y + zaxis.x*v.z;
				T y = xaxis.y*v.x + yaxis.y*v.y + zaxis.y*v.z;
				v.z = xaxis.z*v.x + yaxis.z*v.y + zaxis.z*v.z;
				v.x = x;
				v.y = y;
			}
	};

	// a 2x2 column-major matrix representation of a rotation
	template<typename T>
	class alignas(16) RotationZ {
		public:
			T sin;
			T cos;

			inline RotationZ() = default;
			inline RotationZ(const RotationZ<T> & other) = default;
			~RotationZ() = default;

			__device__
			inline void set(T radians) {
				sincos_intr(radians, sin, cos);
			}

			__device__
			inline void set(T srcx, T srcy, T radians) {

				// just in case ... because it actually happened
				assert (!isnan<T>(srcx));
				assert (!isnan<T>(srcy));
				assert (!isnan<T>(radians));
				assert (!(srcx == 0 && srcy == 0));
				assert (srcx*srcx + srcy*srcy > 0);

				T oolen = rsqrt_intr(srcx*srcx + srcy*srcy);
				T sinsrc = srcy*oolen;
				T cossrc = srcx*oolen;
				T sindst;
				T cosdst;
				sincos_intr(radians, sindst, cosdst);
				sin = sindst*cossrc;
				sin -= cosdst*sinsrc;
				cos = cosdst*cossrc;
				cos += sindst*sinsrc;
			}

			__device__
			inline void mul(Real3<T> & v) {
				T x = cos*v.x - sin*v.y;
				v.y = sin*v.x + cos*v.y;
				v.x = x;
			}
	};
}


#endif //CONFECALC_ROTATION_H
