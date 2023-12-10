
#ifndef CUDACONFECALC_CUDA_H
#define CUDACONFECALC_CUDA_H


#include <cooperative_groups.h>

namespace cg = cooperative_groups;


namespace osprey {

	// use cuda vector types for Real3

	template<typename T>
	struct Real3Map {
		typedef void type;
		const static size_t size;
	};

	template<>
	struct Real3Map<float32_t> {
		typedef float4 type;
		// yes, we really only need float3 here
		// but actually loading exactly 3 floats requires 2 load instructions
		// eg in PTX: ld.global.v2.f32 and ld.global.f32
		// use a float4 instead so we can use 1 load instruction
		// eg in PTX: ld.global.v4.f32
	};

	template<>
	struct Real3Map<float64_t> {
		typedef double3 type;
	};

	template<typename T>
	using Real3 = typename Real3Map<T>::type;

	// these are the sizes and alignments the compiler actually uses
	static_assert(sizeof(Real3<float32_t>) == 16);
	static_assert(alignof(Real3<float32_t>) == 16);

	static_assert(sizeof(Real3<float64_t>) == 24);
	static_assert(alignof(Real3<float64_t>) == 8);


	// add factory methods for Real3

	template<typename T>
	__device__
	inline Real3<T> real3(const T & x, const T & y, const T & z);

	template<>
	__device__
	inline float4 real3<float32_t>(const float32_t & x, const float32_t & y, const float32_t & z) {
		return make_float4(x, y, z, 0.0);
	}

	template<>
	__device__
	inline double3 real3<float64_t>(const float64_t & x, const float64_t & y, const float64_t & z) {
		return make_double3(x, y, z);
	}

	template<typename T>
	__device__
	inline Real3<T> real3(const int & x, const int & y, const int & z);

	template<>
	__device__
	inline float4 real3<float32_t>(const int & x, const int & y, const int & z) {
		return real3(
			static_cast<float32_t>(x),
			static_cast<float32_t>(y),
			static_cast<float32_t>(z)
		);
	}

	template<>
	__device__
	inline double3 real3<float64_t>(const int & x, const int & y, const int & z) {
		return real3(
			static_cast<float64_t>(x),
			static_cast<float64_t>(y),
			static_cast<float64_t>(z)
		);
	}


	// add math functions for vector types, since CUDA apparently doesn't have them in the stdlib ;_;

	template<typename T>
	__device__
	inline T dot(const Real3<T> & a, const Real3<T> & b) {
		return a.x*b.x + a.y*b.y + a.z*b.z;
	}

	template<typename T>
	__device__
	inline T distance_sq(const Real3<T> & a, const Real3<T> & b) {
		T dx = a.x - b.x;
		T dy = a.y - b.y;
		T dz = a.z - b.z;
		return dx*dx + dy*dy + dz*dz;
	}

	template<typename T>
	__device__
	inline bool operator == (const Real3<T> & a, const Real3<T> & b) {
		return a.x == b.x && a.y == b.y && a.z == b.z;
	}

	// nvcc can't find the templated operator for some reason, so explicitly instantiate it here
	__device__
	inline bool operator == (float4 & a, const float4 & b) {
		return operator ==<float32_t>(a, b);
	}
	__device__
	inline bool operator == (double3 & a, const double3 & b) {
		return operator ==<float64_t>(a, b);
	}

	template<typename T>
	__device__
	inline bool operator != (const Real3<T> & a, const Real3<T> & b) {
		return a.x != b.x || a.y != b.y || a.z != b.z;
	}

	// nvcc can't find the templated operator for some reason, so explicitly instantiate it here
	__device__
	inline bool operator != (float4 & a, const float4 & b) {
		return operator !=<float32_t>(a, b);
	}
	__device__
	inline bool operator != (double3 & a, const double3 & b) {
		return operator !=<float64_t>(a, b);
	}

	template<typename T>
	__device__
	inline void operator += (Real3<T> & self, const Real3<T> & other) {
		self.x += other.x;
		self.y += other.y;
		self.z += other.z;
	}

	// nvcc can't find the templated operator for some reason, so explicitly instantiate it here
	__device__
	inline void operator += (float4 & self, const float4 & other) {
		operator +=<float32_t>(self, other);
	}
	__device__
	inline void operator += (double3 & self, const double3 & other) {
		operator +=<float64_t>(self, other);
	}


	template<typename T>
	__device__
	inline void operator -= (Real3<T> & self, const Real3<T> & other) {
		self.x -= other.x;
		self.y -= other.y;
		self.z -= other.z;
	}

	// nvcc can't find the templated operator for some reason, so explicitly instantiate it here
	__device__
	inline void operator -= (float4 & self, const float4 & other) {
		operator -=<float32_t>(self, other);
	}
	__device__
	inline void operator -= (double3 & self, const double3 & other) {
		operator -=<float64_t>(self, other);
	}

	template<typename T>
	__device__
	inline Real3<T> operator - (const Real3<T> & v) {
		return {
			-v.x,
			-v.y,
			-v.z
		};
	}

	// nvcc can't find the templated operator for some reason, so explicitly instantiate it here
	__device__
	inline Real3<float32_t> operator - (const float4 & v) {
		return operator -<float32_t>(v);
	}
	__device__
	inline Real3<float64_t> operator - (const double3 & v) {
		return operator -<float64_t>(v);
	}

	template<typename T>
	__device__
	inline Real3<T> operator + (const Real3<T> & a, const Real3<T> & b) {
		return {
			a.x + b.x,
			a.y + b.y,
			a.z + b.z
		};
	}

	// nvcc can't find the templated operator for some reason, so explicitly instantiate it here
	__device__
	inline float4 operator + (const float4 & a, const float4 & b) {
		return operator +<float32_t>(a, b);
	}
	__device__
	inline double3 operator + (const double3 & a, const double3 & b) {
		return operator +<float64_t>(a, b);
	}

	template<typename T>
	__device__
	inline Real3<T> operator - (const Real3<T> & a, const Real3<T> & b) {
		return {
			a.x - b.x,
			a.y - b.y,
			a.z - b.z
		};
	}

	// nvcc can't find the templated operator for some reason, so explicitly instantiate it here
	__device__
	inline float4 operator - (const float4 & a, const float4 & b) {
		return operator -<float32_t>(a, b);
	}
	__device__
	inline double3 operator - (const double3 & a, const double3 & b) {
		return operator -<float64_t>(a, b);
	}

	template<typename T>
	__device__
	inline Real3<T> cross(const Real3<T> & a, const Real3<T> & b) {
		return {
			a.y*b.z - a.z*b.y,
			a.z*b.x - a.x*b.z,
			a.x*b.y - a.y*b.x
		};
	}

	template<typename T>
	__device__
	inline T len_sq(const Real3<T> & v) {
		return v.x*v.x + v.y*v.y + v.z*v.z;
	}

	template<typename T>
	__device__
	inline T len(const Real3<T> & v) {
		return std::sqrt(len_sq<T>(v));
	}

	template<typename T>
	__device__
	inline void normalize(Real3<T> & v) {
		T invlen = 1.0/len<T>(v);
		v.x *= invlen;
		v.y *= invlen;
		v.z *= invlen;
	}

	template<typename T>
	__device__
	inline bool isnan3(const Real3<T> & v) {
		return isnan<T>(v.x) || isnan<T>(v.y) || isnan<T>(v.z);
	}

	template<typename T>
	__device__
	inline void sincos_intr(T radians, T & sin, T & cos);

	template<>
	__device__
	inline void sincos_intr<float32_t>(float32_t radians, float32_t & sin, float32_t & cos) {
		__sincosf(radians, &sin, &cos);
	}

	template<>
	__device__
	inline void sincos_intr<float64_t>(float64_t radians, float64_t & sin, float64_t & cos) {
		// no intrinsics for double-precision ;_;
		sin = std::sin(radians);
		cos = std::cos(radians);
	}

	template<typename T>
	__device__
	inline T rsqrt_intr(T val);

	template<>
	__device__
	inline float32_t rsqrt_intr<float32_t>(float32_t val) {
		return __frsqrt_rn(val);
	}

	template<>
	__device__
	inline float64_t rsqrt_intr<float64_t>(float64_t val) {
		return __drcp_rn(__dsqrt_rn(val));
	}

	template<typename T>
	__device__
	inline T sqrt_intr(T val);

	template<>
	__device__
	inline float32_t sqrt_intr<float32_t>(float32_t val) {
		return __fsqrt_rn(val);
	}

	template<>
	__device__
	inline float64_t sqrt_intr<float64_t>(float64_t val) {
		return __dsqrt_rn(val);
	}

	template<typename T>
	__device__
	inline T rcp_intr(T val);

	template<>
	__device__
	inline float32_t rcp_intr<float32_t>(float32_t val) {
		return __frcp_rn(val);
	}

	template<>
	__device__
	inline float64_t rcp_intr<float64_t>(float64_t val) {
		return __drcp_rn(val);
	}


	template<typename T>
	__device__
	inline T exp_intr(T val);

	template<>
	__device__
	inline float32_t exp_intr<float32_t>(float32_t val) {
		return __expf(val);
	}

	template<>
	__device__
	inline float64_t exp_intr<float64_t>(float64_t val) {
		// no intrinsic for doubles here ;_;
		return std::exp(val);
	}
}


namespace cuda {

	__host__
	int get_arch();

	__host__
	int get_max_shared_size(int device);

	__host__
	int optimize_threads_void(const void * func, size_t shared_size_static, size_t shared_size_per_thread);

	// pick the greatest number of the threads that keeps occupancy above 0
	template<typename T>
	int optimize_threads(const T & func, size_t shared_size_static, size_t shared_size_per_thread) {
		return optimize_threads_void(reinterpret_cast<const void *>(&func), shared_size_static, shared_size_per_thread);
	}

	template<int A>
	__host__ __device__
	inline int64_t pad_to_alignment(int64_t size);

	template<>
	__host__ __device__
	inline int64_t pad_to_alignment<16>(int64_t size) {
		int8_t lsb = size & 0b1111;
		if (lsb == 0) {
			return size;
		} else {
			return size + 16 - lsb;
		}
	}

	template<int A>
	__host__ __device__
	inline bool is_aligned(int64_t offset);

	template<>
	__host__ __device__
	inline bool is_aligned<16>(int64_t offset) {
		return (offset & 0b1111) == 0;
	}

	template<int A>
	__host__ __device__
	inline bool is_aligned(const void * p) {
		return is_aligned<A>(reinterpret_cast<int64_t>(p));
	}
}

#define CUDACHECK(call)  \
	do { \
		cudaError_t result = call; \
		if (result != cudaSuccess) { \
			auto msg = cudaGetErrorString(result); \
			::std::cerr << "CUDA error @ " __FILE__ ":" S__LINE__ " " << msg << ::std::endl; \
			throw ::std::runtime_error(msg); \
		} \
	} while(0)

#define CUDACHECK_ERR()  CUDACHECK(cudaGetLastError())

#define CUDACHECK_SHAREDSIZE(device, requested_size) \
	do { \
		int max_size = cuda::get_max_shared_size(device); \
		if ((requested_size) > max_size) { \
			::std::cerr << "CUDA error @ " __FILE__ ":" S__LINE__ " requested too much shared memory " << (requested_size) << " but only have " << max_size << ::std::endl; \
			throw ::std::runtime_error("Requested too much shared memory"); \
		} \
	} while(0)

#define PRINTF0(fmt, ...) \
	if (threadIdx.x == 0) { \
		printf(fmt, __VA_ARGS__); \
	}


#endif //CUDACONFECALC_CUDA_H
