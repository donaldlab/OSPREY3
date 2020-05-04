
#ifndef CONFECALC_GLOBAL_H
#define CONFECALC_GLOBAL_H


#include <cstdint>
#include <iostream>
#include <iomanip>
#include <cassert>
#include <iterator>
#include <cstring>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <cstdarg>


// define some types for fixed-size floating point types
static_assert(sizeof(float) == 4, "float type is not 32 bits");
typedef float float32_t;
static_assert(sizeof(double) == 8, "double type is not 64 bits");
typedef double float64_t;


// add some missing std functions

template<typename T>
__host__ __device__
inline T clamp(const T & val, const T & min, const T & max) {
	if (val < min) {
		return min;
	} else if (val > max) {
		return max;
	} else {
		return val;
	}
}

template<typename T>
__host__ __device__
inline void swap(T & a, T & b) {
	T temp = a;
	a = b;
	b = temp;
}

// calling std::isinf(T) or std::isnan(T) directly causes a compiler crash in nvcc V10.2.89:
// Error: Internal Compiler Error (codegen): "there was an error in verifying the lgenfe output!"
// wrapping it in a templated function seems to avoid the problem though
template<typename T>
__host__ __device__
inline bool isinf(const T & val) {
	return std::isinf(val);
}
template<typename T>
__host__ __device__
inline bool isnan(const T & val) {
	return std::isnan(val);
}


#include "config.h"


// HACKS! useful for debugging though
template<int N> constexpr char print_size_as_warning_char = N + 256;
#define WARN_SIZEOF(type) static char print_size_as_warning_var = print_size_as_warning_char<sizeof(type)>


#define ASSERT_MALLOCABLE(type, size) \
	static_assert(std::is_standard_layout<type>(), #type " should have standard layout"); \
	static_assert(sizeof(type) == size, #type " has unexpected size")

#define ASSERT_MALLOCABLE_REALS(type, size_float32, size_float64) \
	ASSERT_MALLOCABLE(type<float32_t>, size_float32); \
	ASSERT_MALLOCABLE(type<float64_t>, size_float64)

// java compatibility-checking macros
#define ASSERT_JAVA_COMPATIBLE(type, size) ASSERT_MALLOCABLE(type, size)
#define ASSERT_JAVA_COMPATIBLE_REALS(type, size_float32, size_float64) ASSERT_MALLOCABLE_REALS(type, size_float32, size_float64)

// other type-checking macros
#define ASSERT_COPYABLE(type) static_assert(std::is_trivially_copyable<type>(), #type " should be trivially copyable")

#define ASSERT_COPYABLE_REALS(type) \
	ASSERT_COPYABLE(type<float32_t>); \
	ASSERT_COPYABLE(type<float64_t>)


#endif //CONFECALC_GLOBAL_H
