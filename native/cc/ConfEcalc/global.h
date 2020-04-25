
#ifndef CONFECALC_GLOBAL_H
#define CONFECALC_GLOBAL_H


#include <cstdint>
#include <iostream>
#include <iomanip>
#include <cassert>
#include <iterator>
#include <cstring>
#include <cmath>


// define some types for fixed-size floating point types
static_assert(sizeof(float) == 4, "float type is not 32 bits");
typedef float float32_t;
static_assert(sizeof(double) == 8, "double type is not 64 bits");
typedef double float64_t;

#include "config.h"


// HACKS! useful for debugging though
template<int N> constexpr char print_size_as_warning_char = N + 256;
#define WARN_SIZEOF(type) static char print_size_as_warning_var = print_size_as_warning_char<sizeof(type)>


// java compatibility-checking macros
#define ASSERT_JAVA_COMPATIBLE(type, size) \
	static_assert(std::is_standard_layout<type>(), #type " should have standard layout"); \
	static_assert(sizeof(type) == size, #type " has unexpected size")

#define ASSERT_JAVA_COMPATIBLE_REALS(type, size_float32, size_float64) \
	ASSERT_JAVA_COMPATIBLE(type<float32_t>, size_float32); \
	ASSERT_JAVA_COMPATIBLE(type<float64_t>, size_float64)

// other type-checking macros
#define ASSERT_COPYABLE(type) static_assert(std::is_trivially_copyable<type>(), #type " should be trivially copyable")

#define ASSERT_COPYABLE_REALS(type) \
	ASSERT_COPYABLE(type<float32_t>); \
	ASSERT_COPYABLE(type<float64_t>)


#endif //CONFECALC_GLOBAL_H
