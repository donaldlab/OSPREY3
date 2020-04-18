
// make it easy to switch the precision of a "real" type
// between double/float at compile-time

#ifndef REAL
	#error "REAL must be either 'f64' or 'f32'"
#endif

// macro-land is complicated...
// need a few extra macros to print define values in pragma messages
#define VAL(x) __STR(x)
#define __STR(x) #x

// map our REAL values to numbers, so we can compare them
#define __REAL_f64 1
#define __REAL_f32 2

// do enough macro-foo to convert the input REAL value into a number
#define __REALID2(v) __REAL_ ## v
#define __REALID(v) __REALID2(v)
#define __REALID_IN __REALID(REAL)

// actually define types for REAL
// and an RL macro to make literals (RL = real literal)
#if __REALID_IN == __REAL_f64
	typedef double real_t;
	#define RL(val) val
#elif __REALID_IN == __REAL_f32
	typedef float real_t;
	#define RL(val) val ## f
#else
	#pragma message "REAL = " VAL(REAL)
	#error "Unknown REAL type"
#endif

