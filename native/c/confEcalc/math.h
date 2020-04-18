
#ifndef _MATH_
#define _MATH_


typedef struct real3 {
	real_t x;
	real_t y;
	real_t z;
	// 4 bytes pad iff real_t is float
} ALIGN_8 real3;

void print_real3(const real3 * p);


#endif // _MATH_

