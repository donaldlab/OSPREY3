
#ifdef cl_khr_fp64
    #pragma OPENCL EXTENSION cl_khr_fp64 : enable
#elif defined(cl_amd_fp64)
    #pragma OPENCL EXTENSION cl_amd_fp64 : enable
#else
    #error "Double precision floating point not supported by OpenCL implementation."
#endif

typedef struct __attribute__((aligned(8))) {
	int numRotatedIndices; // 4
	// pad 4
	double newDihedral; // in radians // 8
} Args; // size = 16

double3 readCoord(global double *coords, int i) {
	int i3 = i*3;
	return (double3)(coords[i3], coords[i3 + 1], coords[i3 + 2]);
}

void writeCoord(global double *coords, int i, double3 val) {
	int i3 = i*3;
	coords[i3] = val.x;
	coords[i3 + 1] = val.y;
	coords[i3 + 2] = val.z;
}

double3 rotateVec(double3 v, double3 x, double3 y, double3 z) {
	return (double3)(
		v.x*x.x + v.y*y.x + v.z*z.x,
		v.x*x.y + v.y*y.y + v.z*z.y,
		v.x*x.z + v.y*y.z + v.z*z.z
	);
}
    
double3 rotateVecInverse(double3 v, double3 x, double3 y, double3 z) {
	return (double3)(
		dot(v, x),
		dot(v, y),
		dot(v, z)
	);
}

kernel void pose(global double *coords, global const int *dihedralIndices, global const int *rotatedIndices, global Args *args) {

	// get the four atom positions: a, b, c, d
	double3 a = readCoord(coords, dihedralIndices[0]);
	double3 b = readCoord(coords, dihedralIndices[1]);
	double3 c = readCoord(coords, dihedralIndices[2]);
	double3 d = readCoord(coords, dihedralIndices[3]);
	
	// translate so everything is centered on b
	a -= b;
	c -= b;
	d -= b;
	
	// build a right orthnormal matrix [rx,ry,rz] where z is bc and ba points along x
	double3 rz = normalize(c);
	double3 rx = normalize(a - c*dot(a, c)/dot(c, c));
	double3 ry = cross(rz, rx);
	
	// use r^{-1} to rotate d into our axis-aligned space
	d = rotateVecInverse(d, rx, ry, rz);
	
	// look at the x,y coords of d to get the dihedral angle
	double2 cossin = normalize(d.xy);
	double currentSin = cossin.y;
	double currentCos = cossin.x;
	
	// get the delta dihedral
	double newSin, newCos;
	newSin = sincos(args->newDihedral, &newCos);
	double deltaSin = newSin*currentCos - newCos*currentSin;
	double deltaCos = newCos*currentCos + newSin*currentSin;
	
	// build the delta rotation matrix about z by angle deltaDihedral
	double3 dx = (double3)( deltaCos, deltaSin, 0 );
	double3 dy = (double3)( -deltaSin, deltaCos, 0 );
	double3 dz = (double3)( 0, 0, 1 );
	
	// apply transformation to every downstream atom (in parallel)
	int index = rotatedIndices[get_global_id(0)];
		
	double3 p = readCoord(coords, index);
	p -= b;
	p = rotateVecInverse(p, rx, ry, rz);
	p = rotateVec(p, dx, dy, dz);
	p = rotateVec(p, rx, ry, rz);
	p += b;
	writeCoord(coords, index, p);
}
