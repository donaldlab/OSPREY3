
#ifdef cl_khr_fp64
    #pragma OPENCL EXTENSION cl_khr_fp64 : enable
#elif defined(cl_amd_fp64)
    #pragma OPENCL EXTENSION cl_amd_fp64 : enable
#else
    #error "Double precision floating point not supported by OpenCL implementation."
#endif


// NOTE: for some reason OpenCL doesn't support the "const" keyword here
constant double Tolerance = 1e-6;


typedef struct __attribute__((aligned(8))) {

	int numRotatedIndices; // @ 0
	int surfStep; // @ 4
	
	double xdmin; // @ 8, 18 doubles
	double xdmax; // @ 16
	double xd; // @ 24
	double fxd; // @ 32
	double fxdmin; // @ 40
	double fxdmax; // @ 48
	double step; // @ 56
	double xdm; // @ 64
	double xdp; // @ 72
	double fxdm; // @ 80
	double fxdp; // @ 88
	double xdstar; // @ 96
	double fxdstar; // @ 104
	double stepScale; // @ 112
	double xdsurfHere; // @ 120
	double fxdsurfHere; // @ 128
	double dihedralDegrees; // @ 136, in degrees
	double internalSolvEnergy; // @ 144
	
	bool isSurfing; // @ 152
	
	// 7 pad
	
} Args; // size = 160


typedef struct __attribute__((aligned(8))) {
	int numPairs;
	int num14Pairs;
	double coulombFactor;
	double scaledCoulombFactor;
	double solvCutoff2;
	bool useDistDepDielec;
	bool useHEs;
	bool useHVdw;
	bool doEnergy;
} ForcefieldArgs;


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

double3 rotateVecZ(double3 v, double sinTheta, double cosTheta) {
	return (double3)(
		v.x*cosTheta - v.y*sinTheta,
		v.x*sinTheta + v.y*cosTheta,
		v.z
	);
}

kernel void pose(global double *coords, global const int *dihedralIndices, global const int *rotatedIndices, global ForcefieldArgs *ffargs, global Args *args) {

	// short circuit if we don't need the energy
	if (!ffargs->doEnergy) {
		return;
	}

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
	newSin = sincos(radians(args->dihedralDegrees), &newCos);
	double deltaSin = newSin*currentCos - newCos*currentSin;
	double deltaCos = newCos*currentCos + newSin*currentSin;
	
	// apply transformation to every downstream atom (in parallel)
	int index = rotatedIndices[get_global_id(0)];
		
	double3 p = readCoord(coords, index);
	p -= b;
	p = rotateVecInverse(p, rx, ry, rz);
	p = rotateVecZ(p, deltaSin, deltaCos);
	p = rotateVec(p, rx, ry, rz);
	p += b;
	writeCoord(coords, index, p);
}

void sumEnergy(global const double* energies, const int numEnergies, global Args *args, local double *scratch) {

	// read the energy buffer in parallel
	int id = get_global_id(0);
	int numThreads = get_global_size(0);
	
	// partition the reads
	int width = (numEnergies + numThreads - 1)/numThreads;
	int start = id*width;
	int stop = min(start + width, numEnergies);
	
	// read energies for our thread
	double energy = 0;
	for (int i=start; i<stop; i++) {
		energy += energies[i];
	}
	scratch[id + 1] = energy;
	
	barrier(CLK_LOCAL_MEM_FENCE);
	
	// do the final sum: read from scratch[1:n+1] and write to scratch[0]
	energy = args->internalSolvEnergy;
	for (int i=0; i<numThreads; i++) {
		energy += scratch[i + 1];
	}
	scratch[0] = energy;
}

double getTolerance(double f) {
	
	// use full tolerance, unless f is very small
	// then scale by the magnitude of f
	return Tolerance * fmax(1.0, fabs(f));
}

kernel void search0(global const double *energies, const int numEnergies, global ForcefieldArgs *ffargs, global Args *args, local double *scratch) {

	// sum the energies
	if (ffargs->doEnergy) {
		sumEnergy(energies, numEnergies, args, scratch);
	}
	if (get_global_id(0) != 0) {
		return;
	}
	
	// read energy result
	args->fxd = scratch[0];
	
	// get the min edge
	ffargs->doEnergy = true;
	args->dihedralDegrees = args->xdmin;
}

kernel void search1(global const double *energies, const int numEnergies, global ForcefieldArgs *ffargs, global Args *args, local double *scratch) {

	// sum the energies
	sumEnergy(energies, numEnergies, args, scratch);
	if (get_global_id(0) != 0) {
		return;
	}
	
	// read energy result
	args->fxdmin = scratch[0];
	
	// get the max edge
	ffargs->doEnergy = true;
	args->dihedralDegrees = args->xdmax;
}

kernel void search2(global const double *energies, const int numEnergies, global ForcefieldArgs *ffargs, global Args *args, local double *scratch) {

	// sum the energies
	if (ffargs->doEnergy) {
		sumEnergy(energies, numEnergies, args, scratch);
	}
	if (get_global_id(0) != 0) {
		return;
	}
	
	// read energy result
	args->fxdmax = scratch[0];
	
	// get the minus neighbor
	args->xdm = args->xd - args->step;
	
	// init next energy calc
	ffargs->doEnergy = args->xdm >= args->xdmin;
	args->dihedralDegrees = args->xdm;
}

kernel void search3(global const double *energies, const int numEnergies, global ForcefieldArgs *ffargs, global Args *args, local double *scratch) {

	// sum the energies
	if (ffargs->doEnergy) {
		sumEnergy(energies, numEnergies, args, scratch);
	}
	if (get_global_id(0) != 0) {
		return;
	}
	
	// read energy result
	args->fxdm = INFINITY;
	if (ffargs->doEnergy) {
		args->fxdm = scratch[0];
	}
	
	// get the plus neighbor
	args->xdp = args->xd + args->step;
	
	// init next energy calc
	ffargs->doEnergy = args->xdp <= args->xdmax;
	args->dihedralDegrees = args->xdp;
}

kernel void search4(global const double *energies, const int numEnergies, global ForcefieldArgs *ffargs, global Args *args, local double *scratch) {

	// sum the energies
	if (ffargs->doEnergy) {
		sumEnergy(energies, numEnergies, args, scratch);
	}
	if (get_global_id(0) != 0) {
		return;
	}
	
	// read energy result
	args->fxdp = INFINITY;
	if (ffargs->doEnergy) {
		args->fxdp = scratch[0];
	}
	
	// fit a quadratic to the objective function, locally:
	// q(x) = fx + a*(x - xd)^2 + b*(x - xd)
	// a*step^2 + b*step = fxp - fx
	// a*step^2 - b*step = fxm - fx
	
	// solve for a to determine the shape
	double a = (args->fxdp + args->fxdm - 2*args->fxd)/(2*args->step*args->step);
	args->xdstar = 0;
	if (a <= 0 || a == NAN || a == INFINITY) {
		
		// negative a means quadratic is concave down, I think
		// infinite or nan a means we're hitting a constraint or impossible conformation
		// so just minimize over the endpoints of the interval
		if (args->fxdm < args->fxdp) {
			args->xdstar = args->xdm;
		} else {
			args->xdstar = args->xdp;
		}
		
	} else {
		
		// positive a means quadratic is concave up, I think
		// solve for the b param
		double b = (args->fxdp - args->fxd)/args->step - a*args->step;
		
		// then minimize the quadratic to get the minimum x:
		// 2*a*(x - xd) + b = 0
		args->xdstar = args->xd - b/2/a;
	}

	// clamp xdstar to the range
	if (args->xdstar < args->xdmin) {
		args->xdstar = args->xdmin;
	}
	if (args->xdstar > args->xdmax) {
		args->xdstar = args->xdmax;
	}
	
	// init next energy calc
	ffargs->doEnergy = true;
	args->dihedralDegrees = args->xdstar;
}

kernel void search5(global const double *energies, const int numEnergies, global ForcefieldArgs *ffargs, global Args *args, local double *scratch) {

	// sum the energies
	if (ffargs->doEnergy) {
		sumEnergy(energies, numEnergies, args, scratch);
	}
	if (get_global_id(0) != 0) {
		return;
	}
	
	// read energy result
	args->fxdstar = scratch[0];
	
	// if we went downhill, use growing steps
	if (args->fxdstar < args->fxd + getTolerance(args->fxd)) {
		args->stepScale = 2;
	} else {
		args->stepScale = 0.5;
	}
	
	// surf along f locally to try to find better minimum
	args->isSurfing = true;
	args->surfStep = 0;
	args->xdsurfHere = args->xdstar;
	args->fxdsurfHere = args->fxdstar;
	ffargs->doEnergy = false;
}

kernel void surf(global const double *energies, const int numEnergies, global ForcefieldArgs *ffargs, global Args *args, local double *scratch) {

	// sum the energies
	if (ffargs->doEnergy) {
		sumEnergy(energies, numEnergies, args, scratch);
	}
	if (get_global_id(0) != 0) {
		return;
	}
	
	ffargs->doEnergy = false;
		
	if (args->isSurfing && args->surfStep > 0) {
			
		double xdsurfNext = args->dihedralDegrees;
		double fxdsurfNext = scratch[0];
		
		// did we improve the min enough to keep surfing?
		if (fxdsurfNext < args->fxdsurfHere - getTolerance(args->fxdsurfHere)) {
		
			// yeah, keep going
			args->xdsurfHere = xdsurfNext;
			args->fxdsurfHere = fxdsurfNext;
			
		} else {
			
			// nope, stop surfing
			args->isSurfing = false;
		}
	}
		
	if (args->isSurfing) {
	
		// surf a step
		double xdsurfNext = args->xd + args->stepScale*(args->xdsurfHere - args->xd);
		
		// did we step off the min?
		if (xdsurfNext < args->xdmin) {
			
			// if the min is better, go there instead
			if (args->fxdmin < args->fxdstar) {
				args->xdsurfHere = args->xdmin;
			}
			
			args->isSurfing = false;
		
		// did we step off the max?
		} else if (xdsurfNext > args->xdmax) {
			
			// if the max is better, go there instead
			if (args->fxdmax < args->fxdstar) {
				args->xdsurfHere = args->xdmax;
			}
			
			args->isSurfing = false;
			
		} else {
		
			// init energy calc for new surf spot
			args->surfStep++;
			ffargs->doEnergy = true;
			args->dihedralDegrees = xdsurfNext;
		}
	}
}

kernel void search6(global const double *energies, const int numEnergies, global ForcefieldArgs *ffargs, global Args *args, local double *scratch) {

	// did the quadratic step help at all?
	if (args->fxdstar < args->fxd) {
		
		// yeah, keep it!
		
	} else {
		
		// nope, the original spot was lower
		args->xdstar = args->xd;
		args->fxdstar = args->fxd;
	}
	
	// did surfing help at all?
	if (args->fxdsurfHere < args->fxdstar) {
		
		// yeah, use the surf spot
		args->xdstar = args->xdsurfHere;
		args->fxdstar = args->fxdsurfHere;
	}
	
	// try to jump over walls arbitrarily
	// look in a 1-degree step for a better minimum
	
	// NOTE: skipping this can make minimization a bit faster,
	// but skipping this causes a noticeable rise in final energies too
	// it's best to keep doing it I think
	
	// init next energy calc
	double xdm = args->xdstar - 1;
	ffargs->doEnergy = xdm >= args->xdmin;
	args->dihedralDegrees = xdm;
}

kernel void search7(global const double *energies, const int numEnergies, global ForcefieldArgs *ffargs, global Args *args, local double *scratch) {

	// sum the energies
	if (ffargs->doEnergy) {
		sumEnergy(energies, numEnergies, args, scratch);
	}
	if (get_global_id(0) != 0) {
		return;
	}
	
	double xdstarOld = args->xdstar;

	// read energy result
	if (ffargs->doEnergy) {
		double energy = scratch[0];
		if (energy < args->fxdstar - getTolerance(args->fxdstar)) {
			args->xdstar = args->dihedralDegrees;
			args->fxdstar = energy;
		}
	}
	
	// init next energy calc
	double xdp = xdstarOld + 1;
	ffargs->doEnergy = xdp <= args->xdmax;
	args->dihedralDegrees = xdp;
}

kernel void search8(global const double *energies, const int numEnergies, global ForcefieldArgs *ffargs, global Args *args, local double *scratch) {

	// sum the energies
	if (ffargs->doEnergy) {
		sumEnergy(energies, numEnergies, args, scratch);
	}
	if (get_global_id(0) != 0) {
		return;
	}
	
	// read energy result
	if (ffargs->doEnergy) {
		double energy = scratch[0];
		if (energy < args->fxdstar - getTolerance(args->fxdstar)) {
			args->xdstar = args->dihedralDegrees;
			args->fxdstar = energy;
		}
	}
	
	// update the step to the one we actually took
	args->step = args->xdstar - args->xd;
	
	// init one last pose
	ffargs->doEnergy = true;
	args->dihedralDegrees = args->xdstar;
}
