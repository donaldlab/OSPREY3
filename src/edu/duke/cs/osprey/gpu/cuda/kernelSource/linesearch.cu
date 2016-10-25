
/* compile with
	nvcc -cubin -dlink
	-gencode arch=compute_35,code=sm_35
	-gencode arch=compute_50,code=sm_50
	-gencode arch=compute_52,code=sm_52
	"kernelSource/linesearch.cu" -o "kernelBinaries/linesearch.bin"  -L /usr/lib/x86_64-linux-gnu
	
	see: http://docs.nvidia.com/cuda/maxwell-compatibility-guide/index.html#building-maxwell-compatible-apps-using-cuda-6-0
*/

// CUDA quick reference:
// http://www.icl.utk.edu/~mgates3/docs/cuda.html
// http://docs.nvidia.com/cuda/cuda-math-api/group__CUDA__MATH__DOUBLE.html#group__CUDA__MATH__DOUBLE


#include <vector_types.h>
#include <stdio.h>
#include <math.h>

const double Tolerance = 1e-6;
const double OneDegree = 0.017453293; // in radians


//===========================
//   forcefield stuff
//===========================


typedef struct __align__(8) {
	int numPairs; // @ 0
	int num14Pairs; // @ 4
	double coulombFactor; // @ 8
	double scaledCoulombFactor; // @ 16
	double solvCutoff2; // @ 24
	bool useDistDepDielec; // @ 32
	bool useHEs; // @ 33
	bool useHVdw; // @ 34
	bool doEnergy; // @ 35 // TODO: get rid of this
} ForcefieldArgs;
// sizeof = 36


__device__ int getAtomIndex(int flags) {
	return abs(flags) - 1;
}

__device__ bool isHydrogen(int flags) {
	return flags > 0;
}


__global__ void energyKernel(
	const double *coords,
	const int *atomFlags,
	const double *precomputed,
	const int *subsetTable,
	const ForcefieldArgs *args,
	double *blockEnergies
) {

	extern __shared__ double scratch[];

	// start with zero energy
	double energy = 0;
	
	int globalId = blockIdx.x*blockDim.x + threadIdx.x;
	
	// which atom pair are we calculating?
	if (globalId < args->numPairs) {
		int i = subsetTable[globalId];
		
		// read atom flags and calculate all the things that use the atom flags in this scope
		bool bothHeavy;
		double r2 = 0;
		{
			int atom1Flags, atom2Flags;
			{
				int i2 = i*2;
				atom1Flags = atomFlags[i2];
				atom2Flags = atomFlags[i2 + 1];
			}
			
			bothHeavy = !isHydrogen(atom1Flags) && !isHydrogen(atom2Flags);
			
			// calculate the squared radius
			int atom1Index3 = getAtomIndex(atom1Flags)*3;
			int atom2Index3 = getAtomIndex(atom2Flags)*3;
			double d;
			d = coords[atom1Index3] - coords[atom2Index3];
			r2 += d*d;
			d = coords[atom1Index3 + 1] - coords[atom2Index3 + 1];
			r2 += d*d;
			d = coords[atom1Index3 + 2] - coords[atom2Index3 + 2];
			r2 += d*d;
		}
		
		int i9 = i*9;
		
		// calculate electrostatics
		if (bothHeavy || args->useHEs) {
		
			double esEnergy = 1;
		
			{
				bool is14Pair = globalId < args->num14Pairs;
				esEnergy *= is14Pair ? args->scaledCoulombFactor : args->coulombFactor;
			}
			
			{
				double charge = precomputed[i9 + 2];
				esEnergy *= charge;
			}
			
			{
				esEnergy /= args->useDistDepDielec ? r2 : sqrt(r2);
			}
			
			energy += esEnergy;
		}
		
		// calculate vdw
		if (bothHeavy || args->useHVdw) {
			
			double Aij, Bij;
			{
				Aij = precomputed[i9];
				Bij = precomputed[i9 + 1];
			}
			
			// compute vdw
			double r6 = r2*r2*r2;
			double r12 = r6*r6;
			energy += Aij/r12 - Bij/r6;
		}
		
		// calculate solvation
		if (bothHeavy && r2 < args->solvCutoff2) {
				
			double r = sqrt(r2);
			{
				double lambda1 = precomputed[i9 + 3];
				double radius1 = precomputed[i9 + 4];
				double alpha1 = precomputed[i9 + 5];
				double Xij = (r - radius1)/lambda1;
				energy -= alpha1*exp(-Xij*Xij)/r2;
			}
			{
				double lambda2 = precomputed[i9 + 6];
				double radius2 = precomputed[i9 + 7];
				double alpha2 = precomputed[i9 + 8];
				double Xji = (r - radius2)/lambda2;
				energy -= alpha2*exp(-Xji*Xji)/r2;
			}
		}
	}
	
	// compute the energy sum in SIMD-style
	// see url for a tutorial on GPU reductions:
	// http://developer.amd.com/resources/articles-whitepapers/opencl-optimization-case-study-simple-reductions/

	int localId = threadIdx.x;
	scratch[localId] = energy;
	
	__syncthreads();
	
	for (int offset = 1; offset < blockDim.x; offset <<= 1) {
	
		// sum this level of the reduction tree
		int mask = (offset << 1) - 1;
		if ((localId & mask) == 0) {
			scratch[localId] += scratch[localId + offset];
		}
		
		__syncthreads();
	}
	
	// finally, if we're the 0 thread, write the summed energy for this work group
	if (localId == 0) {
		blockEnergies[blockIdx.x] = scratch[0];
	}
}

__device__ double energy(int numBlocks, int blockThreads, const double *coords, const int *atomFlags, const double *precomputed, const int *subsetTable, const ForcefieldArgs *ffargs, double *blockEnergies) {

	// launch the energy kernel
	int sharedBytes = blockThreads*sizeof(double);
	energyKernel<<<numBlocks, blockThreads, sharedBytes>>>(coords, atomFlags, precomputed, subsetTable, ffargs, blockEnergies);
	cudaDeviceSynchronize();
	
	// sum block energies
	double energy = 0;
	for (int i=0; i<numBlocks; i++) {
		energy += blockEnergies[i];
	}
	return energy;
}


//===========================
//   dihedral posing stuff
//===========================

// jesus, CUDA... no libraries for vector math? what gives??

__device__ void set(double2 &v, double x, double y) {
	v.x = x;
	v.y = y;
}

__device__ void set(double3 &v, double x, double y, double z) {
	v.x = x;
	v.y = y;
	v.z = z;
}

__device__ void sub(double3 &a, double3 &b) {
	a.x -= b.x;
	a.y -= b.y;
	a.z -= b.z;
}

__device__ void add(double3 &a, double3 &b) {
	a.x += b.x;
	a.y += b.y;
	a.z += b.z;
}

__device__ double dot(double2 &a, double2 &b) {
	return a.x*b.x + a.y*b.y;
}

__device__ double dot(double3 &a, double3 &b) {
	return a.x*b.x + a.y*b.y + a.z*b.z;
}

__device__ void cross(double3 &out, double3 &a, double3 &b) {
	out.x = a.y*b.z - a.z*b.y;
	out.y = a.z*b.x - a.x*b.z;
	out.z = a.x*b.y - a.y*b.x;
}

__device__ double lengthSq(double2 &v) {
	return dot(v, v);
}

__device__ double lengthSq(double3 &v) {
	return dot(v, v);
}

__device__ double length(double2 &v) {
	return sqrt(lengthSq(v));
}

__device__ double length(double3 &v) {
	return sqrt(lengthSq(v));
}

__device__ void negate(double3 &v) {
	v.x = -v.x;
	v.y = -v.y;
	v.z = -v.z;
}

__device__ void mult(double3 &v, double c) {
	v.x *= c;
	v.y *= c;
	v.z *= c;
}

__device__ void div(double3 &v, double c) {
	v.x /= c;
	v.y /= c;
	v.z /= c;
}

__device__ void normalize(double2 &v) {
	double l = length(v);
	v.x /= l;
	v.y /= l;
}

__device__ void normalize(double3 &v) {
	double l = length(v);
	v.x /= l;
	v.y /= l;
	v.z /= l;
}

__device__ void rotateVec(double3 &v, double3 &x, double3 &y, double3 &z) {
	set(v,
		v.x*x.x + v.y*y.x + v.z*z.x,
		v.x*x.y + v.y*y.y + v.z*z.y,
		v.x*x.z + v.y*y.z + v.z*z.z
	);
}

__device__ void rotateVecInverse(double3 &v, double3 &x, double3 &y, double3 &z) {
	set(v,
		dot(v, x),
		dot(v, y),
		dot(v, z)
	);
}

__device__ void rotateVecZ(double3 &v, double &sinTheta, double &cosTheta) {
	double vx = v.x*cosTheta - v.y*sinTheta;
	double vy = v.x*sinTheta + v.y*cosTheta;
	v.x = vx;
	v.y = vy;
}

__device__ double3 readCoord(double *coords, int i) {
	int i3 = i*3;
	return make_double3(coords[i3], coords[i3 + 1], coords[i3 + 2]);
}

__device__ void writeCoord(double *coords, int i, double3 &val) {
	int i3 = i*3;
	coords[i3] = val.x;
	coords[i3 + 1] = val.y;
	coords[i3 + 2] = val.z;
}

__device__ void pose(double *coords, const int *dihedralIndices, const int numRotatedIndices, const int *rotatedIndices, const double dihedralRadians) {

	// get the four atom positions: a, b, c, d
	double3 a = readCoord(coords, dihedralIndices[0]);
	double3 b = readCoord(coords, dihedralIndices[1]);
	double3 c = readCoord(coords, dihedralIndices[2]);
	double3 d = readCoord(coords, dihedralIndices[3]);
	
	// translate so everything is centered on b
	sub(a, b);
	sub(c, b);
	sub(d, b);
	
	// build a right orthonormal matrix [rx,ry,rz] where z is bc and ba points along x
	double3 rz = c;
	normalize(rz);
	
	double3 rx = c;
	mult(rx, dot(a, c)/dot(c, c));
	negate(rx);
	add(rx, a);
	normalize(rx);
	
	double3 ry;
	cross(ry, rz, rx);
	
	// use r^{-1} to rotate d into our axis-aligned space
	rotateVecInverse(d, rx, ry, rz);
	
	// look at the x,y coords of d to get the dihedral angle
	double2 cossin = make_double2(d.x, d.y);
	normalize(cossin);
	double currentSin = cossin.y;
	double currentCos = cossin.x;
	
	// get the delta dihedral
	double newSin, newCos;
	sincos(dihedralRadians, &newSin, &newCos);
	double deltaSin = newSin*currentCos - newCos*currentSin;
	double deltaCos = newCos*currentCos + newSin*currentSin;
	
	// there are only 10s of atoms rotate
	// it's probably not worth doing a kernel launch here
	for (int i=0; i<numRotatedIndices; i++) {
		int index = rotatedIndices[i];
		double3 p = readCoord(coords, index);
		sub(p, b);
		rotateVecInverse(p, rx, ry, rz);
		rotateVecZ(p, deltaSin, deltaCos);
		rotateVec(p, rx, ry, rz);
		add(p, b);
		writeCoord(coords, index, p);
	}
}


//===========================
//   main stuff
//===========================

typedef struct __align__(8) {
	double xdmin;
	double xdmax;
	double xd;
	double step;
} LineSearchArgs;
// sizeof = 32

typedef struct {
	double *coords;
	const int *dihedralIndices;
	const int numRotatedIndices;
	const int *rotatedIndices;
	const int energyBlockThreads;
	const int energyNumBlocks;
	const int *atomFlags;
	const double *precomputed;
	const int *subsetTable;
	const ForcefieldArgs *ffargs;
	double *blockEnergies;
} PoseAndEnergyArgs;

typedef struct __align__(8) {
	double xdstar;
	double fxdstar;
} Results;
// sizeof = 16


__device__ int calcNumBlocks(int numThreads, int blockThreads) {
	return (numThreads + blockThreads - 1)/blockThreads;
}

__device__ void pose(const PoseAndEnergyArgs &args, double dihedralRadians) {
	pose(
		args.coords,
		args.dihedralIndices,
		args.numRotatedIndices,
		args.rotatedIndices,
		dihedralRadians
	);
}

__device__ double poseAndEnergy(const PoseAndEnergyArgs &args, double dihedralRadians) {
	pose(args, dihedralRadians);
	return energy(
		args.energyNumBlocks,
		args.energyBlockThreads,
		args.coords,
		args.atomFlags,
		args.precomputed,
		args.subsetTable,
		args.ffargs,
		args.blockEnergies
	);
}

__device__ double getTolerance(double f) {
	
	// use full tolerance, unless f is very small
	// then scale by the magnitude of f
	return Tolerance * fmax(1.0, fabs(f));
}

// TEMP
__device__ double toDegrees(double radians) {
	return radians*180/M_PI;
}

extern "C" __global__ void calc(
	double *coords,
	const int *dihedralIndices,
	const int numRotatedIndices,
	const int *rotatedIndices,
	const int *atomFlags,
	const double *precomputed,
	const int *subsetTable,
	const ForcefieldArgs *ffargs,
	const LineSearchArgs *lsargs,
	Results *out
) {

	// calculate energy launch sizes
	int energyBlockThreads = 128;
	int energyNumBlocks = calcNumBlocks(ffargs->numPairs, energyBlockThreads);
	
	// make pose and energy args
	PoseAndEnergyArgs peargs = {
		coords, dihedralIndices, numRotatedIndices, rotatedIndices,
		energyBlockThreads, energyNumBlocks,
		atomFlags, precomputed, subsetTable, ffargs,
		NULL
	};
	
	// allocate global memory for energies
	peargs.blockEnergies = (double *)malloc(peargs.energyNumBlocks*sizeof(double));
	
	// put some vars on the stack
	double xdmin = lsargs->xdmin;
	double xdmax = lsargs->xdmax;
	double xd = lsargs->xd;
	double step = lsargs->step;
	
	double fxdmin = NAN;
	double fxdmax = NAN;
		
	// get initial energy
	double fxd = poseAndEnergy(peargs, xd);
	
	// get the positive (p) neighbor
	double xdp = xd + step;
	double fxdp = INFINITY;
	if (xdp <= xdmax) {
		fxdp = poseAndEnergy(peargs, xdp);
	}
	
	// get the negative (n) neighbor
	double fxdm = INFINITY;
	double xdm = xd - step;
	if (xdm >= xdmin) {
		fxdm = poseAndEnergy(peargs, xdm);
	}
	
	// fit a quadratic to the objective function, locally:
	// q(x) = fx + a*(x - xd)^2 + b*(x - xd)
	// a*step^2 + b*step = fxp - fx
	// a*step^2 - b*step = fxm - fx
	
	// solve for a to determine the shape
	double a = (fxdp + fxdm - 2*fxd)/(2*step*step);
	double xdstar = 0;
	if (a <= 0 || a == NAN || a == INFINITY) {
		
		// negative a means quadratic is concave down, I think
		// infinite or nan a means we're hitting a constraint or impossible conformation
		// so just minimize over the endpoints of the interval
		if (fxdm < fxdp) {
			xdstar = xdm;
		} else {
			xdstar = xdp;
		}
		
	} else {
		
		// positive a means quadratic is concave up, I think
		// solve for the b param
		double b = (fxdp - fxd)/step - a*step;
		
		// then minimize the quadratic to get the minimum x:
		// 2*a*(x - xd) + b = 0
		xdstar = xd - b/2/a;
	}

	// clamp xdstar to the range
	if (xdstar < xdmin) {
		xdstar = xdmin;
	}
	if (xdstar > xdmax) {
		xdstar = xdmax;
	}
	
	double fxdstar = poseAndEnergy(peargs, xdstar);
	
	// TEMP
	int numUpSurfs = 0;
	int numDownSurfs = 0;
	
	/* TEMP
	printf("min:%8.3f,%14.6f   max:%8.3f,%14.6f   init:%8.3f,%14.6f   minus:%8.3f,%14.6f   plus:%8.3f,%14.6f   opt:%8.3f,%14.6f   surfs=[%2d,%2d]\n",
		toDegrees(xdmin), fxdmin,
		toDegrees(xdmax), fxdmax,
		toDegrees(xd), fxd,
		toDegrees(xdm), fxdm,
		toDegrees(xdp), fxdp,
		toDegrees(xdstar), fxdstar,
		numUpSurfs, numDownSurfs
	);
	*/
	
	// did we go downhill?
	if (fxdstar < fxd) {
		
		// surf along f locally to try to find better minimum
		double xdsurfHere = xdstar;
		double fxdsurfHere = fxdstar;
		while (true) {
		
			// TEMP
			numUpSurfs++;
			
			// take a step twice as far as we did last time
			double xdsurfNext = xd + 2*(xdsurfHere - xd);
			
			// did we step off the min?
			if (xdsurfNext < xdmin) {
				
				// if the min is better, go there instead
				if (isnan(fxdmin)) {
					fxdmin = poseAndEnergy(peargs, xdmin);
				}
				if (fxdmin < fxdstar) {
					xdsurfHere = xdmin;
					fxdsurfHere = fxdmin;
				}
				
				break;
			
			// did we step off the max?
			} else if (xdsurfNext > xdmax) {
				
				// if the max is better, go there instead
				if (isnan(fxdmax)) {
					fxdmax = poseAndEnergy(peargs, xdmax);
				}
				if (fxdmax < fxdstar) {
					xdsurfHere = xdmax;
					fxdsurfHere = fxdmax;
				}
				
				break;
			}
			
			double fxdsurfNext = poseAndEnergy(peargs, xdsurfNext);
			
			/* TEMP
			printf("surf:%8.3f,%14.6f   surfs=[%2d,%2d]\n",
				toDegrees(xdsurfNext), fxdsurfNext,
				numUpSurfs, numDownSurfs
			);
			*/
			
			// did we improve the min enough to keep surfing?
			if (fxdsurfNext < fxdsurfHere - getTolerance(fxdsurfHere)) {
			
				// yeah, keep going
				xdsurfHere = xdsurfNext;
				fxdsurfHere = fxdsurfNext;
				
			} else {
				
				// nope, stop surfing
				break;
			}
		}
		
		// update the minimum estimate so far
		xdstar = xdsurfHere;
		fxdstar = fxdsurfHere;
		
	// did we go significantly uphill?
	} else if (fxdstar > fxd + Tolerance) {
		
		// try to surf back downhill
		double xdsurfHere = xdstar;
		double fxdsurfHere = fxdstar;
		while (true) {
		
			// TEMP
			numDownSurfs++;
			
			// cut the step in half
			double xdsurfNext = xd + (xdsurfHere - xd)/2;
			double fxdsurfNext = poseAndEnergy(peargs, xdsurfNext);
			
			// did we improve the min enough to keep surfing?
			if (fxdsurfNext < fxdsurfHere - getTolerance(fxdsurfHere)) {
			
				// yeah, keep going
				xdsurfHere = xdsurfNext;
				fxdsurfHere = fxdsurfNext;
				
			} else {
				
				// nope, stop surfing
				break;
			}
		}
		
		// did the quadratic step help at all?
		if (fxdstar < fxd) {
			
			// yeah, keep it!
			
		} else {
			
			// nope, the original spot was lower
			xdstar = xd;
			fxdstar = fxd;
		}
		
		// did surfing help at all?
		if (fxdsurfHere < fxdstar) {
			
			// yeah, use the surf spot
			xdstar = xdsurfHere;
			fxdstar = fxdsurfHere;
		}
	}
	
	/* TEMP
	printf("min:%8.3f,%14.6f   max:%8.3f,%14.6f   init:%8.3f,%14.6f   minus:%8.3f,%14.6f   plus:%8.3f,%14.6f   opt:%8.3f,%14.6f   surfs=[%2d,%2d]\n",
		toDegrees(xdmin), fxdmin,
		toDegrees(xdmax), fxdmax,
		toDegrees(xd), fxd,
		toDegrees(xdm), fxdm,
		toDegrees(xdp), fxdp,
		toDegrees(xdstar), fxdstar,
		numUpSurfs, numDownSurfs
	);
	*/
	
	// try to jump over walls arbitrarily
	// look in a 1-degree step for a better minimum
	
	// NOTE: skipping this can make minimization a bit faster,
	// but skipping this causes a noticeable rise in final energies too
	// it's best to keep doing it I think
	
	xdm = xdstar - OneDegree;
	if (xdm >= xdmin) {
		fxdm = poseAndEnergy(peargs, xdm);
		if (fxdm < fxdstar) {
			xdstar = xdm;
			fxdstar = fxdm;
		}
	}
	
	xdp = xdstar + OneDegree;
	if (xdp <= xdmax) {
		fxdp = poseAndEnergy(peargs, xdp);
		if (fxdp < fxdstar) {
			xdstar = xdp;
			fxdstar = fxdp;
		}
	}
	
	/* TEMP
	printf("min:%8.3f,%14.6f   max:%8.3f,%14.6f   init:%8.3f,%14.6f   minus:%8.3f,%14.6f   plus:%8.3f,%14.6f   opt:%8.3f,%14.6f   surfs=[%2d,%2d]\n",
		toDegrees(xdmin), fxdmin,
		toDegrees(xdmax), fxdmax,
		toDegrees(xd), fxd,
		toDegrees(xdm), fxdm,
		toDegrees(xdp), fxdp,
		toDegrees(xdstar), fxdstar,
		numUpSurfs, numDownSurfs
	);
	*/
	
	// one last pose
	pose(peargs, xdstar);
		
	// set results
	out->xdstar = xdstar;
	out->fxdstar = fxdstar;
	
	// cleanup
	free(peargs.blockEnergies);
}
