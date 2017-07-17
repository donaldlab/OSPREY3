
/* compile with:

nvcc -fatbin -O2
	-gencode=arch=compute_20,code=sm_20
	-gencode=arch=compute_30,code=sm_30
	-gencode=arch=compute_35,code=sm_35
	-gencode=arch=compute_50,code=sm_50
	-gencode=arch=compute_52,code=sm_52
	-gencode=arch=compute_60,code=sm_60
	-gencode=arch=compute_61,code=sm_61
	-gencode=arch=compute_62,code=sm_62
	-gencode=arch=compute_62,code=compute_62
	"kernelSource/ccd.cu" -o "kernelBinaries/ccd.bin"
	
	See Maxwell compatibility guide for more info:
	http://docs.nvidia.com/cuda/maxwell-compatibility-guide/index.html#building-maxwell-compatible-apps-using-cuda-6-0
*/


// use same settings as CCDMinimizer on the java side
const int MaxIterations = 30;
const double ConvergenceThreshold = 0.001;
const double Tolerance = 1e-6;
const double OneDegree = M_PI/180.0;
const double InitialStep = OneDegree*0.25;


typedef struct __align__(8) {
	int numPairs;
	int num14Pairs;
	double coulombFactor;
	double scaledCoulombFactor;
	double solvCutoff2;
	double internalSolvationEnergy;
	bool useDistDepDielec;
	bool useHEs;
	bool useHVdw;
	bool useEEF1;
} ForcefieldArgs;
// sizeof = 48

typedef struct __align__(8) {
	int subsetTableOffset;
	int numPairs;
	int num14Pairs;
	int rotatedIndicesOffset;
	int numRotatedIndices;
	int firstModifiedCoord;
	int lastModifiedCoord;
	// 4 bytes space
	double internalSolvationEnergy;
} DofArgs;
// sizeof = 48

typedef struct {
	const double *coords;
	const int *atomFlags;
	const double *precomputed;
	const ForcefieldArgs *ffargs;
	const int *subsetTable;
	const int *dihedralIndices;
	const int *rotatedIndices;
	const DofArgs *dofdargs;
	double *modifiedCoords;
	double *threadEnergies;
} DofPoseAndEnergyArgs;

typedef struct {
	double step;
	double xdstar;
	double fxdstar;
} LinesearchOut;

typedef struct __align__(8) {
	double xd;
	double xdmin;
	double xdmax;
} XdAndBounds;
// sizeof = 24


// dayum, CUDA... no libraries for vector math? what gives??

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

__device__ int getAtomIndex(int flags) {
	return abs(flags) - 1;
}

__device__ bool isHydrogen(int flags) {
	return flags > 0;
}

__device__ double getCoord(const double *coords, const double *modifiedCoords, int firstModifiedCoord, int lastModifiedCoord, int coordIndex) {
	if (modifiedCoords != NULL && coordIndex >= firstModifiedCoord && coordIndex <= lastModifiedCoord) {
		return modifiedCoords[coordIndex - firstModifiedCoord];
	} else {
		return coords[coordIndex];
	}
}

__device__ double calcPairEnergy(
	const double *coords,
	const int *atomFlags,
	const double *precomputed,
	const ForcefieldArgs *args,
	const int i,
	const bool is14Pair,
	const double *modifiedCoords,
	const int firstModifiedCoord,
	const int lastModifiedCoord
) {

	// start with zero energy
	double energy = 0;
	
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
		d = getCoord(coords, modifiedCoords, firstModifiedCoord, lastModifiedCoord, atom1Index3)
			- getCoord(coords, modifiedCoords, firstModifiedCoord, lastModifiedCoord, atom2Index3);
		r2 += d*d;
		d = getCoord(coords, modifiedCoords, firstModifiedCoord, lastModifiedCoord, atom1Index3 + 1)
			- getCoord(coords, modifiedCoords, firstModifiedCoord, lastModifiedCoord, atom2Index3 + 1);
		r2 += d*d;
		d = getCoord(coords, modifiedCoords, firstModifiedCoord, lastModifiedCoord, atom1Index3 + 2)
			- getCoord(coords, modifiedCoords, firstModifiedCoord, lastModifiedCoord, atom2Index3 + 2);
		r2 += d*d;
	}
	
	int i9 = i*9;
	
	// calculate electrostatics
	if (bothHeavy || args->useHEs) {
	
		double esEnergy = is14Pair ? args->scaledCoulombFactor : args->coulombFactor;
		
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
	if (args->useEEF1 && bothHeavy && r2 < args->solvCutoff2) {
			
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
	
	return energy;
}

__device__ void blockSum(double threadEnergy, double *threadEnergies) {

	// compute the energy sum in SIMD-style
	// see url for a tutorial on GPU reductions:
	// http://developer.amd.com/resources/articles-whitepapers/opencl-optimization-case-study-simple-reductions/

	threadEnergies[threadIdx.x] = threadEnergy;
	
	__syncthreads();
	
	for (int offset = 1; offset < blockDim.x; offset <<= 1) {
	
		// sum this level of the reduction tree
		int mask = (offset << 1) - 1;
		if ((threadIdx.x & mask) == 0) {
			int pos = threadIdx.x + offset;
			if (pos < blockDim.x) {
				threadEnergies[threadIdx.x] += threadEnergies[pos];
			}
		}
		
		__syncthreads();
	}
}

__device__ double calcFullEnergy(const double *coords, const int *atomFlags, const double *precomputed, const ForcefieldArgs *ffargs, double *threadEnergies) {

	// add up the pairwise energies
	double energy = 0;
	for (int i = threadIdx.x; i < ffargs->numPairs; i += blockDim.x) {
		energy += calcPairEnergy(
			coords, atomFlags, precomputed, ffargs,
			i, i < ffargs->num14Pairs,
			NULL, 0, 0
		);
	}
	blockSum(energy, threadEnergies);
	energy = threadEnergies[0];
	
	if (ffargs->useEEF1) {
		energy += ffargs->internalSolvationEnergy;
	}
	
	return energy;
}

__device__ void pose(const DofPoseAndEnergyArgs &args, double dihedralRadians) {

	double *coords = args.modifiedCoords;

	// get the four atom positions: a, b, c, d
	double3 a = readCoord(coords, args.dihedralIndices[0]);
	double3 b = readCoord(coords, args.dihedralIndices[1]);
	double3 c = readCoord(coords, args.dihedralIndices[2]);
	double3 d = readCoord(coords, args.dihedralIndices[3]);
	
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
	
	// modify the atoms in parallel
	if (threadIdx.x < args.dofdargs->numRotatedIndices) {
		int index = args.rotatedIndices[threadIdx.x]; 
		double3 p = readCoord(coords, index);
		sub(p, b);
		rotateVecInverse(p, rx, ry, rz);
		rotateVecZ(p, deltaSin, deltaCos);
		rotateVec(p, rx, ry, rz);
		add(p, b);
		writeCoord(coords, index, p);
	}
	__syncthreads();
}

__device__ void copyCoordsGtoS(const DofArgs *dofdargs, const double *coords, double *modifiedCoords) {

	int numModifiedCoords = dofdargs->lastModifiedCoord - dofdargs->firstModifiedCoord + 1;
	
	for (int i = threadIdx.x; i < numModifiedCoords; i += blockDim.x) {
		modifiedCoords[i] = coords[dofdargs->firstModifiedCoord + i];
	}
	__syncthreads();
}

__device__ void copyCoordsStoG(const DofArgs *dofdargs, double *coords, const double *modifiedCoords) {

	int numModifiedCoords = dofdargs->lastModifiedCoord - dofdargs->firstModifiedCoord + 1;
	
	for (int i = threadIdx.x; i < numModifiedCoords; i += blockDim.x) {
		coords[dofdargs->firstModifiedCoord + i] = modifiedCoords[i];
	}
	__syncthreads();
}

__device__ double poseAndCalcDofEnergy(const DofPoseAndEnergyArgs &args, double dihedralRadians) {
	
	// pose the coords
	pose(args, dihedralRadians);
	
	// add up the pairwise energies
	double energy = 0;
	for (int i = threadIdx.x; i < args.dofdargs->numPairs; i += blockDim.x) {
		bool is14Pair = i < args.dofdargs->num14Pairs;
		energy += calcPairEnergy(
			args.coords, args.atomFlags, args.precomputed, args.ffargs,
			args.subsetTable[i], is14Pair,
			args.modifiedCoords, args.dofdargs->firstModifiedCoord, args.dofdargs->lastModifiedCoord
		);
	}
	
	blockSum(energy, args.threadEnergies);
	
	return args.threadEnergies[0] + args.dofdargs->internalSolvationEnergy;
}

__device__ double getTolerance(double f) {
	
	// scale abs(f) by tolerance, unless f is very small
	return Tolerance * fmax(1.0, fabs(f));
}

__device__ LinesearchOut linesearch(const DofPoseAndEnergyArgs &args, const double xd, const double xdmin, const double xdmax, const double step) {

	// sample energy at the starting point
	double fxd = poseAndCalcDofEnergy(args, xd);
	
	// get the positive (p) neighbor
	double fxdp = INFINITY;
	{
		double xdp = xd + step;
		if (xdp <= xdmax) {
			fxdp = poseAndCalcDofEnergy(args, xdp);
		}
	}
	
	// get the negative (n) neighbor
	double fxdm = INFINITY;
	{
		double xdm = xd - step;
		if (xdm >= xdmin) {
			fxdm = poseAndCalcDofEnergy(args, xdm);
		}
	}
	
	// fit a quadratic to the objective function, locally:
	// q(x) = fx + a*(x - xd)^2 + b*(x - xd)
	// a*step^2 + b*step = fxp - fx
	// a*step^2 - b*step = fxm - fx
	
	// solve for a to determine the shape
	double xdstar = 0;
	{
		double shape = fxdp + fxdm - 2*fxd;
		const double ShapeEpsilon = 1e-12;
		if (shape < -ShapeEpsilon || shape == NAN || shape == INFINITY) {
			
			// negative shape means quadratic is concave down
			// infinite or nan a means we're hitting a constraint or impossible conformation
			// so just minimize over the endpoints of the interval
			if (fxdm < fxdp) {
				xdstar = xd - step;
			} else {
				xdstar = xd + step;
			}
		
		} else if (shape <= ShapeEpsilon) {
		
			// flat here, don't step
			xdstar = xd;
			
		} else {
			
			// positive shape means quadratic is concave up
			// step to the optimum
			xdstar = xd + (fxdm - fxdp)*step/2/shape;
		}
	}

	// clamp xdstar to the range
	if (xdstar < xdmin) {
		xdstar = xdmin;
	}
	if (xdstar > xdmax) {
		xdstar = xdmax;
	}
	
	double fxdstar = poseAndCalcDofEnergy(args, xdstar);
	
	// did we go downhill?
	if (fxdstar < fxd) {
	
		double fxdmin = NAN;
		double fxdmax = NAN;
		
		// surf along f locally to try to find better minimum
		double xdsurfHere = xdstar;
		double fxdsurfHere = fxdstar;
		while (true) {
		
			// take a step twice as far as we did last time
			double xdsurfNext = xd + 2*(xdsurfHere - xd);
			
			// did we step off the min?
			if (xdsurfNext < xdmin) {
				
				// if the min is better, go there instead
				if (isnan(fxdmin)) {
					fxdmin = poseAndCalcDofEnergy(args, xdmin);
				}
				if (fxdmin < fxdsurfHere) {
					xdsurfHere = xdmin;
					fxdsurfHere = fxdmin;
				}
				
				break;
			
			// did we step off the max?
			} else if (xdsurfNext > xdmax) {
				
				// if the max is better, go there instead
				if (isnan(fxdmax)) {
					fxdmax = poseAndCalcDofEnergy(args, xdmax);
				}
				if (fxdmax < fxdsurfHere) {
					xdsurfHere = xdmax;
					fxdsurfHere = fxdmax;
				}
				
				break;
			}
			
			double fxdsurfNext = poseAndCalcDofEnergy(args, xdsurfNext);
			
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
		
			// cut the step in half
			double xdsurfNext = xd + (xdsurfHere - xd)/2;
			double fxdsurfNext = poseAndCalcDofEnergy(args, xdsurfNext);
			
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
	
	// compute the step taken before wall jumping
	LinesearchOut out;
	out.step = xdstar - xd;
	
	// try to jump over walls arbitrarily
	// look in a 1-degree step for a better minimum
	
	// NOTE: skipping this can make minimization a bit faster,
	// but skipping this causes a noticeable rise in final energies too
	// it's best to keep doing it I think
	
	double xdm = xdstar - OneDegree;
	double xdp = xdstar + OneDegree;
	if (xdm >= xdmin) {
		fxdm = poseAndCalcDofEnergy(args, xdm);
		if (fxdm < fxdstar) {
			xdstar = xdm;
			fxdstar = fxdm;
		}
	}
	if (xdp <= xdmax) {
		fxdp = poseAndCalcDofEnergy(args, xdp);
		if (fxdp < fxdstar) {
			xdstar = xdp;
			fxdstar = fxdp;
		}
	}
	
	// one last pose
	pose(args, xdstar);
	
	// set outputs
	out.xdstar = xdstar;
	out.fxdstar = fxdstar;
	return out;
}

__device__ void copyx(const double *src, double *dest, int size) {
	for (int i = threadIdx.x; i < size; i += blockDim.x) {
		dest[i] = src[i];
	}
	__syncthreads();
}

extern "C" __global__ void ccd(
	double *coords,
	const int *atomFlags,
	const double *precomputed,
	const ForcefieldArgs *ffargs,
	const int *subsetTables,
	const int *dihedralIndices,
	const int *rotatedIndices,
	const DofArgs *dofargs,
	const int maxNumModifiedCoords,
	const XdAndBounds *xAndBounds,
	const int numDofs,
	double *out // size is numDofs + 1
) {

	// partition shared memory
	extern __shared__ unsigned char shared[];
	double *threadEnergies = (double *)shared;
	double *modifiedCoords = threadEnergies + blockDim.x;
	double *nextx = modifiedCoords + maxNumModifiedCoords;
	double *firstSteps = nextx + numDofs;
	double *lastSteps = firstSteps + numDofs;
	
	// partition out memory
	double *outfx = out;
	double *outx = out + 1; // size is numDofs
	
	// build the poseAndCalcDofEnergy() args
	// at least, the parts that are independent of the dof
	DofPoseAndEnergyArgs args = {
		coords, atomFlags, precomputed, ffargs,
		NULL, NULL, NULL, NULL,
		modifiedCoords, threadEnergies
	};
	
	// make a copy of x (in parallel)
	double *herex = outx;
	for (int d = threadIdx.x; d < numDofs; d += blockDim.x) {
		herex[d] = xAndBounds[d].xd;
	}
	__syncthreads();
	
	// init the step sizes
	for (int d = threadIdx.x; d < numDofs; d+= blockDim.x) {
		firstSteps[d] = OneDegree;
		lastSteps[d] = OneDegree;
	}
	__syncthreads();
	
	// get the initial energy
	double herefx = calcFullEnergy(coords, atomFlags, precomputed, ffargs, threadEnergies);
	
	for (int iter=0; iter<MaxIterations; iter++) {
	
		copyx(herex, nextx, numDofs);
		
		// for each dimension...
		for (int d=0; d<numDofs; d++) {
		
			// get the dof info
			args.dofdargs = dofargs + d;
			args.subsetTable = subsetTables + args.dofdargs->subsetTableOffset;
			args.dihedralIndices = dihedralIndices + d*4;
			args.rotatedIndices = rotatedIndices + args.dofdargs->rotatedIndicesOffset;
			
			// copy the coords we need to modify to shared mem
			copyCoordsGtoS(args.dofdargs, coords, modifiedCoords);

			double xd = nextx[d];
			double xdmin = xAndBounds[d].xdmin;
			double xdmax = xAndBounds[d].xdmax;
			
			// get the step size, try to make it adaptive (based on historical steps if possible; else on step #)
			double step;
			{
				double firstStep = firstSteps[d];
				double lastStep = lastSteps[d];
				if (fabs(lastStep) > Tolerance && fabs(firstStep) > Tolerance) {
					step = InitialStep*fabs(lastStep/firstStep);
				} else {
					step = InitialStep/pow(iter + 1.0, 3.0);
				}
				
				// make sure the step isn't so big that the quadratic approximation is worthless
				while (xdmax > xdmin && xd - step < xdmin && xd + step > xdmax) {
					step /= 2;
				}
			}
			
			// do line search
			LinesearchOut lsout = linesearch(args, xd, xdmin, xdmax, step);
			
			// update x and the step
			if (threadIdx.x == 0) {
			
				// update step tracking
				if (iter == 0) {
					firstSteps[d] = lsout.step;
				}
				lastSteps[d] = lsout.step;
				
				// update nextxd
				nextx[d] = lsout.xdstar;
			}
			__syncthreads();
	
			// update the global protein pose
			copyCoordsStoG(args.dofdargs, coords, modifiedCoords);
		}
		
		// evaluate the whole energy function
		double nextfx = calcFullEnergy(coords, atomFlags, precomputed, ffargs, threadEnergies);
		double improvement = herefx - nextfx;
		
		if (improvement > 0) {
		
			// take the step
			copyx(nextx, herex, numDofs);
			herefx = nextfx;
			
			if (improvement < ConvergenceThreshold) {
				break;
			}
			
		} else {
			break;
		}
	}
	
	// update outputs
	*outfx = herefx;
}
