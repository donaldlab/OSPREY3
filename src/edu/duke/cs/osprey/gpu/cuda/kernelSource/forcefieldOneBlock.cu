
// compile with:
// nvcc -fatbin -arch=compute_20 "kernelSource/subForcefields.cu" -o "kernelBinaries/subForcefields.bin"

#include <stdio.h>


typedef struct __align__(8) {
	int numPairs;
	int num14Pairs;
	double coulombFactor;
	double scaledCoulombFactor;
	double solvCutoff2;
	bool useDistDepDielec;
	bool useHEs;
	bool useHVdw;
} ForcefieldArgs;
// sizeof = 40

typedef struct __align__(8) {
	int subsetTableOffset;
	int numPairs;
	int num14Pairs;
	int rotatedIndicesOffset;
	int numRotatedIndices;
	int firstModifiedCoord;
	int lastModifiedCoord;
	// 4 bytes space
} DofArgs;
// sizeof = 32


typedef struct {
	int threadId;
	int numThreads;
	int blockId;
	int numBlocks;
	int blockThreads;
	int blockFirstPair;
	int blockLastPair;
	int threadFirstPair;
	int threadStride;
} Partition;

__device__ int divUp(int a, int b) {
	// ie.,  ceil(a/b)
	return (a + b - 1)/b;
}

__device__ void makePartition(Partition &p, int numPairs) {

	p.threadId = threadIdx.x;
	p.numThreads = blockDim.x;
	p.blockId = blockIdx.x;
	p.numBlocks = gridDim.x;
	
	// partition atom pairs among blocks
	// NOTE: keep adjacent pairs in the same block, to improve memory locality
	p.blockThreads = divUp(numPairs, p.numBlocks);
	p.blockFirstPair = p.blockId*p.blockThreads;
	p.blockLastPair = min(numPairs, (p.blockId + 1)*p.blockThreads);
	
	// partition atom pairs among threads
	// NOTE: stagger atom pairs so we get efficient memory access (all threads share the same cache)
	// ie, thread 1 does pair pair 1, thread 2 does atom pair 2, etc
	p.threadFirstPair = p.blockFirstPair + p.threadId;
	p.threadStride = p.numThreads;
}

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
	
	return energy;
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
	__syncthreads();
}

__device__ void blockSum(double threadEnergy, double *threadEnergies, double *blockEnergies) {

	int threadId = threadIdx.x;
	int numThreads = blockDim.x;
	
	// compute the energy sum in SIMD-style
	// see url for a tutorial on GPU reductions:
	// http://developer.amd.com/resources/articles-whitepapers/opencl-optimization-case-study-simple-reductions/

	threadEnergies[threadId] = threadEnergy;
	
	__syncthreads();
	
	for (int offset = 1; offset < numThreads; offset <<= 1) {
	
		// sum this level of the reduction tree
		int mask = (offset << 1) - 1;
		if ((threadId & mask) == 0) {
			int pos = threadId + offset;
			if (pos < numThreads) {
				threadEnergies[threadId] += threadEnergies[pos];
			}
		}
		
		__syncthreads();
	}
	
	// finally, if we're the 0 thread, write the summed energy for this work group
	if (threadId == 0) {
		blockEnergies[blockIdx.x] = threadEnergies[0];
	}
}

// NOTE: can't declare as pointer, must be array
extern __shared__ unsigned char shared[];


extern "C" __global__ void calcEnergy(
	const double *coords,
	const int *atomFlags,
	const double *precomputed,
	const ForcefieldArgs *ffargs,
	double *out
) {

	// partition work among blocks/threads
	Partition p;
	makePartition(p, ffargs->numPairs);
	
	// partition shared memory
	double *threadEnergies = (double *)shared;
	
	// add up the pairwise energies
	double energy = 0;
	for (int i = p.threadFirstPair; i < p.blockLastPair; i += p.threadStride) {
		energy += calcPairEnergy(
			coords, atomFlags, precomputed, ffargs,
			i, i < ffargs->num14Pairs,
			NULL, 0, 0
		);
	}
	blockSum(energy, threadEnergies, out);
}

extern "C" __global__ void calcDofEnergy(
	const double *coords,
	const int *atomFlags,
	const double *precomputed,
	const ForcefieldArgs *ffargs,
	const int *subsetTables,
	const DofArgs *dofargs,
	const int d,
	double *out
) {

	// get our dof
	const DofArgs *dofdargs = dofargs + d;
	const int *subsetTable = subsetTables + dofdargs->subsetTableOffset;
	
	// partition work among blocks/threads
	Partition p;
	makePartition(p, dofdargs->numPairs);
	
	// partition shared memory
	double *threadEnergies = (double *)shared;
	
	// add up the pairwise energies
	double energy = 0;
	for (int i = p.threadFirstPair; i < p.blockLastPair; i += p.threadStride) {
		energy += calcPairEnergy(
			coords, atomFlags, precomputed, ffargs,
			subsetTable[i], i < dofdargs->num14Pairs,
			NULL, 0, 0
		);
	}
	blockSum(energy, threadEnergies, out);
}

extern "C" __global__ void poseAndCalcDofEnergy(
	const double *coords,
	const int *atomFlags,
	const double *precomputed,
	const ForcefieldArgs *ffargs,
	const int *subsetTables,
	const int *dihedralIndices,
	const int *rotatedIndices,
	const DofArgs *dofargs,
	const int d,
	const double dihedralRadians,
	double *out
) {

	// get our dof
	const DofArgs *dofdargs = dofargs + d;
	const int *subsetTable = subsetTables + dofdargs->subsetTableOffset;
	
	// partition work among blocks/threads
	Partition p;
	makePartition(p, dofdargs->numPairs);
	
	// partition shared memory
	double *threadEnergies = (double *)shared;
	double *modifiedCoords = (double *)(shared + p.numThreads*sizeof(double));
	
	// copy the coords we need to modify to shared mem
	int numModifiedCoords = dofdargs->lastModifiedCoord - dofdargs->firstModifiedCoord + 1;
	for (int i = p.threadId; i < numModifiedCoords; i += p.blockThreads) {
		modifiedCoords[i] = coords[dofdargs->firstModifiedCoord + i];
	}
	__syncthreads();
	
	// pose the coords
	const int *dihedralIndicesD = dihedralIndices + d;
	const int *rotatedIndicesD = rotatedIndices + d;
	pose(modifiedCoords, dihedralIndicesD, dofdargs->numRotatedIndices, rotatedIndicesD, dihedralRadians);
	
	// add up the pairwise energies
	double energy = 0;
	for (int i = p.threadFirstPair; i < p.blockLastPair; i += p.threadStride) {
		energy += calcPairEnergy(
			coords, atomFlags, precomputed, ffargs,
			subsetTable[i], i < dofdargs->num14Pairs,
			modifiedCoords, dofdargs->firstModifiedCoord, dofdargs->lastModifiedCoord
		);
	}
	blockSum(energy, threadEnergies, out);
}

