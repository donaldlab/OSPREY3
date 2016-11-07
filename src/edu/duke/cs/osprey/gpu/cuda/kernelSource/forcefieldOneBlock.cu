
// compile with:
// nvcc -fatbin -arch=compute_20 "kernelSource/subForcefields.cu" -o "kernelBinaries/subForcefields.bin"

#include <stdio.h>


typedef struct __align__(8) {
	double coulombFactor;
	double scaledCoulombFactor;
	double solvCutoff2;
	bool useDistDepDielec;
	bool useHEs;
	bool useHVdw;
} ForcefieldArgs;
// sizeof = 32


__device__ int divUp(int a, int b) {
	// ie.,  ceil(a/b)
	return (a + b - 1)/b;
}

__device__ int getAtomIndex(int flags) {
	return abs(flags) - 1;
}

__device__ bool isHydrogen(int flags) {
	return flags > 0;
}

__device__ double calcPairEnergy(
	const double *coords,
	const int *atomFlags,
	const double *precomputed,
	const ForcefieldArgs *args,
	const int i,
	const bool is14Pair
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

__device__ void blockSum(double *scratch, double *out) {

	int threadId = threadIdx.x;
	int numThreads = blockDim.x;
	
	// compute the energy sum in SIMD-style
	// see url for a tutorial on GPU reductions:
	// http://developer.amd.com/resources/articles-whitepapers/opencl-optimization-case-study-simple-reductions/

	__syncthreads();
	
	for (int offset = 1; offset < numThreads; offset <<= 1) {
	
		// sum this level of the reduction tree
		int mask = (offset << 1) - 1;
		if ((threadId & mask) == 0) {
			int pos = threadId + offset;
			if (pos < numThreads) {
				scratch[threadId] += scratch[pos];
			}
		}
		
		__syncthreads();
	}
	
	// finally, if we're the 0 thread, write the summed energy for this work group
	if (threadId == 0) {
		int blockId = blockIdx.x;
		out[blockId] = scratch[0];
	}
}

extern "C" __global__ void calcEnergy(
	const double *coords,
	const int *atomFlags,
	const double *precomputed,
	const ForcefieldArgs *args,
	const int *subsetTable,
	const int numPairs,
	const int num14Pairs,
	double *out
) {

	extern __shared__ double threadEnergies[]; // NOTE: can't declare as pointer, must be array
	
	int threadId = threadIdx.x;
	int numThreads = blockDim.x;
	int blockId = blockIdx.x;
	int numBlocks = gridDim.x;
	
	double energy = 0;
	
	// partition atom pairs among blocks
	// NOTE: keep adjacent pairs in the same block, to improve memory locality
	int blockThreads = divUp(numPairs, numBlocks);
	int firstBlockI = blockId*blockThreads;
	int lastBlockI = min(numPairs, (blockId + 1)*blockThreads);
	
	// partition atom pairs among threads
	// NOTE: stagger atom pairs so we get efficient memory access (all threads share the same cache)
	// ie, thread 1 does pair pair 1, thread 2 does atom pair 2, etc
	int i = firstBlockI + threadId;
	while (i < lastBlockI) {
		energy += calcPairEnergy(coords, atomFlags, precomputed, args, subsetTable[i], i < num14Pairs);
		i += numThreads;
	}
	
	threadEnergies[threadId] = energy;
	
	// sum energies from all threads
	blockSum(threadEnergies, out);
}
