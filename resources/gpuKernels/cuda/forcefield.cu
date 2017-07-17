
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
	"kernelSource/forcefield.cu" -o "kernelBinaries/forcefield.bin"
	
	See Maxwell compatibility guide for more info:
	http://docs.nvidia.com/cuda/maxwell-compatibility-guide/index.html#building-maxwell-compatible-apps-using-cuda-6-0
*/


typedef struct __align__(8) {
	int numPairs; // @ 0
	int num14Pairs; // @ 4
	double coulombFactor; // @ 8
	double scaledCoulombFactor; // @ 16
	double solvCutoff2; // @ 24
	bool useDistDepDielec; // @ 32
	bool useHEs; // @ 33
	bool useHVdw; // @ 34
	bool useSubset; // @ 35
	bool useEEF1; // @ 36
	// 3 bytes pad
} ForcefieldArgs;
// sizeof = 40


__device__ int getAtomIndex(int flags) {
	return abs(flags) - 1;
}

__device__ bool isHydrogen(int flags) {
	return flags > 0;
}

extern "C" __global__ void calc(
	const double *coords,
	const int *atomFlags,
	const double *precomputed,
	const int *subsetTable,
	const ForcefieldArgs *args,
	double *out
) {

	extern __shared__ double scratch[];

	// start with zero energy
	double energy = 0;
	
	int globalId = blockIdx.x*blockDim.x + threadIdx.x;
	
	// which atom pair are we calculating?
	if (globalId < args->numPairs) {
	
		int i = globalId;
		
		// are we using the subset?
		if (args->useSubset) {
			i = subsetTable[i];
		}
		
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
		out[blockIdx.x] = scratch[0];
	}
}
