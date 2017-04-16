
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
	"kernelSource/residueForcefield.cu" -o "kernelBinaries/residueForcefield.bin"
	
	See Maxwell compatibility guide for more info:
	http://docs.nvidia.com/cuda/maxwell-compatibility-guide/index.html#building-maxwell-compatible-apps-using-cuda-6-0
*/


typedef struct __align__(8) {
	long flags; // useHEs, useHvdW, distDepDielect, useEEF1
	double coulombFactor;
	double scaledCoulombFactor;
	long offsets[];
} Header;

typedef struct __align__(8) {
	long numAtomPairs;
	long offset1;
	long offset2;
	double weight;
	double offset;
} ResPair;

typedef struct __align__(8) {
	short offset1;
	short offset2;
	int flags; // isHeavyPair, is14Bonded
	double charge;
	double Aij;
	double Bij;
	double radius1;
	double lambda1;
	double alpha1;
	double radius2;
	double lambda2;
	double alpha2;
} AtomPair;


const double SolvCutoff = 9.0;
const double SolvCutoff2 = SolvCutoff*SolvCutoff;


typedef unsigned char byte;


extern "C" __global__ void calc(
	const double *coords,
	const byte *data,
	const int numIndices,
	const int *indices,
	double *out
) {

	extern __shared__ double scratch[];
	
	const Header &header = *(Header *)&data[0];
	
	// unpack the flags
	bool useEEF1, distDepDielect, useHvdW, useHEs;
	{
		long flags = header.flags;
		useEEF1 = (flags & 0x1) == 0x1;
		flags >>= 1;
		distDepDielect = (flags & 0x1) == 0x1;
		flags >>= 1;
		useHvdW = (flags & 0x1) == 0x1;
		flags >>= 1;
		useHEs = (flags & 0x1) == 0x1;
	}
	
	// start with 0 energy
	double energy = 0;
	
	int indexi = 0;
	int numAtomPairs = 0;
	for (int n=0; true; n++) {
	
		int atomPairIndex = blockDim.x*n + threadIdx.x;
	
		// find our res pair and atom pair
		ResPair *resPair = NULL;
		AtomPair *atomPair = NULL;
		int atomPairIndexInResPair;
		for (; indexi<numIndices; indexi++) {
		
			// seek to the res pair
			int offset = header.offsets[indices[indexi]];
			resPair = (ResPair *)&data[offset];
			offset += sizeof(ResPair);
			
			// is our atom pair in this res pair?
			atomPairIndexInResPair = atomPairIndex - numAtomPairs;
			if (atomPairIndexInResPair < resPair->numAtomPairs) {
			
				// yup, found it!
				atomPair = (AtomPair *)&data[offset + atomPairIndexInResPair*sizeof(AtomPair)];
				break;	
			}
			
			numAtomPairs += resPair->numAtomPairs;
		}
		
		// stop if we ran out of atom pairs
		if (atomPair == NULL) {
			break;
		}
		
		double resPairEnergy = 0;
		
		// get the radius
		double r2, r;
		{
			int offset1 = resPair->offset1 + atomPair->offset1;
			int offset2 = resPair->offset2 + atomPair->offset2;
			double d = coords[offset1] - coords[offset2];
			r2 = d*d;
			d = coords[offset1 + 1] - coords[offset2 + 1];
			r2 += d*d;
			d = coords[offset1 + 2] - coords[offset2 + 2];
			r2 += d*d;
			r = sqrt(r2);
		}
		
		// unpack the flags
		bool isHeavyPair, is14Bonded;
		{
			int flags = atomPair->flags;
			isHeavyPair = (flags & 0x1) == 0x1;
			flags >>= 1;
			is14Bonded = (flags & 0x1) == 0x1;
		}
		
		// electrostatics
		if (isHeavyPair || useHEs) {
			if (is14Bonded) {
				if (distDepDielect) {
					resPairEnergy += header.scaledCoulombFactor*atomPair->charge/r2;
				} else {
					resPairEnergy += header.scaledCoulombFactor*atomPair->charge/r;
				}
			} else {
				if (distDepDielect) {
					resPairEnergy += header.coulombFactor*atomPair->charge/r2;
				} else {
					resPairEnergy += header.coulombFactor*atomPair->charge/r;
				}
			}
		}
		
		// van der Waals
		if (isHeavyPair || useHvdW) {
			double r6 = r2*r2*r2;
			double r12 = r6*r6;
			resPairEnergy += atomPair->Aij/r12 - atomPair->Bij/r6;
		}
		
		// solvation
		if (useEEF1 && isHeavyPair && r2 < SolvCutoff2) {
			double Xij = (r - atomPair->radius1)/atomPair->lambda1;
			double Xji = (r - atomPair->radius2)/atomPair->lambda2;
			resPairEnergy -= (atomPair->alpha1*exp(-Xij*Xij) + atomPair->alpha2*exp(-Xji*Xji))/r2;
		}
		
		// apply weight and offset
		energy += resPairEnergy*resPair->weight;
		if (atomPairIndexInResPair == 0) {
			energy += resPair->offset;
		}
	}
	
	// compute the energy sum in SIMD-style
	// see url for a tutorial on GPU reductions:
	// http://developer.amd.com/resources/articles-whitepapers/opencl-optimization-case-study-simple-reductions/

	scratch[threadIdx.x] = energy;
	
	__syncthreads();
	
	for (int offset = 1; offset < blockDim.x; offset <<= 1) {
	
		// sum this level of the reduction tree
		int mask = (offset << 1) - 1;
		if ((threadIdx.x & mask) == 0) {
			scratch[threadIdx.x] += scratch[threadIdx.x + offset];
		}
		
		__syncthreads();
	}
	
	// finally, if we're the 0 thread, write the summed energy for this work group
	if (threadIdx.x == 0) {
		*out = scratch[0];
	}
}
