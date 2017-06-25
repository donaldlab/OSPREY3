
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
	unsigned long flags; // useHEs, useHvdW, distDepDielect, useEEF1
	double coulombFactor;
	double scaledCoulombFactor;
} Header;

typedef struct __align__(8) {
	unsigned long numAtomPairs;
	unsigned long offset1;
	unsigned long offset2;
	double weight;
	double offset;
} ResPair;

/* res pairs also have atom pairs after them in struct-of-arrays layout:
typedef struct __align__(8) {

	unsigned long flags; // bit isHeavyPair, bit is14Bonded, 6 bits space, 3 bytes space, short offset1, short offset2
	double charge;
	double Aij;
	double Bij;
	
	// if EEF1 == true
	double radius1;
	double lambda1;
	double alpha1;
	double radius2;
	double lambda2;
	double alpha2;
	
} AtomPair;
*/


typedef unsigned char byte;


class Data {
public:

	const Header & header;
	
	__device__ Data(const byte * const data) :
		m_data(data),
		header(*(Header *)&m_data[0]),
		m_resPairOffsets((unsigned long *)&m_data[sizeof(Header)])
	{
		// nothing else to do
	}
	
	__device__ const ResPair & getResPair(const unsigned int i) const {
		return *(ResPair *)&m_data[m_resPairOffsets[i]];
	}
	
	__device__ unsigned long getAtomPairFlags(const ResPair & resPair, const unsigned int i) const {
		const byte * const p = (byte *)&resPair;
		const unsigned long * const a = (unsigned long *)&p[sizeof(ResPair)];
		return a[i];
	}
	
	__device__ double getAtomPairCharge(const ResPair & resPair, const unsigned int i) const {
		return getAtomPairDouble(resPair, i, 0);
	}
	
	__device__ double getAtomPairAij(const ResPair & resPair, const unsigned int i) const {
		return getAtomPairDouble(resPair, i, 1);
	}
	
	__device__ double getAtomPairBij(const ResPair & resPair, const unsigned int i) const {
		return getAtomPairDouble(resPair, i, 2);
	}
	
	__device__ double getAtomPairRadius1(const ResPair & resPair, const unsigned int i) const {
		return getAtomPairDouble(resPair, i, 3);
	}
	
	__device__ double getAtomPairLambda1(const ResPair & resPair, const unsigned int i) const {
		return getAtomPairDouble(resPair, i, 4);
	}
	
	__device__ double getAtomPairAlpha1(const ResPair & resPair, const unsigned int i) const {
		return getAtomPairDouble(resPair, i, 5);
	}
	
	__device__ double getAtomPairRadius2(const ResPair & resPair, const unsigned int i) const {
		return getAtomPairDouble(resPair, i, 6);
	}
	
	__device__ double getAtomPairLambda2(const ResPair & resPair, const unsigned int i) const {
		return getAtomPairDouble(resPair, i, 7);
	}
	
	__device__ double getAtomPairAlpha2(const ResPair & resPair, const unsigned int i) const {
		return getAtomPairDouble(resPair, i, 8);
	}
	
private:
	const byte * const m_data;
	const unsigned long * const m_resPairOffsets;
	
	__device__ double getAtomPairDouble(const ResPair & resPair, const unsigned int i, const unsigned int pos) const {
		const byte * const p = (byte *)&resPair;
		const double * const a = (double *)&p[sizeof(ResPair) + sizeof(long)*resPair.numAtomPairs + sizeof(double)*pos*resPair.numAtomPairs];
		return a[i];
	}
};

const double SolvCutoff = 9.0;
const double SolvCutoff2 = SolvCutoff*SolvCutoff;


extern "C" __global__ void calc(
	const double *coords,
	const byte *rawdata,
	const int numIndices,
	const int *indices,
	double *out
) {

	// parse data
	const Data data(rawdata);
	
	// unpack the flags
	bool useEEF1, distDepDielect, useHvdW, useHEs;
	{
		long flags = data.header.flags;
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
	
	unsigned int indexi = 0;
	unsigned int numAtomPairs = 0;
	for (int n = threadIdx.x; true; n += blockDim.x) {
	
		// find our res pair and atom pair
		const ResPair *resPair;
		int atomPairIndex = -1;
		
		for (; indexi<numIndices; indexi++) {
		
			unsigned int resPairIndex = indices[indexi];
			resPair = &data.getResPair(resPairIndex);
			
			// is our atom pair in this res pair?
			atomPairIndex = n - numAtomPairs;
			if (atomPairIndex < resPair->numAtomPairs) {
			
				// yup, found it!
				break;	
			}
			
			numAtomPairs += resPair->numAtomPairs;
			atomPairIndex = -1;
		}
		
		// stop if we ran out of atom pairs
		if (atomPairIndex < 0) {
			break;
		}
		
		// unpack the flags
		bool isHeavyPair, is14Bonded;
		unsigned long offset1, offset2;
		{
			unsigned long flags = data.getAtomPairFlags(*resPair, atomPairIndex);
			offset2 = (flags & 0xffff) + resPair->offset2;
			flags >>= 16;
			offset1 = (flags & 0xffff) + resPair->offset1;
			flags >>= 46;
			isHeavyPair = (flags & 0x1) == 0x1;
			flags >>= 1;
			is14Bonded = (flags & 0x1) == 0x1;
		}
		
		double resPairEnergy = 0;
		
		// get the radius
		double r2, r;
		{
			double d = coords[offset1] - coords[offset2];
			r2 = d*d;
			d = coords[offset1 + 1] - coords[offset2 + 1];
			r2 += d*d;
			d = coords[offset1 + 2] - coords[offset2 + 2];
			r2 += d*d;
			r = sqrt(r2);
		}
		
		// electrostatics
		if (isHeavyPair || useHEs) {
		
			double charge = data.getAtomPairCharge(*resPair, atomPairIndex);
			
			if (is14Bonded) {
				if (distDepDielect) {
					resPairEnergy += data.header.scaledCoulombFactor*charge/r2;
				} else {
					resPairEnergy += data.header.scaledCoulombFactor*charge/r;
				}
			} else {
				if (distDepDielect) {
					resPairEnergy += data.header.coulombFactor*charge/r2;
				} else {
					resPairEnergy += data.header.coulombFactor*charge/r;
				}
			}
		}
		
		// van der Waals
		if (isHeavyPair || useHvdW) {
		
			double Aij = data.getAtomPairAij(*resPair, atomPairIndex);
			double Bij = data.getAtomPairBij(*resPair, atomPairIndex);
			
			double r6 = r2*r2*r2;
			double r12 = r6*r6;
			resPairEnergy += Aij/r12 - Bij/r6;
		}
		
		// solvation
		if (useEEF1 && isHeavyPair && r2 < SolvCutoff2) {
		
			double radius1 = data.getAtomPairRadius1(*resPair, atomPairIndex);
			double lambda1 = data.getAtomPairLambda1(*resPair, atomPairIndex);
			double alpha1 = data.getAtomPairAlpha1(*resPair, atomPairIndex);
			double radius2 = data.getAtomPairRadius2(*resPair, atomPairIndex);
			double lambda2 = data.getAtomPairLambda2(*resPair, atomPairIndex);
			double alpha2 = data.getAtomPairAlpha2(*resPair, atomPairIndex);
			
			double Xij = (r - radius1)/lambda1;
			double Xji = (r - radius2)/lambda2;
			resPairEnergy -= (alpha1*exp(-Xij*Xij) + alpha2*exp(-Xji*Xji))/r2;
		}
		
		// apply weight and offset
		if (atomPairIndex == 0) {
			resPairEnergy += resPair->offset;
		}
		energy += resPairEnergy*resPair->weight;
	}
	
	// compute the energy sum in SIMD-style
	// see url for a tutorial on GPU reductions:
	// http://developer.amd.com/resources/articles-whitepapers/opencl-optimization-case-study-simple-reductions/

	extern __shared__ double scratch[];
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
