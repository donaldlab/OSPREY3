
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


// TEMP
#include <stdio.h>


typedef unsigned char byte;


// use same settings as CCDMinimizer on the java side
// TODO: use these
//const int MaxIterations = 30;
//const double ConvergenceThreshold = 0.001;
//const double Tolerance = 1e-6;
const double OneDegree = M_PI/180.0;
//const double InitialStep = OneDegree*0.25;
const double SolvCutoff = 9.0;
const double SolvCutoff2 = SolvCutoff*SolvCutoff;


typedef struct __align__(8) {
	unsigned int flags; // useHEs, useHvdW, distDepDielect, useEEF1
	unsigned int numDihedrals;
	unsigned int numResPairs;
	unsigned int maxNumAtoms;
	double coulombFactor;
	double scaledCoulombFactor;
} Header;

typedef struct __align__(8) {
	unsigned int resIndex;
	unsigned int numAtoms;
	unsigned int angleAtomOffsets[4];
	unsigned int numRotatedAtoms;
	unsigned int numResPairs;
} Dihedral;

typedef struct __align__(8) {
	unsigned long numAtomPairs;
	unsigned long atomsOffset1;
	unsigned long atomsOffset2;
	double weight;
	double offset;
} ResPair;

typedef struct __align__(8) {
	unsigned long flags; // bit isHeavyPair, bit is14Bonded, 6 bits space, 3 bytes space, short offset1, short offset2
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

class Data {
public:

	const Header & header;
	
	__device__ Data(const byte * const data) :
		m_data(data),
		header(*(Header *)&m_data[0]),
		m_dihedralOffsets((unsigned long *)&m_data[sizeof(Header)]),
		m_resPairOffsets((unsigned long *)&m_data[sizeof(Header) + sizeof(long)*header.numDihedrals])
	{
		// nothing else to do
	}
	
	__device__ const Dihedral & getDihedral(int d) const {
		return *(Dihedral *)&m_data[m_dihedralOffsets[d]];
	}

	__device__ const ResPair & getResPair(int i) const {
		return *(ResPair *)&m_data[m_resPairOffsets[i]];
	}
	
	__device__ const AtomPair & getAtomPair(const ResPair & resPair, int i) const {
		const byte * const p = (byte *)&resPair;
		return *(AtomPair *)&p[sizeof(ResPair) + sizeof(AtomPair)*i];
	}

private:
	const byte * const m_data;
	const unsigned long * const m_dihedralOffsets;
	const unsigned long * const m_resPairOffsets;
};


__device__ double blockSum(double threadEnergy, double * const threadEnergies) {

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
	
	return threadEnergies[0];
}

__device__ double calcFullEnergy(const Data & data, const double * const coords, double * const threadEnergies) {

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
	
	int resPairIndex = 0;
	int numAtomPairs = 0;
	for (int n=0; true; n++) {
	
		int atomPairIndex = blockDim.x*n + threadIdx.x;
	
		// find our res pair and atom pair
		const ResPair *resPair = NULL;
		const AtomPair *atomPair = NULL;
		int atomPairIndexInResPair;
		for (; resPairIndex<data.header.numResPairs; resPairIndex++) {
		
			resPair = &data.getResPair(resPairIndex);
			
			// is our atom pair in this res pair?
			atomPairIndexInResPair = atomPairIndex - resPair->numAtomPairs;
			if (atomPairIndexInResPair < resPair->numAtomPairs) {
			
				// yup, found it!
				atomPair = &data.getAtomPair(*resPair, atomPairIndexInResPair);
				break;	
			}
			
			numAtomPairs += resPair->numAtomPairs;
		}
		
		// stop if we ran out of atom pairs
		if (atomPair == NULL) {
			break;
		}
		
		// unpack the flags
		bool isHeavyPair, is14Bonded;
		unsigned long offset1, offset2;
		{
			unsigned long flags = atomPair->flags;
			offset2 = (flags & 0xffff) + resPair->atomsOffset2;
			flags >>= 16;
			offset1 = (flags & 0xffff) + resPair->atomsOffset1;
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
			if (is14Bonded) {
				if (distDepDielect) {
					resPairEnergy += data.header.scaledCoulombFactor*atomPair->charge/r2;
				} else {
					resPairEnergy += data.header.scaledCoulombFactor*atomPair->charge/r;
				}
			} else {
				if (distDepDielect) {
					resPairEnergy += data.header.coulombFactor*atomPair->charge/r2;
				} else {
					resPairEnergy += data.header.coulombFactor*atomPair->charge/r;
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
	
	return blockSum(energy, threadEnergies);
}

__device__ void copyx(const double * const src, double * const dest, const int size) {
	for (int i = threadIdx.x; i < size; i += blockDim.x) {
		dest[i] = src[i];
	}
	__syncthreads();
}

extern "C" __global__ void ccd(
	const byte * const rawdata,
	const double * const coords,
	double * const out
) {

	// parse data
	Data data(rawdata);
	
	/* TEMP
	if (threadIdx.x == 0) {
		const Dihedral &dihedral = data.getDihedral(0);
		const ResPair &resPair = data.getResPair(0);
		printf("data:            %016p\n", rawdata);
		printf("header:          %016p\n", &data.header);
		printf("dihedral[0]:     %016p\n", &data.getDihedral(0));
		printf("resPair[0]:      %016p\n", &data.getResPair(0));
		printf("atomPair[0][0]:  %016p\n", &data.getAtomPair(data.getResPair(0), 0));
		printf("dihedral[0].numResPairs: %d\n", data.getDihedral(0).numResPairs);
		printf("resPair[0].numAtomPairs: %d\n", data.getResPair(0).numAtomPairs);
	}
	*/
	
	// partition shared memory
	extern __shared__ byte shared[];
	double * const threadEnergies = (double *)shared;
	double * const resCoords = threadEnergies + blockDim.x;
	double * const nextx = resCoords + data.header.maxNumAtoms*3;
	double * const firstSteps = nextx + data.header.numDihedrals;
	double * const lastSteps = firstSteps + data.header.numDihedrals;
	
	// partition out memory
	double &outfx = out[data.header.numDihedrals];
	double * const outx = out;
	
	// init the step sizes
	for (int d = threadIdx.x; d < data.header.numDihedrals; d += blockDim.x) {
		firstSteps[d] = OneDegree;
		lastSteps[d] = OneDegree;
	}
	__syncthreads();
	
	// get the initial energy
	double herefx = calcFullEnergy(data, coords, threadEnergies);
	
	// TEMP
	if (threadIdx.x == 0) {
		printf("herefx = %12.6f\n", herefx);
	}
	
	// TODO
	//double * const herex = outx;
	
	// TEMP
	outfx = 5.0;
	int i = 0;
	for (; i<data.header.numDihedrals; i++) {
		outx[i] = 4.2;
	}
}