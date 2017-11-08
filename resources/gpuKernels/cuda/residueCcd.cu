
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
	-maxrregcount 64
	"kernelSource/ccd.cu" -o "kernelBinaries/ccd.bin"
	
	See Maxwell compatibility guide for more info:
	http://docs.nvidia.com/cuda/maxwell-compatibility-guide/index.html#building-maxwell-compatible-apps-using-cuda-6-0
*/


typedef unsigned char byte;


// use same settings as CCDMinimizer on the java side
const int MaxIterations = 30;
const double ConvergenceThreshold = 0.001;
const double Tolerance = 1e-6;
const double OneDegree = M_PI/180.0;
const double InitialStep = OneDegree*0.25;
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
	unsigned long atomsOffset;
	unsigned short angleAtomOffsets[4];
	unsigned int numRotatedAtoms;
	unsigned int numResPairs;
	double xdmin;
	double xdmax;
} Dihedral;

typedef struct __align__(8) {
	unsigned long numAtomPairs;
	unsigned long atomsOffset1;
	unsigned long atomsOffset2;
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

__device__ unsigned int roundUpToMultiple(unsigned int val, unsigned int base) {
	int mod = val % base;
	if (mod == 0) {
		return val;
	}
	return val + base - mod;
}

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
	
	__device__ const Dihedral & getDihedral(const unsigned int d) const {
		return *(Dihedral *)&m_data[m_dihedralOffsets[d]];
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
	
	__device__ unsigned short getRotatedAtomOffset(const Dihedral & dihedral, const unsigned int i) const {
		const byte * const p = (byte *)&dihedral;
		const unsigned short * const atomOffsets = (unsigned short *)&(p[sizeof(Dihedral)]);
		return atomOffsets[i];
	}
	
	__device__ const ResPair & getResPair(const Dihedral & dihedral, const unsigned int i) const {
		const byte * const p = (byte *)&dihedral;
		const unsigned int * const resPairIndices = (unsigned int *)&(p[sizeof(Dihedral) + sizeof(short)*roundUpToMultiple(dihedral.numRotatedAtoms, 4)]);
		return getResPair(resPairIndices[i]);
	}

private:
	const byte * const m_data;
	const unsigned long * const m_dihedralOffsets;
	const unsigned long * const m_resPairOffsets;
	
	__device__ double getAtomPairDouble(const ResPair & resPair, const unsigned int i, const unsigned int pos) const {
		const byte * const p = (byte *)&resPair;
		const double * const a = (double *)&p[sizeof(ResPair) + sizeof(long)*resPair.numAtomPairs + sizeof(double)*pos*resPair.numAtomPairs];
		return a[i];
	}
};

typedef struct {
	double step;
	double xdstar;
	double fxdstar;
} LinesearchOut;


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


__device__ double calcEnergy(const double * const coords, const Data & data, const Dihedral * const dihedral, double * const threadEnergies) {

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
	
	unsigned int numResPairs;
	if (dihedral == NULL) {
		numResPairs = data.header.numResPairs;
	} else {
		numResPairs = dihedral->numResPairs;
	}
	
	unsigned int resPairIndex = 0;
	unsigned int numAtomPairs = 0;
	for (int n = threadIdx.x; true; n += blockDim.x) {
	
		// find our res pair and atom pair
		const ResPair *resPair;
		int atomPairIndex = -1;
		
		for (; resPairIndex<numResPairs; resPairIndex++) {
		
			if (dihedral == NULL) {
				resPair = &data.getResPair(resPairIndex);
			} else {
				resPair = &data.getResPair(*dihedral, resPairIndex);
			}
			
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
	
	return blockSum(energy, threadEnergies);
}

__device__ void copyx(const double * const src, double * const dest, const int size) {
	for (int i = threadIdx.x; i < size; i += blockDim.x) {
		dest[i] = src[i];
	}
	__syncthreads();
}


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

__device__ double3 readCoord(const double * const coords, const Dihedral & dihedral, const unsigned int i) {
	unsigned long offset = dihedral.atomsOffset + i;
	return make_double3(coords[offset], coords[offset + 1], coords[offset + 2]);
}

__device__ void writeCoord(double * const coords, const Dihedral & dihedral, unsigned int i, const double3 &val) {
	unsigned long offset = dihedral.atomsOffset + i;
	coords[offset] = val.x;
	coords[offset + 1] = val.y;
	coords[offset + 2] = val.z;
}

__device__ void pose(double * const coords, const Data & data, const Dihedral & dihedral, double dihedralRadians) {

	// get the four atom positions: a, b, c, d
	double3 a = readCoord(coords, dihedral, dihedral.angleAtomOffsets[0]);
	double3 b = readCoord(coords, dihedral, dihedral.angleAtomOffsets[1]);
	double3 c = readCoord(coords, dihedral, dihedral.angleAtomOffsets[2]);
	double3 d = readCoord(coords, dihedral, dihedral.angleAtomOffsets[3]);
	
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
	if (threadIdx.x < dihedral.numRotatedAtoms) {
		unsigned short index = data.getRotatedAtomOffset(dihedral, threadIdx.x);
		double3 p = readCoord(coords, dihedral, index);
		sub(p, b);
		rotateVecInverse(p, rx, ry, rz);
		rotateVecZ(p, deltaSin, deltaCos);
		rotateVec(p, rx, ry, rz);
		add(p, b);
		writeCoord(coords, dihedral, index, p);
	}
	__syncthreads();
}

__device__ double poseAndCalcDihedralEnergy(double * const coords, const Data & data, const Dihedral & dihedral, double * const threadEnergies, double dihedralRadians) {
	pose(coords, data, dihedral, dihedralRadians);
	return calcEnergy(coords, data, &dihedral, threadEnergies);
}

__device__ double getTolerance(double f) {
	
	// scale abs(f) by tolerance, unless f is very small
	return Tolerance * fmax(1.0, fabs(f));
}

__device__ LinesearchOut linesearch(double * const coords, const Data & data, const Dihedral & dihedral, double * const threadEnergies, const double xd, const double step) {

	// sample energy at the starting point
	double fxd = poseAndCalcDihedralEnergy(coords, data, dihedral, threadEnergies, xd);
	
	// get the positive (p) neighbor
	double fxdp = INFINITY;
	{
		double xdp = xd + step;
		if (xdp <= dihedral.xdmax) {
			fxdp = poseAndCalcDihedralEnergy(coords, data, dihedral, threadEnergies, xdp);
		}
	}
	
	// get the negative (n) neighbor
	double fxdm = INFINITY;
	{
		double xdm = xd - step;
		if (xdm  >= dihedral.xdmin) {
			fxdm = poseAndCalcDihedralEnergy(coords, data, dihedral, threadEnergies, xdm);
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
	if (xdstar < dihedral.xdmin) {
		xdstar = dihedral.xdmin;
	}
	if (xdstar > dihedral.xdmax) {
		xdstar = dihedral.xdmax;
	}
	
	double fxdstar = poseAndCalcDihedralEnergy(coords, data, dihedral, threadEnergies, xdstar);
	
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
			if (xdsurfNext < dihedral.xdmin) {
				
				// if the min is better, go there instead
				if (isnan(fxdmin)) {
					fxdmin = poseAndCalcDihedralEnergy(coords, data, dihedral, threadEnergies, dihedral.xdmin);
				}
				if (fxdmin < fxdsurfHere) {
					xdsurfHere = dihedral.xdmin;
					fxdsurfHere = fxdmin;
				}
				
				break;
			
			// did we step off the max?
			} else if (xdsurfNext > dihedral.xdmax) {
				
				// if the max is better, go there instead
				if (isnan(fxdmax)) {
					fxdmax = poseAndCalcDihedralEnergy(coords, data, dihedral, threadEnergies, dihedral.xdmax);
				}
				if (fxdmax < fxdsurfHere) {
					xdsurfHere = dihedral.xdmax;
					fxdsurfHere = fxdmax;
				}
				
				break;
			}
			
			double fxdsurfNext = poseAndCalcDihedralEnergy(coords, data, dihedral, threadEnergies, xdsurfNext);
			
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
			double fxdsurfNext = poseAndCalcDihedralEnergy(coords, data, dihedral, threadEnergies, xdsurfNext);
			
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
	if (xdm >= dihedral.xdmin) {
		fxdm = poseAndCalcDihedralEnergy(coords, data, dihedral, threadEnergies, xdm);
		if (fxdm < fxdstar) {
			xdstar = xdm;
			fxdstar = fxdm;
		}
	}
	if (xdp <= dihedral.xdmax) {
		fxdp = poseAndCalcDihedralEnergy(coords, data, dihedral, threadEnergies, xdp);
		if (fxdp < fxdstar) {
			xdstar = xdp;
			fxdstar = fxdp;
		}
	}
	
	// one last pose
	pose(coords, data, dihedral, xdstar);
	
	// set outputs
	out.xdstar = xdstar;
	out.fxdstar = fxdstar;
	return out;
}

extern "C" __global__ void ccd(
	const byte * const rawdata,
	double * const coords,
	double * const xin,
	double * const out
) {

	// parse data
	const Data data(rawdata);
	
	// partition shared memory
	extern __shared__ byte shared[];
	double * const threadEnergies = (double *)shared;
	double * const nextx = threadEnergies + blockDim.x;
	double * const firstSteps = nextx + data.header.numDihedrals;
	double * const lastSteps = firstSteps + data.header.numDihedrals;
	
	// partition out memory
	double & herefx = out[data.header.numDihedrals];
	double * const herex = out;
	
	// init the step sizes
	for (int d = threadIdx.x; d < data.header.numDihedrals; d += blockDim.x) {
		firstSteps[d] = OneDegree;
		lastSteps[d] = OneDegree;
	}
	__syncthreads();
	
	// get the initial energy
	herefx = calcEnergy(coords, data, NULL, threadEnergies);
	
	// init starting x
	for (int d = threadIdx.x; d < data.header.numDihedrals; d += blockDim.x) {
		herex[d] = xin[d];
	}
	__syncthreads();
	
	// for each iteraction of CCD...
	for (int iter=0; iter<MaxIterations; iter++) {

		copyx(herex, nextx, data.header.numDihedrals);
		
		// for each dimension/dihedral...
		for (int d=0; d<data.header.numDihedrals; d++) {
			const Dihedral & dihedral = data.getDihedral(d);
			
			double xd = nextx[d];
			
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
				while (dihedral.xdmax > dihedral.xdmin && xd - step < dihedral.xdmin && xd + step > dihedral.xdmax) {
					step /= 2;
				}
			}
			
			// do line search
			LinesearchOut lsout = linesearch(coords, data, dihedral, threadEnergies, xd, step);
			
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
		}
		
		// evaluate the whole energy function
		double nextfx = calcEnergy(coords, data, NULL, threadEnergies);
		double improvement = herefx - nextfx;
		
		if (improvement > 0) {
		
			// take the step
			copyx(nextx, herex, data.header.numDihedrals);
			herefx = nextfx;
			
			if (improvement < ConvergenceThreshold) {
				break;
			}
			
		} else {
			break;
		}
	}
}