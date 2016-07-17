
#ifdef cl_khr_fp64
    #pragma OPENCL EXTENSION cl_khr_fp64 : enable
#elif defined(cl_amd_fp64)
    #pragma OPENCL EXTENSION cl_amd_fp64 : enable
#else
    #error "Double precision floating point not supported by OpenCL implementation."
#endif


int getAtomIndex(int flags) {
	return abs(flags) - 1;
}

bool isHydrogen(int flags) {
	return flags > 0;
}

double calcSquaredRadius(global const double *coords, const int atom1Index, const int atom2Index) {
	int atom1Index3 = atom1Index*3;
	int atom2Index3 = atom2Index*3;
	double rx = coords[atom1Index3] - coords[atom2Index3];
	double ry = coords[atom1Index3 + 1] - coords[atom2Index3 + 1];
	double rz = coords[atom1Index3 + 2] - coords[atom2Index3 + 2];
	return rx*rx + ry*ry + rz*rz;
}

double calcElectrostatics(const double coulombFactor, const double r, const double charge) {
	return coulombFactor*charge/r;
}

bool useDistDepDielec(const int flags) {
	return (flags & 0x1) == 1;
}

bool useHEs(const int flags) {
	return ((flags >> 1) & 0x1) == 1;
}

bool useHVdw(const int flags) {
	return ((flags >> 2) & 0x1) == 1;
}

kernel void calc(
	global const double *coords, global const int *atomFlags, global const double *precomputed, global double *out,
	const int numPairs, const int num14Pairs, const double coulombFactor, const double scaledCoulombFactor,
	const double solvCutoff2, const double solvScale, const int flags
) {

	// NOTE: looks like we don't have enough gpu registers to inline everything
	// so need to call out to subroutines
	// and try to keep temp variables in small scopes
	// that's why the code looks all weird
	// it's optimized to minimize register usage at the expense of re-computing values sometimes

	// NOTE: can't have bools in kernel args for some reason, so the int flags encode 0=false, 1=true
	// pos 0 is useDistDepDielec
	// pos 1 is useHEs
	// pos 2 is useHVdw
	// packing the bools into one int saves registers too
	
	// start with zero energy
	out[get_global_id(0)] = 0;
	
	// which atom pair are we calculating?
	if (get_global_id(0) >= numPairs) {
		return;
	}
	
	// read atom flags and calculate all the things that use the atom flags in this scope
	bool bothHeavy;
	double r2;
	{
		int atom1Flags, atom2Flags;
		{
			int i2 = get_global_id(0)*2;
			atom1Flags = atomFlags[i2];
			atom2Flags = atomFlags[i2 + 1];
		}
		
		bothHeavy = !isHydrogen(atom1Flags) && !isHydrogen(atom2Flags);
		r2 = calcSquaredRadius(coords, getAtomIndex(atom1Flags), getAtomIndex(atom2Flags));
	}
	
	// calculate electrostatics
	if (bothHeavy || useHEs(flags)) {
	
		bool is14Pair = get_global_id(0) < num14Pairs;
		int i9 = get_global_id(0)*9;
		double charge = precomputed[i9 + 2];
		
		out[get_global_id(0)] += calcElectrostatics(
			is14Pair ? scaledCoulombFactor : coulombFactor,
			useDistDepDielec(flags) ? r2 : sqrt(r2),
			charge
		);
	}
	
	// calculate vdw
	if (bothHeavy || useHVdw(flags)) {
		
		double Aij, Bij;
		{
			int i9 = get_global_id(0)*9;
			Aij = precomputed[i9];
			Bij = precomputed[i9 + 1];
		}
		
		// compute vdw
		double r6 = r2*r2*r2;
		double r12 = r6*r6;
		out[get_global_id(0)] += Aij/r12 - Bij/r6;
	}
	
	// calculate solvation
	if (bothHeavy && r2 < solvCutoff2) {
			
		double r = sqrt(r2);
		{
			int i9 = get_global_id(0)*9;
			double lambda1 = precomputed[i9 + 3];
			double radius1 = precomputed[i9 + 4];
			double alpha1 = precomputed[i9 + 5];
			double Xij = (r - radius1)/lambda1;
			out[get_global_id(0)] -= alpha1*exp(-Xij*Xij)/r2*solvScale;
		}
		{
			int i9 = get_global_id(0)*9;
			double lambda2 = precomputed[i9 + 6];
			double radius2 = precomputed[i9 + 7];
			double alpha2 = precomputed[i9 + 8];
			double Xji = (r - radius2)/lambda2;
			out[get_global_id(0)] -= alpha2*exp(-Xji*Xji)/r2*solvScale;
		}
	}
}
