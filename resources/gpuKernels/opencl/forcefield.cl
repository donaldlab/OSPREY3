/*
** This file is part of OSPREY 3.0
** 
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
** 
** OSPREY relies on grants for its development, and since visibility
** in the scientific literature is essential for our success, we
** ask that users of OSPREY cite our papers. See the CITING_OSPREY
** document in this distribution for more information.
** 
** Contact Info:
**    Bruce Donald
**    Duke University
**    Department of Computer Science
**    Levine Science Research Center (LSRC)
**    Durham
**    NC 27708-0129
**    USA
**    e-mail: www.cs.duke.edu/brd/
** 
** <signature of Bruce Donald>, Mar 1, 2018
** Bruce Donald, Professor of Computer Science
*/

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

typedef struct __attribute__((aligned(8))) {
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
} ForcefieldArgs;
// sizeof = 40

kernel void calc(
	global const double *coords, global const int *atomFlags, global const double *precomputed, global const int *subsetTable, global double *out,
	global const ForcefieldArgs *args,
	local double *scratch
) {

	// NOTE: looks like we're running severely short on gpu registers
	// if the compiler says we used 29 registers, everything seems to work fine
	// if we use too many registers though, we get CL_OUT_OF_RESOURCES errors
	// not sure how many we're allowed to use though
	// so try to contain intermediate calculations in the smallest possible scopes
	// that's why the code looks all weird
	// it's optimized to minimize register usage at the expense of readability
	// we also re-compute values sometimes to save registers
	
	// NOTE: CL_OUT_OF_RESOURCES can also be thrown when you misuse a pointer (ie like a segfault)

	// start with zero energy
	double energy = 0;
	
	// which atom pair are we calculating?
	if (get_global_id(0) < args->numPairs) {
	
		int i = get_global_id(0);
		
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
				bool is14Pair = get_global_id(0) < args->num14Pairs;
				esEnergy *= is14Pair ? args->scaledCoulombFactor : args->coulombFactor;
			}
			
			{
				double charge = precomputed[i9 + 2];
				esEnergy *= charge;
			}
			
			{
				//esEnergy /= useDistDepDielec(flags) ? r2 : sqrt(r2);
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

	int locali = get_local_id(0);
	scratch[locali] = energy;
	
	barrier(CLK_LOCAL_MEM_FENCE);
	
	int localsize = get_local_size(0);
	for (int offset = 1; offset < localsize; offset <<= 1) {
	
		// sum this level of the reduction tree
		int mask = (offset << 1) - 1;
		if ((locali & mask) == 0) {
			scratch[locali] += scratch[locali + offset];
		}
		
		barrier(CLK_LOCAL_MEM_FENCE);
	}
	
	// finally, if we're the 0 thread, write the summed energy for this work group
	if (locali == 0) {
		out[get_group_id(0)] = scratch[0];
	}
}
