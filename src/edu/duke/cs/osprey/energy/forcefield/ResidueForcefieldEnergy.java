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

package edu.duke.cs.osprey.energy.forcefield;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.energy.ResidueInteractions;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams.SolvationForcefield;
import edu.duke.cs.osprey.energy.forcefield.ResPairCache.ResPair;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.structure.Residues;

public class ResidueForcefieldEnergy implements EnergyFunction.DecomposableByDof {
	
	private static final long serialVersionUID = -4768384219061898745L;
	
	public final ResPairCache resPairCache;
	public final ResidueInteractions inters;
	public final Residues residues;
	
	public final ResPair[] resPairs;
	public final boolean isBroken;
	
	private double coulombFactor;
	private double scaledCoulombFactor;
	
	public ResidueForcefieldEnergy(ResPairCache resPairCache, ResidueInteractions inters, Molecule mol) {
		this(resPairCache, inters, mol.residues);
	}
	
	public ResidueForcefieldEnergy(ResPairCache resPairCache, ResidueInteractions inters, Residues residues) {
		
		this.resPairCache = resPairCache;
		this.inters = inters;
		// NOTE: don't filter residues using the interactions here... for some reason this is breaking EPIC
		//this.residues = inters.filter(residues);
		this.residues = residues;
		
		// compute solvation info if needed
		SolvationForcefield.ResiduesInfo solvInfo = null;
		if (resPairCache.ffparams.solvationForcefield != null) {
			solvInfo = resPairCache.ffparams.solvationForcefield.makeInfo(resPairCache.ffparams, this.residues);
		}
		
		// map the residue numbers to residues
		resPairs = new ResPair[inters.size()];
		int index = 0;
		for (ResidueInteractions.Pair pair : inters) {
			resPairs[index++] = resPairCache.get(this.residues, pair, solvInfo);
		}
		
		// is this a broken conformation?
		for (ResPair pair : resPairs) {
			if (pair.res1.confProblems.size() + pair.res2.confProblems.size() > 0) {
				isBroken = true;
				
				// we're done here, no need to analyze broken conformations
				return;
			}
		}
		isBroken = false;

		// pre-compute some constants needed by getEnergy()
		coulombFactor = ForcefieldParams.coulombConstant/resPairCache.ffparams.dielectric;
		scaledCoulombFactor = coulombFactor*resPairCache.ffparams.forcefld.coulombScaling;
	}

	public ResidueForcefieldEnergy makeSubset(ResidueInteractions.Pair pair) {
		return makeSubset(new ResidueInteractions(pair));
	}

	public ResidueForcefieldEnergy makeSubset(ResidueInteractions inters) {
		return new ResidueForcefieldEnergy(resPairCache, inters, residues);
	}
	
	@Override
	public double getEnergy() {
		return getEnergy(resPairs);
	}
	
	private double getEnergy(ResPair[] resPairs) {

		// NOTE: this function gets hammered a lot! Performance is super important here,
		// and even pedantic optimizations can make a big difference.
		// don't make changes here unless you're carefully profiling too
		
		// check broken-ness first. easy peasy
		if (isBroken) {
			return Double.POSITIVE_INFINITY;
		}
		
		// copy stuff to the stack/registers, to improve CPU cache performance
		boolean useHEs = resPairCache.ffparams.hElect;
		boolean useHvdW = resPairCache.ffparams.hVDW;
		double coulombFactor = this.coulombFactor;
		double scaledCoulombFactor = this.scaledCoulombFactor;
		boolean distDepDielect = resPairCache.ffparams.distDepDielect;
		boolean useEEF1 = resPairCache.ffparams.solvationForcefield == SolvationForcefield.EEF1;
		
		double energy = 0;
		
		for (int i=0; i<resPairs.length; i++) {
			ResPair pair = resPairs[i];
			
			// copy pair values/references to the stack/registers
			// so we don't have to touch the pair object while inside the atom pair loop
			double[] coords1 = pair.res1.coords;
			double[] coords2 = pair.res2.coords;
			int numAtomPairs = pair.info.numAtomPairs;
			long[] flags = pair.info.flags;
			double[] precomputed = pair.info.precomputed;
			
			double resPairEnergy = 0;
			
			// for each atom pair...
			int pos = 0;
			for (int j=0; j<numAtomPairs; j++) {
				
				// read the flags
				// NOTE: this is efficient, but destructive to the val
				long atomPairFlags = flags[j];
				int atomOffset2 = (int)(atomPairFlags & 0xffff);
				atomPairFlags >>= 16;
				int atomOffset1 = (int)(atomPairFlags & 0xffff);
				atomPairFlags >>= 46;
				boolean isHeavyPair = (atomPairFlags & 0x1) == 0x1;
				atomPairFlags >>= 1;
				boolean is14Bonded = (atomPairFlags & 0x1) == 0x1;
				
				// get the radius
				double r2;
				double r;
				{
					// read atom coords
					double x1 = coords1[atomOffset1];
					double y1 = coords1[atomOffset1 + 1];
					double z1 = coords1[atomOffset1 + 2];
					double x2 = coords2[atomOffset2];
					double y2 = coords2[atomOffset2 + 1];
					double z2 = coords2[atomOffset2 + 2];
					
					double d;
					
					d = x1 - x2;
					r2 = d*d;
					d = y1 - y2;
					r2 += d*d;
					d = z1 - z2;
					r2 += d*d;
					r = Math.sqrt(r2);
				}
				
				// read the bit flags
				
				// electrostatics
				if (isHeavyPair || useHEs) {
					double charge = precomputed[pos++];
					if (is14Bonded) {
						if (distDepDielect) {
							resPairEnergy += scaledCoulombFactor*charge/r2;
						} else {
							resPairEnergy += scaledCoulombFactor*charge/r;
						}
					} else {
						if (distDepDielect) {
							resPairEnergy += coulombFactor*charge/r2;
						} else {
							resPairEnergy += coulombFactor*charge/r;
						}
					}
				} else {
					pos++;
				}
				
				// van der Waals
				if (isHeavyPair || useHvdW) {
					
					double Aij = precomputed[pos++];
					double Bij = precomputed[pos++];
					
					// compute vdw
					double r6 = r2*r2*r2;
					double r12 = r6*r6;
					resPairEnergy += Aij/r12 - Bij/r6;
					
				} else {
					pos += 2;
				}
				
				// solvation
				if (useEEF1) {
					if (isHeavyPair && r2 < ForcefieldParams.solvCutoff2) {
							
						double radius1 = precomputed[pos++];
						double lambda1 = precomputed[pos++];
						double alpha1 = precomputed[pos++];
						double radius2 = precomputed[pos++];
						double lambda2 = precomputed[pos++];
						double alpha2 = precomputed[pos++];
						
						// compute solvation energy
						double Xij = (r - radius1)/lambda1;
						double Xji = (r - radius2)/lambda2;
						resPairEnergy -= (alpha1*Math.exp(-Xij*Xij) + alpha2*Math.exp(-Xji*Xji))/r2;
						
					} else {
						pos += 6;
					}
				}
			}
			
			// apply weights and offsets
			energy += (resPairEnergy + pair.offset + pair.solvEnergy)*pair.weight;
		}
		
		return energy;
	}

	public double getElectrostaticsEnergy() {
		return getElectrostaticsEnergy(resPairs);
	}

	public double getElectrostaticsEnergy(ResPair[] resPairs) {

		// NOTE: this function isn't hammered like getEnergy() is, so performance isn't important here

		// check broken-ness first. easy peasy
		if (isBroken) {
			return Double.POSITIVE_INFINITY;
		}

		double energy = 0.0;

		for (ResPair pair : resPairs) {

			double resPairEnergy = 0;

			// for each atom pair...
			int pos = 0;
			for (int j=0; j<pair.info.numAtomPairs; j++) {

				// read the flags
				long atomPairFlags = pair.info.flags[j];
				int atomOffset2 = (int)(atomPairFlags & 0xffff);
				atomPairFlags >>= 16;
				int atomOffset1 = (int)(atomPairFlags & 0xffff);
				atomPairFlags >>= 46;
				boolean isHeavyPair = (atomPairFlags & 0x1) == 0x1;
				atomPairFlags >>= 1;
				boolean is14Bonded = (atomPairFlags & 0x1) == 0x1;

				// get the radius
				double r2 = getR2(pair, atomOffset1, atomOffset2);
				double r = Math.sqrt(r2);

				// compute the electrostatics
				if (isHeavyPair || resPairCache.ffparams.hElect) {
					double charge = pair.info.precomputed[pos++];
					if (is14Bonded) {
						if (resPairCache.ffparams.distDepDielect) {
							resPairEnergy += scaledCoulombFactor*charge/r2;
						} else {
							resPairEnergy += scaledCoulombFactor*charge/r;
						}
					} else {
						if (resPairCache.ffparams.distDepDielect) {
							resPairEnergy += coulombFactor*charge/r2;
						} else {
							resPairEnergy += coulombFactor*charge/r;
						}
					}
				} else {
					pos++;
				}

				// skip the van der Waals
				pos += 2;

				// skip the solvation precomputed values
				if (resPairCache.ffparams.solvationForcefield == SolvationForcefield.EEF1) {
					pos += 6;
				}
			}

			// apply weight, but not offset
			energy += resPairEnergy*pair.weight;
		}

		return energy;
	}

	public double getVanDerWaalsEnergy() {
		return getVanDerWaalsEnergy(resPairs);
	}

	public double getVanDerWaalsEnergy(ResPair[] resPairs) {

		// NOTE: this function isn't hammered like getEnergy() is, so performance isn't important here

		// check broken-ness first. easy peasy
		if (isBroken) {
			return Double.POSITIVE_INFINITY;
		}

		double energy = 0.0;

		for (ResPair pair : resPairs) {

			double resPairEnergy = 0;

			// for each atom pair...
			int pos = 0;
			for (int j=0; j<pair.info.numAtomPairs; j++) {

				// read the flags
				long atomPairFlags = pair.info.flags[j];
				int atomOffset2 = (int)(atomPairFlags & 0xffff);
				atomPairFlags >>= 16;
				int atomOffset1 = (int)(atomPairFlags & 0xffff);
				atomPairFlags >>= 46;
				boolean isHeavyPair = (atomPairFlags & 0x1) == 0x1;

				// skip the electrostatics precomputed values
				pos += 1;

				if (isHeavyPair || resPairCache.ffparams.hVDW) {

					// get the radius
					double r2 = getR2(pair, atomOffset1, atomOffset2);

					double Aij = pair.info.precomputed[pos++];
					double Bij = pair.info.precomputed[pos++];

					// compute vdw
					double r6 = r2*r2*r2;
					double r12 = r6*r6;
					resPairEnergy += Aij/r12 - Bij/r6;

				} else {
					pos += 2;
				}

				// skip the solvation precomputed values
				if (resPairCache.ffparams.solvationForcefield == SolvationForcefield.EEF1) {
					pos += 6;
				}
			}

			// apply weight, but not offset
			energy += resPairEnergy*pair.weight;
		}

		return energy;
	}

	public double getSolvationEnergy() {
		return getSolvationEnergy(resPairs);
	}

	public double getSolvationEnergy(ResPair[] resPairs) {

		// NOTE: this function isn't hammered like getEnergy() is, so performance isn't important here

		// check broken-ness first. easy peasy
		if (isBroken) {
			return Double.POSITIVE_INFINITY;
		}

		double energy = 0.0;

		for (ResPair pair : resPairs) {

			double resPairEnergy = 0;

			// for each atom pair...
			int pos = 0;
			for (int j=0; j<pair.info.numAtomPairs; j++) {

				// read the flags
				// NOTE: this is efficient, but destructive to the val
				long atomPairFlags = pair.info.flags[j];
				int atomOffset2 = (int)(atomPairFlags & 0xffff);
				atomPairFlags >>= 16;
				int atomOffset1 = (int)(atomPairFlags & 0xffff);
				atomPairFlags >>= 46;
				boolean isHeavyPair = (atomPairFlags & 0x1) == 0x1;
				atomPairFlags >>= 1;
				boolean is14Bonded = (atomPairFlags & 0x1) == 0x1;

				// get the radius
				double r2 = getR2(pair, atomOffset1, atomOffset2);
				double r = Math.sqrt(r2);

				// skip the vdW and electrostatics precomputed values
				pos += 3;

				// solvation
				if (resPairCache.ffparams.solvationForcefield == SolvationForcefield.EEF1) {
					if (isHeavyPair && r2 < ForcefieldParams.solvCutoff2) {

						double radius1 = pair.info.precomputed[pos++];
						double lambda1 = pair.info.precomputed[pos++];
						double alpha1 = pair.info.precomputed[pos++];
						double radius2 = pair.info.precomputed[pos++];
						double lambda2 = pair.info.precomputed[pos++];
						double alpha2 = pair.info.precomputed[pos++];

						// compute solvation energy
						double Xij = (r - radius1)/lambda1;
						double Xji = (r - radius2)/lambda2;
						resPairEnergy -= (alpha1*Math.exp(-Xij*Xij) + alpha2*Math.exp(-Xji*Xji))/r2;

					} else {
						pos += 6;
					}
				}
			}

			// apply weight, but not offset
			energy += (resPairEnergy + pair.solvEnergy)*pair.weight;
		}

		return energy;
	}

	public double getOffsetsEnergy() {
		return getOffsetsEnergy(resPairs);
	}

	public double getOffsetsEnergy(ResPair[] resPairs) {

		// NOTE: this function isn't hammered like getEnergy() is, so performance isn't important here

		// check broken-ness first. easy peasy
		if (isBroken) {
			return Double.POSITIVE_INFINITY;
		}

		double energy = 0.0;

		for (ResPair pair : resPairs) {

			// add the regular offsets
			energy += pair.offset*pair.weight;
		}

		return energy;
	}

	private ResPair findResPair(ResidueInteractions.Pair interPair) {
		for (ResPair resPair : resPairs) {
			if (interPair.resNum1.equalsIgnoreCase(resPair.res1.getPDBResNumber())
				&& interPair.resNum2.equalsIgnoreCase(resPair.res2.getPDBResNumber())) {
				return resPair;
			}
		}
		return null;
	}

	private double getR2(ResPair pair, int atomOffset1, int atomOffset2) {

		// read atom coords
		double x1 = pair.res1.coords[atomOffset1];
		double y1 = pair.res1.coords[atomOffset1 + 1];
		double z1 = pair.res1.coords[atomOffset1 + 2];
		double x2 = pair.res2.coords[atomOffset2];
		double y2 = pair.res2.coords[atomOffset2 + 1];
		double z2 = pair.res2.coords[atomOffset2 + 2];

		// compute r2
		double dx = x1 - x2;
		double dy = y1 - y2;
		double dz = z1 - z2;
		return dx*dx + dy*dy + dz*dz;
	}
	
	public ResPair[] makeResPairsSubset(Residue res) {
	
		// pass 1: count
		int num = 0;
		for (ResPair resPair : resPairs) {
			if (resPair.res1 == res || resPair.res2 == res) {
				num++;
			}
		}
		
		// pass 2: collect
		ResPair[] pairs = new ResPair[num];
		num = 0;
		for (ResPair resPair : resPairs) {
			if (resPair.res1 == res || resPair.res2 == res) {
				pairs[num++] = resPair;
			}
		}
		return pairs;
	}
	
	public int[] makeResPairIndicesSubset(Residue res) {
		
		// pass 1: count
		int num = 0;
		for (ResPair resPair : resPairs) {
			if (resPair.res1 == res || resPair.res2 == res) {
				num++;
			}
		}
		
		// pass 2: collect
		int[] indices = new int[num];
		num = 0;
		for (int i=0; i<resPairs.length; i++) {
			if (resPairs[i].res1 == res || resPairs[i].res2 == res) {
				indices[num++] = i;
			}
		}
		return indices;
	}
	
	@Override
	public List<EnergyFunction> decomposeByDof(Molecule mol, List<DegreeOfFreedom> dofs) {
		
		class Subset implements EnergyFunction {
			
			private static final long serialVersionUID = 4664215035458391734L;
			
			private ResPair[] resPairs;
			
			@Override
			public double getEnergy() {
				return ResidueForcefieldEnergy.this.getEnergy(resPairs);
			}
		}
		
		Map<Residue,Subset> cache = new HashMap<>();
		
		List<EnergyFunction> efuncs = new ArrayList<>();
		for (DegreeOfFreedom dof : dofs) {
			Residue res = dof.getResidue();
			
			if (res == null) {
				
				// no res, just use the whole efunc
				efuncs.add(this);
				
			} else {
				
				// make a subset energy function
				Subset subset = cache.get(res);
				if (subset == null) {
					subset = new Subset();
					subset.resPairs = makeResPairsSubset(res);
					cache.put(res, subset);
				}
				efuncs.add(subset);
			}
		}
		
		return efuncs;
	}

	public static class Vdw extends ResidueForcefieldEnergy {

		public Vdw(ResPairCache resPairCache, ResidueInteractions inters, Molecule mol) {
			super(resPairCache, inters, mol);
		}

		public Vdw(ResPairCache resPairCache, ResidueInteractions inters, Residues residues) {
			super(resPairCache, inters, residues);
		}

		public Vdw(ResidueForcefieldEnergy efunc) {
			this(efunc.resPairCache, efunc.inters, efunc.residues);
		}

		@Override
		public double getEnergy() {
			return super.getVanDerWaalsEnergy();
		}

		@Override
		public List<EnergyFunction> decomposeByDof(Molecule mol, List<DegreeOfFreedom> dofs) {

			class Subset implements EnergyFunction {

				private static final long serialVersionUID = 4664215035458391734L;

				private ResPair[] resPairs;

				@Override
				public double getEnergy() {
					return Vdw.this.getVanDerWaalsEnergy(resPairs);
				}
			}

			Map<Residue,Subset> cache = new HashMap<>();

			List<EnergyFunction> efuncs = new ArrayList<>();
			for (DegreeOfFreedom dof : dofs) {
				Residue res = dof.getResidue();

				if (res == null) {

					// no res, just use the whole efunc
					efuncs.add(this);

				} else {

					// make a subset energy function
					Subset subset = cache.get(res);
					if (subset == null) {
						subset = new Subset();
						subset.resPairs = makeResPairsSubset(res);
						cache.put(res, subset);
					}
					efuncs.add(subset);
				}
			}

			return efuncs;
		}
	}
}
