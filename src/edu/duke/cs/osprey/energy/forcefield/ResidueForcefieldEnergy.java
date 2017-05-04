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

public class ResidueForcefieldEnergy implements EnergyFunction.DecomposableByDof {
	
	private static final long serialVersionUID = -4768384219061898745L;
	
	public final ResPairCache resPairCache;
	public final ResidueInteractions inters;
	public final List<Residue> residues;
	
	public final ResPair[] resPairs;
	
	private boolean isBroken;
	
	private double coulombFactor;
	private double scaledCoulombFactor;
	
	public ResidueForcefieldEnergy(ResPairCache resPairCache, ResidueInteractions inters, Molecule mol) {
		this(resPairCache, inters, mol.residues);
	}
	
	public ResidueForcefieldEnergy(ResPairCache resPairCache, ResidueInteractions inters, List<Residue> residues) {
		
		this.resPairCache = resPairCache;
		this.inters = inters;
		this.residues = residues;
		
		// compute solvation info if needed
		SolvationForcefield.ResiduesInfo solvInfo = null;
		if (resPairCache.ffparams.solvationForcefield != null) {
			solvInfo = resPairCache.ffparams.solvationForcefield.makeInfo(resPairCache.ffparams, residues);
		}
		
		// map the residue numbers to residues
		resPairs = new ResPair[inters.size()];
		int index = 0;
		for (ResidueInteractions.Pair pair : inters) {
			resPairs[index++] = resPairCache.get(residues, pair, solvInfo);
		}
		
		// is this a broken conformation?
		isBroken = false;
		for (ResPair pair : resPairs) {
			if (pair.res1.confProblems.size() + pair.res2.confProblems.size() > 0) {
				isBroken = true;
				
				// we're done here, no need to analyze broken conformations
				return;
			}
		}
		
		// pre-compute some constants needed by getEnergy()
		coulombFactor = ForcefieldParams.coulombConstant/resPairCache.ffparams.dielectric;
		scaledCoulombFactor = coulombFactor*resPairCache.ffparams.forcefld.coulombScaling;
	}
	
	@Override
	public double getEnergy() {
		return getEnergy(resPairs);
	}
	
	private double getEnergy(ResPair[] resPairs) {
		
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
			energy += resPairEnergy*pair.weight;
			energy += pair.offset;
		}
		
		return energy;
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
}
