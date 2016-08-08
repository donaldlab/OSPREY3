package edu.duke.cs.osprey.ematrix;

import java.util.ArrayList;
import java.util.List;

import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.energy.EnergyFunctionGenerator;
import edu.duke.cs.osprey.minimization.CCDMinimizer;
import edu.duke.cs.osprey.minimization.MoleculeModifierAndScorer;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.Residue;

public class SimpleEnergyCalculator {
	
	public static enum ShellDistribution {
		
		AllOnSingles {
			@Override
			public double getSingleWeight(int numPos) {
				return 1.0;
			}
		},
		HalfSinglesHalfPairs {
			@Override
			public double getSingleWeight(int numPos) {
				return 0.5;
			}
		},
		Even {
			@Override
			public double getSingleWeight(int numPos) {
				return 1.0/numPos;
			}
		},
		AllOnPairs {
			@Override
			public double getSingleWeight(int numPos) {
				return 0.0;
			}
		};
		
		public abstract double getSingleWeight(int numPos);
	}
	
	// NOTE: empirical testing didn't reveal any distribution to be consistently better than
	// the others on a variety of designs.
	// It probably depends on the ratio of flexible positions to shell positions
	// anyway, AllOnSingles is waaay faster since it doesn't compute shell energies for pairwise terms
	// also, the other energy calculator uses this distribution, so we should match to be consistent by default
	public static ShellDistribution DefaultDist = ShellDistribution.AllOnSingles;
	
	public static class Result {
		
		private double[] minDofValues;
		private double energy;
		
		public Result(double[] minDofValues, double energy) {
			this.minDofValues = minDofValues;
			this.energy = energy;
		}
		
		public double[] getDofValues() {
			return minDofValues;
		}
		
		public double getEnergy() {
			return energy;
		}
	}
	
	private EnergyFunctionGenerator efuncGen;
	private ConfSpace confSpace;
	private ArrayList<Residue> shellResidues;
	private ShellDistribution dist;
	
	public SimpleEnergyCalculator(EnergyFunctionGenerator efuncGen, ConfSpace confSpace, ArrayList<Residue> shellResidues) {
		this(efuncGen, confSpace, shellResidues, DefaultDist);
	}
	
	public SimpleEnergyCalculator(EnergyFunctionGenerator efuncGen, ConfSpace confSpace, ArrayList<Residue> shellResidues, ShellDistribution dist) {
		this.efuncGen = efuncGen;
		this.confSpace = confSpace;
		this.shellResidues = shellResidues;
		this.dist = dist;
	}
	
	public EnergyFunctionGenerator getEnergyFunctionGenerator() {
		return efuncGen;
	}
	
	public ShellDistribution getShellDistribution() {
		return dist;
	}
	
	public ConfSpace getConfSpace() {
		return confSpace;
	}
	
	public EnergyFunction getSingleEfunc(int pos) {
		return getSingleEfunc(pos, null);
	}
	
	public EnergyFunction getSingleEfunc(int pos, Molecule mol) {
		double singleWeight = dist.getSingleWeight(confSpace.numPos);
		return efuncGen.intraAndDistributedShellEnergy(getResidue(pos, mol), getResidues(shellResidues, mol), confSpace.numPos, singleWeight);
	}
	
	public Result calcSingle(int pos, int rc) {
		return calcSingle(pos, rc, null);
	}
	
	public Result calcSingle(int pos, int rc, Molecule mol) {
		return calc(getSingleEfunc(pos, mol), new RCTuple(pos, rc), mol);
	}
	
	public EnergyFunction getPairEfunc(int pos1, int pos2) {
		return getPairEfunc(pos1, pos2, null);
	}
	
	public EnergyFunction getPairEfunc(int pos1, int pos2, Molecule mol) {
		double singleWeight = dist.getSingleWeight(confSpace.numPos);
		return efuncGen.resPairAndDistributedShellEnergy(getResidue(pos1, mol), getResidue(pos2, mol), getResidues(shellResidues, mol), confSpace.numPos, singleWeight);
	}
	
	public Result calcPair(int pos1, int rc1, int pos2, int rc2) {
		return calcPair(pos1, rc1, pos2, rc2, null);
	}
	
	public Result calcPair(int pos1, int rc1, int pos2, int rc2, Molecule mol) {
		return calc(getPairEfunc(pos1, pos2, mol), new RCTuple(pos1, rc1, pos2, rc2), mol);
	}
	
	public Result calc(EnergyFunction efunc, RCTuple tuple, Molecule mol) {
		
		double[] minDofValues = null;
		
		// put molecule in correct conformation for rcs
		MoleculeModifierAndScorer mof = new MoleculeModifierAndScorer(efunc, confSpace, tuple, mol);
		
		// optimize the degrees of freedom, if needed
		if (mof.getNumDOFs() > 0) {
			CCDMinimizer ccdMin = new CCDMinimizer(mof, true);
			DoubleMatrix1D minDofVec = ccdMin.minimize();
			minDofValues = minDofVec.toArray();
			mof.setDOFs(minDofVec);
		}

		// calculate the energy
		return new Result(minDofValues, efunc.getEnergy());
	}
	
	private Residue getResidue(int pos, Molecule mol) {
		
		Residue res = confSpace.posFlex.get(pos).res;
		
		// if we're using a separate molecule, match this residue to the one in the molecule
		if (mol != null) {
			res = mol.getResByPDBResNumber(res.getPDBResNumber());
		}
		
		return res;
	}
	
	private List<Residue> getResidues(List<Residue> residues, Molecule mol) {
		
		// no molecule? just return the original residues
		if (mol == null) {
			return residues;
		}
		
		// but if there is a molecule, then match the residues to it
		List<Residue> matched = new ArrayList<>(residues.size());
		for (Residue res : residues) {
			matched.add(mol.getResByPDBResNumber(res.getPDBResNumber()));
		}
		return matched;
	}
}
