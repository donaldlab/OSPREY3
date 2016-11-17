package edu.duke.cs.osprey.ematrix;

import java.util.ArrayList;
import java.util.List;

import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.confspace.ParameterizedMoleculeCopy;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.energy.EnergyFunctionGenerator;
import edu.duke.cs.osprey.minimization.Minimizer;
import edu.duke.cs.osprey.minimization.MoleculeModifierAndScorer;
import edu.duke.cs.osprey.minimization.SimpleCCDMinimizer;
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
	
	private EnergyFunctionGenerator efuncGen;
	private ConfSpace confSpace;
	private List<Residue> shellResidues;
	private ShellDistribution dist;
	
	public SimpleEnergyCalculator(EnergyFunctionGenerator efuncGen, ConfSpace confSpace, List<Residue> shellResidues) {
		this(efuncGen, confSpace, shellResidues, DefaultDist);
	}
	
	public SimpleEnergyCalculator(EnergyFunctionGenerator efuncGen, ConfSpace confSpace, List<Residue> shellResidues, ShellDistribution dist) {
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
	
	public EnergyFunction getSingleEfunc(int pos, ParameterizedMoleculeCopy pmol) {
		double singleWeight = dist.getSingleWeight(confSpace.numPos);
		return efuncGen.intraAndDistributedShellEnergy(getResidue(pos, pmol), getResidues(shellResidues, pmol), confSpace.numPos, singleWeight);
	}
	
	public Minimizer.Result calcSingle(int pos, int rc) {
		return calcSingle(pos, rc, null);
	}
	
	public Minimizer.Result calcSingle(int pos, int rc, ParameterizedMoleculeCopy pmol) {
		EnergyFunction efunc = getSingleEfunc(pos, pmol);
		Minimizer.Result result = calc(efunc, new RCTuple(pos, rc), pmol);
		cleanup(efunc);
		return result;
	}
	
	public EnergyFunction getPairEfunc(int pos1, int pos2) {
		return getPairEfunc(pos1, pos2, null);
	}
	
	public EnergyFunction getPairEfunc(int pos1, int pos2, ParameterizedMoleculeCopy pmol) {
		double singleWeight = dist.getSingleWeight(confSpace.numPos);
		return efuncGen.resPairAndDistributedShellEnergy(getResidue(pos1, pmol), getResidue(pos2, pmol), getResidues(shellResidues, pmol), confSpace.numPos, singleWeight);
	}
	
	public Minimizer.Result calcPair(int pos1, int rc1, int pos2, int rc2) {
		return calcPair(pos1, rc1, pos2, rc2, null);
	}
	
	public Minimizer.Result calcPair(int pos1, int rc1, int pos2, int rc2, ParameterizedMoleculeCopy pmol) {
		EnergyFunction efunc = getPairEfunc(pos1, pos2, pmol);
		Minimizer.Result result = calc(efunc, new RCTuple(pos1, rc1, pos2, rc2), pmol);
		cleanup(efunc);
		return result;
	}
	
	public Minimizer.Result calc(EnergyFunction efunc, RCTuple tuple, ParameterizedMoleculeCopy pmol) {
		
		// put molecule in correct conformation for rcs
		MoleculeModifierAndScorer mof = new MoleculeModifierAndScorer(efunc, confSpace, tuple, pmol);
		
		// optimize the degrees of freedom, if needed
		if (mof.getNumDOFs() > 0) {
			// TODO: use configurable minimizer
			SimpleCCDMinimizer minimizer = new SimpleCCDMinimizer();
			minimizer.init(mof);
			return minimizer.minimize();
		}

		// otherwise, just evaulate the energy function
		return new Minimizer.Result(null, efunc.getEnergy());
	}
	
	private void cleanup(EnergyFunction efunc) {
		
		// cleanup the energy function if needed
		if (efunc instanceof EnergyFunction.NeedsCleanup) {
			((EnergyFunction.NeedsCleanup)efunc).cleanup();
		}
	}
	
	private Residue getResidue(int pos, ParameterizedMoleculeCopy pmol) {
		
		Residue res = confSpace.posFlex.get(pos).res;
		
		// if we're using a separate molecule, match this residue to the one in the molecule
		if (pmol != null) {
			res = pmol.getCopiedMolecule().residues.get(res.indexInMolecule);
		}
		
		return res;
	}
	
	private List<Residue> getResidues(List<Residue> residues, ParameterizedMoleculeCopy pmol) {
		
		// no molecule? just return the original residues
		if (pmol == null) {
			return residues;
		}
		
		// but if there is a molecule, then match the residues to it
		List<Residue> matched = new ArrayList<>(residues.size());
		for (Residue res : residues) {
			matched.add(pmol.getCopiedMolecule().residues.get(res.indexInMolecule));
		}
		return matched;
	}
}
