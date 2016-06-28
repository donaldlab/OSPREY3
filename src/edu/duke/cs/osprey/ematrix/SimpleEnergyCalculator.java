package edu.duke.cs.osprey.ematrix;

import java.util.ArrayList;

import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.energy.EnergyFunctionGenerator;
import edu.duke.cs.osprey.minimization.CCDMinimizer;
import edu.duke.cs.osprey.minimization.MoleculeModifierAndScorer;
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
	
	public SimpleEnergyCalculator(EnergyFunctionGenerator efuncGen, ConfSpace confSpace, ArrayList<Residue> shellResidues) {
		this.efuncGen = efuncGen;
		this.confSpace = confSpace;
		this.shellResidues = shellResidues;
	}
	
	public EnergyMatrix calcEnergyMatrix(ShellDistribution dist) {
		
		EnergyMatrix emat = new EnergyMatrix(confSpace, 0);
		
		for (int pos1=0; pos1<emat.getNumPos(); pos1++) {
			
			// singles
			System.out.println(String.format("calculating single energy %d...", pos1));
			for (int rc1=0; rc1<emat.getNumConfAtPos(pos1); rc1++) {
				emat.setOneBody(pos1, rc1, calcSingle(pos1, rc1, dist).getEnergy());
			}
				
			// pairwise
			for (int pos2=0; pos2<pos1; pos2++) {
				System.out.println(String.format("calculating pair energy %d,%d...", pos1, pos2));
				for (int rc1=0; rc1<emat.getNumConfAtPos(pos1); rc1++) {
					for (int rc2=0; rc2<emat.getNumConfAtPos(pos2); rc2++) {
						emat.setPairwise(pos1, rc1, pos2, rc2, calcPair(pos1, rc1, pos2, rc2, dist).getEnergy());
					}
				}
			}
		}
		
		return emat;
	}
	
	public EnergyFunction getSingleEfunc(int pos, int rc, ShellDistribution dist) {
		double singleWeight = dist.getSingleWeight(confSpace.numPos);
		if (singleWeight == 1.0) {
			return efuncGen.intraAndShellEnergy(getResidue(pos), shellResidues);
		} else {
			return efuncGen.intraAndDistributedShellEnergy(getResidue(pos), shellResidues, confSpace.numPos, singleWeight); 
		}
	}
	
	public Result calcSingle(int pos, int rc, ShellDistribution dist) {
		return calc(getSingleEfunc(pos, rc, dist), new RCTuple(pos, rc));
	}

	public EnergyFunction getPairEfunc(int pos1, int rc1, int pos2, int rc2, ShellDistribution dist) {
		double singleWeight = dist.getSingleWeight(confSpace.numPos);
		if (singleWeight == 1.0) {
			return efuncGen.resPairEnergy(getResidue(pos1), getResidue(pos2));
		} else {
			return efuncGen.resPairAndDistributedShellEnergy(getResidue(pos1), getResidue(pos2), shellResidues, confSpace.numPos, singleWeight); 
		}
	}
	
	public Result calcPair(int pos1, int rc1, int pos2, int rc2, ShellDistribution dist) {
		return calc(getPairEfunc(pos1, rc1, pos2, rc2, dist), new RCTuple(pos1, rc1, pos2, rc2));
	}
	
	private Result calc(EnergyFunction efunc, RCTuple tuple) {
		
		double[] minDofValues = null;
		
		// optimize the degrees of freedom, if needed
		MoleculeModifierAndScorer mof = new MoleculeModifierAndScorer(efunc, confSpace, tuple);
		if (mof.getNumDOFs() > 0) {
			CCDMinimizer ccdMin = new CCDMinimizer(mof, true);
			DoubleMatrix1D minDofVec = ccdMin.minimize();
			minDofValues = minDofVec.toArray();
			mof.setDOFs(minDofVec);
		}

		// calculate the energy
		return new Result(minDofValues, efunc.getEnergy());
	}
	
	private Residue getResidue(int pos) {
		return confSpace.posFlex.get(pos).res;
	}
}
