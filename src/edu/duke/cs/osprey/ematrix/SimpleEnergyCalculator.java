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

	public Result calcSingle(int pos, int rc) {
		
		// make the energy function
		Residue res = getResidue(pos);
		EnergyFunction efunc = efuncGen.intraAndShellEnergy(res, shellResidues);

		return calc(efunc, new RCTuple(pos, rc));
	}
	
	public Result calcPair(int pos1, int rc1, int pos2, int rc2) {
		
		// make the energy function
		Residue res1 = getResidue(pos1);
		Residue res2 = getResidue(pos2);
		EnergyFunction efunc = efuncGen.resPairEnergy(res1, res2);
		
		return calc(efunc, new RCTuple(pos1, rc1, pos2, rc2));
	}
	
	private Result calc(EnergyFunction efunc, RCTuple tuple) {
		
		// optimize the degrees of freedom using the min efunc, if needed
		MoleculeModifierAndScorer mof = new MoleculeModifierAndScorer(efunc, confSpace, tuple);
		double[] minDofValues = null;
		if (mof.getNumDOFs() > 0) {
			CCDMinimizer ccdMin = new CCDMinimizer(mof, true);
			DoubleMatrix1D minDofVec = ccdMin.minimize();
			minDofValues = minDofVec.toArray();
			mof.setDOFs(minDofVec);
		}

		// calculate the energy using the eval efunc
		return new Result(minDofValues, efunc.getEnergy());
	}
	
	private Residue getResidue(int pos) {
		return confSpace.posFlex.get(pos).res;
	}
}
