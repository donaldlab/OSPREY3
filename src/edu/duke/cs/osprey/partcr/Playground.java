package edu.duke.cs.osprey.partcr;

import java.io.File;
import java.util.ArrayList;

import edu.duke.cs.osprey.TestBase;
import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.astar.conf.order.AStarOrder;
import edu.duke.cs.osprey.astar.conf.order.StaticScoreHMeanAStarOrder;
import edu.duke.cs.osprey.astar.conf.scoring.AStarScorer;
import edu.duke.cs.osprey.astar.conf.scoring.MPLPPairwiseHScorer;
import edu.duke.cs.osprey.astar.conf.scoring.PairwiseGScorer;
import edu.duke.cs.osprey.astar.conf.scoring.mplp.NodeUpdater;
import edu.duke.cs.osprey.confspace.PositionConfSpace;
import edu.duke.cs.osprey.confspace.RC;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.EnergyMatrixCalculator;
import edu.duke.cs.osprey.ematrix.epic.EPICSettings;
import edu.duke.cs.osprey.energy.EnergyFunctionGenerator;
import edu.duke.cs.osprey.energy.MultiTermEnergyFunction;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.HashCalculator;
import edu.duke.cs.osprey.tools.ObjectIO;

public class Playground extends TestBase {
	
	public static class EmatConfig {
		
		String pdbPath;
		int firstRes;
		int numRes;
		boolean doMinimize;
		
		@Override
		public boolean equals(Object other) {
			if (other instanceof EmatConfig) {
				return equals((EmatConfig)other);
			}
			return false;
		}
		
		public boolean equals(EmatConfig other) {
			return this.pdbPath.equals(other.pdbPath)
				&& this.firstRes == other.firstRes
				&& this.numRes == other.numRes
				&& this.doMinimize == other.doMinimize;
		}
		
		@Override
		public int hashCode() {
			return HashCalculator.combineHashes(
				this.pdbPath.hashCode(),
				Integer.hashCode(this.firstRes),
				Integer.hashCode(this.numRes),
				Boolean.hashCode(this.doMinimize)
			);
		}
	}
	
	public static void main(String[] args)
	throws Exception {
		
		initDefaultEnvironment();
		
		MultiTermEnergyFunction.setNumThreads(4);
		
		EmatConfig ematConfig = new EmatConfig();
		ematConfig.pdbPath = "test/DAGK/2KDC.P.forOsprey.pdb";
		ematConfig.firstRes = 80;
		ematConfig.numRes = 20;
		ematConfig.doMinimize = true;
		
		// make the search problem
		ArrayList<String> flexRes = new ArrayList<>();
		ArrayList<ArrayList<String>> allowedAAs = new ArrayList<>();
		for (int i=0; i<ematConfig.numRes; i++) {
			flexRes.add(Integer.toString(ematConfig.firstRes + i + 1));
			allowedAAs.add(new ArrayList<String>());
		}
		boolean addWt = true;
		boolean useEpic = false;
		boolean useTupleExpansion = false;
		boolean useEllipses = false;
		boolean useERef = false;
		boolean addResEntropy = false;
		boolean addWtRots = true;
		ArrayList<String[]> moveableStrands = new ArrayList<String[]>();
		ArrayList<String[]> freeBBZones = new ArrayList<String[]>();
		SearchProblem search = new SearchProblem(
			"test", ematConfig.pdbPath, 
			flexRes, allowedAAs, addWt, ematConfig.doMinimize, useEpic, new EPICSettings(), useTupleExpansion,
			new DEEPerSettings(), moveableStrands, freeBBZones, useEllipses, useERef, addResEntropy, addWtRots
		);
		
		// compute the energy matrix
		File ematFile = new File(String.format("/tmp/emat.%d.dat", ematConfig.hashCode()));
		if (ematFile.exists()) {
			System.out.println("\nReading energy matrix...");
			search.emat = (EnergyMatrix)ObjectIO.readObject(ematFile.getAbsolutePath(), true);
		}
		if (search.emat == null) {
			System.out.println("\nComputing energy matrix...");
			EnergyMatrixCalculator emCalc = new EnergyMatrixCalculator(search.confSpace, search.shellResidues, useERef, addResEntropy);
			emCalc.calcPEM();
			search.emat = emCalc.getEMatrix();
			ObjectIO.writeObject(search.emat, ematFile.getAbsolutePath());
		}
		
		// calculate an "identity" pruning matrix (ie, no pruning)
		search.pruneMat = new PruningMatrix(search.confSpace, search.emat.getPruningInterval());
		
		// get the min bound conformation
		RCs rcs = new RCs(search.pruneMat);
		AStarOrder order = new StaticScoreHMeanAStarOrder();
		AStarScorer hscorer = new MPLPPairwiseHScorer(new NodeUpdater(), search.emat, 1, 0.0001);
		ConfAStarTree tree = new ConfAStarTree(order, new PairwiseGScorer(search.emat), hscorer, rcs);
		tree.initProgress();
		int[] minBoundConf = tree.nextConf();
		
		// get the min bound energy
		double minBoundEnergy = search.emat.getInternalEnergy(new RCTuple(minBoundConf));
		System.out.println("min bound energy: " + minBoundEnergy);
		
		// get the minimized energy and bound error
		System.out.println("Minimizing...");
		double minimizedEnergy = search.minimizedEnergy(minBoundConf);
		double boundError = minimizedEnergy - minBoundEnergy;
		System.out.println("minimized energy: " + minimizedEnergy);
		System.out.println("Bound error:      " + boundError);
		
		// find the pos with the biggest error
		int bestPos = -1;
		double bestErr = 0;
		for (int pos=0; pos<search.confSpace.numPos; pos++) {
			int res = ematConfig.firstRes + pos;
			double posBoundEnergy = getPosBoundEnergy(search.emat, minBoundConf, pos);
			double posMinEnergy = getPosMinEnergy(search, pos);
			double posErr = posMinEnergy - posBoundEnergy;
			System.out.println(String.format("%2d: res=%3d, bound=%10f, energy=%10f, err=%10f",
				pos, res, posBoundEnergy, posMinEnergy, posErr
			));
			
			if (posErr > bestErr) {
				bestErr = posErr;
				bestPos = pos;
			}
		}
		System.out.println("Pos with largest error: " + bestPos + ", err=" + bestErr);
		
		// TODO: partition that rc!
		RC bestRc = search.confSpace.posFlex.get(bestPos).RCs.get(minBoundConf[bestPos]);
	}
	
	private static double getPosBoundEnergy(EnergyMatrix emat, int[] conf, int pos1) {
		
		int rc1 = conf[pos1];
		double energy = emat.getOneBody(pos1, rc1);
		
		for (int pos2=0; pos2<emat.getNumPos(); pos2++) {
			if (pos1 != pos2) {
				int rc2 = conf[pos2];
				energy += emat.getPairwise(pos1, rc1, pos2, rc2)/2;
			}
		}
		
		return energy;
	}
	
	private static double getPosMinEnergy(SearchProblem search, int pos1) {
		
		EnergyFunctionGenerator egen = EnvironmentVars.curEFcnGenerator;
		double energy = 0;
		
		// single energy
		Residue res1 = search.confSpace.posFlex.get(pos1).res;
		energy += egen.singleResEnergy(res1).getEnergy();
		
		// single to shell energy
		for (Residue shellResidue : search.shellResidues) {
			energy += egen.resPairEnergy(res1, shellResidue).getEnergy();
		}
		
		// pairwise energy
		for (PositionConfSpace flexPos : search.confSpace.posFlex) {
			Residue res2 = flexPos.res;
			if (res1 != res2) {
				energy += egen.resPairEnergy(res1, res2).getEnergy()/2;
			}
		}
		
		return energy;
	}
}
