package edu.duke.cs.osprey.partcr;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import edu.duke.cs.osprey.TestBase;
import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.astar.conf.order.AStarOrder;
import edu.duke.cs.osprey.astar.conf.order.StaticScoreHMeanAStarOrder;
import edu.duke.cs.osprey.astar.conf.scoring.AStarScorer;
import edu.duke.cs.osprey.astar.conf.scoring.MPLPPairwiseHScorer;
import edu.duke.cs.osprey.astar.conf.scoring.PairwiseGScorer;
import edu.duke.cs.osprey.astar.conf.scoring.mplp.NodeUpdater;
import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.confspace.PositionConfSpace;
import edu.duke.cs.osprey.confspace.RC;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.EnergyMatrixCalculator;
import edu.duke.cs.osprey.ematrix.SimpleEnergyCalculator;
import edu.duke.cs.osprey.ematrix.epic.EPICSettings;
import edu.duke.cs.osprey.energy.EnergyFunctionGenerator;
import edu.duke.cs.osprey.energy.MultiTermEnergyFunction;
import edu.duke.cs.osprey.pruning.PruningControl;
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
		
		// config
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
		
		// keep track of improvements over time
		List<Double> minBoundEnergies = new ArrayList<>();
		List<Double> minimizedEnergies = new ArrayList<>();
		
		int numPartCrIters = 20;
		for (int i=0; i<numPartCrIters; i++) {
			
			// do DEE
			int pruningInterval = 1; // don't set this to 0, pruning is too good, only 6 conformations in the tree
			boolean typeDep = false;
			double boundsThreshold = 100;
			int algoOption = 1;
			boolean useFlags = true;
			boolean useTriples = false;
			boolean useDacs = false;
			boolean useTupExp = false;
			double stericThreshold = 100;
			search.pruneMat = new PruningMatrix(search.confSpace, pruningInterval);
			//
			new PruningControl(
				search, search.pruneMat.getPruningInterval(), typeDep, boundsThreshold,
				algoOption, useFlags, useTriples, useDacs, useEpic, useTupExp, stericThreshold
			).prune();
			//
			
			// get the min bound conformation
			System.out.println("Finding min bound GMEC among " + search.confSpace.getNumConformations().doubleValue() + " conformations...");
			RCs rcs = new RCs(search.pruneMat);
			AStarOrder order = new StaticScoreHMeanAStarOrder();
			AStarScorer hscorer = new MPLPPairwiseHScorer(new NodeUpdater(), search.emat, 1, 0.0001);
			ConfAStarTree tree = new ConfAStarTree(order, new PairwiseGScorer(search.emat), hscorer, rcs);
			tree.initProgress();
			int[] minBoundConf = tree.nextConf();
			
			// get the min bound and minimized energies
			System.out.println("Minimizing...");
			double minBoundEnergy = search.emat.getInternalEnergy(new RCTuple(minBoundConf));
			double minimizedEnergy = search.minimizedEnergy(minBoundConf);
			double boundError = minimizedEnergy - minBoundEnergy;
			System.out.println("min bound energy: " + minBoundEnergy);
			System.out.println("minimized energy: " + minimizedEnergy);
			System.out.println("Bound error:      " + boundError);
			
			// on the first iteration, enumerate lower bounds to check against PartCR
			if (i == 0) {
				System.out.println("lower bounds");
				System.out.println(minBoundEnergy);
				for (int j=1; j<numPartCrIters; j++) {
					int[] conf = tree.nextConf();
					double bound = search.emat.getInternalEnergy(new RCTuple(conf));
					System.out.println(bound);
				}
			}
			
			minBoundEnergies.add(minBoundEnergy);
			minimizedEnergies.add(minimizedEnergy);
			
			// on the last iteration, don't partition
			if (i == numPartCrIters - 1) {
				break;
			}
			
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
			
			// get the rc at that pos
			RC bestRc = search.confSpace.posFlex.get(bestPos).RCs.get(minBoundConf[bestPos]);
			//System.out.println(dumpRc(bestRc, "best rc"));
			
			// NEXTTIME: the new A* implementation is super fast
			// which means we can bombard it with tons of new RCs and it won't slow down much
			// which means we should really aggressive with RC partioning
			// since A* time is basically zero compared to energy calculation time
			// how about partition the top 20% of worst errors instead of just the top 1?
			
			// partition the rc and update the conf space
			List<RC> splitRCs = new NAryRCSplitter().split(bestRc);
			List<Integer> rcMap = search.confSpace.posFlex.get(bestPos).replaceRC(bestRc, splitRCs);
			
			// update energy matrix
			System.out.println("updating energy matrix...");
			search.emat = updateEnergyMatrix(search.confSpace, search.emat, search.shellResidues, bestPos, rcMap);
		}
		
		// final report
		System.out.println("final improvements over " + numPartCrIters + " iterations:");
		System.out.println("min bound energies,minimized energies");
		for (int i=0; i<minBoundEnergies.size(); i++) {
			System.out.print(minBoundEnergies.get(i));
			System.out.print(",");
			System.out.print(minimizedEnergies.get(i));
			System.out.println();
		}
	}
	
	private static EnergyMatrix updateEnergyMatrix(ConfSpace confSpace, EnergyMatrix oldEmat, ArrayList<Residue> shellResidues, int updatedPos, List<Integer> rcMap) {
		
		// copy as many energies as possible from the old matrix
		// to a new matrix that is sized for the new conf space (with the added rotamers)
		EnergyMatrix newEmat = new EnergyMatrix(confSpace, oldEmat.getPruningInterval());
		
		SimpleEnergyCalculator ecalc = new SimpleEnergyCalculator(
			EnvironmentVars.curEFcnGenerator,
			confSpace,
			shellResidues
		);
		
    	for (int pos1=0; pos1<confSpace.numPos; pos1++) {
			for (int rc1=0; rc1<newEmat.getNumConfAtPos(pos1); rc1++) {

				// get the old rc1, if any
				Integer oldRc1 = rc1;
				if (pos1 == updatedPos) {
					oldRc1 = rcMap.get(rc1);
				}
				
				// one-body
				if (oldRc1 == null) {
					newEmat.setOneBody(pos1, rc1, ecalc.calcSingle(pos1, rc1));
				} else {
					newEmat.setOneBody(pos1, rc1, oldEmat.getOneBody(pos1, oldRc1));
				}
				
				// pairwise
				for (int pos2=0; pos2<pos1; pos2++) {
					for (int rc2=0; rc2<newEmat.getNumConfAtPos(pos2); rc2++) {
						
						// get the old rc2, if any
						Integer oldRc2 = rc2;
						if (pos2 == updatedPos) {
							oldRc2 = rcMap.get(rc2);
						}
						
						if (oldRc1 == null || oldRc2 == null) {
							newEmat.setPairwise(pos1, rc1, pos2, rc2, ecalc.calcPair(pos1, rc1, pos2, rc2));
						} else {
							newEmat.setPairwise(pos1, rc1, pos2, rc2, oldEmat.getPairwise(pos1, oldRc1, pos2, oldRc2));
						}
					}
				}
			}	
    	}
    	
    	return newEmat;
	}
	
	private static String dumpRc(RC rc, String title) {
		StringBuilder buf = new StringBuilder();
		buf.append(title + ":");
		for (int i=0; i<rc.DOFs.size(); i++) {
			buf.append(String.format("\t%s: [%f,%f]",
				rc.DOFs.get(i).getClass().getSimpleName(),
				rc.DOFmin.get(i),
				rc.DOFmax.get(i)
			));
		}
		return buf.toString();
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
