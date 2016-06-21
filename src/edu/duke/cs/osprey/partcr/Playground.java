package edu.duke.cs.osprey.partcr;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.TreeMap;

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
import edu.duke.cs.osprey.tools.ObjectIO;

public class Playground extends TestBase {
	
	public static void main(String[] args)
	throws Exception {
		
		// NEXTTIME:
		// bounds aren't improving fast enough after partitioning, need to figure out why
		// either the splitting is wrong somehow
		// or we're picking the wrong RCs to split
		// or the splits aren't being incorporated into the search problem correctly
		
		// config
		initDefaultEnvironment();
		
		MultiTermEnergyFunction.setNumThreads(4);
		
		// make the search problem
		//String aaNames = "ALA VAL LEU ILE PHE TYR TRP CYS MET SER THR LYS ARG HIE HID ASP GLU ASN GLN GLY";
		String aaNames = "ALA VAL LEU ILE GLU ASN GLN GLY";
		String flexRes = "38 39 40 41 42 43 44";
		ArrayList<String> flexResList = new ArrayList<>(Arrays.asList(flexRes.split(" ")));
		ArrayList<ArrayList<String>> allowedAAs = new ArrayList<>();
		for (int i=0; i<flexResList.size(); i++) {
			allowedAAs.add(new ArrayList<>(Arrays.asList(aaNames.split(" "))));
		}
		boolean doMinimize = true;
		boolean addWt = true;
		boolean useEpic = false;
		boolean useTupleExpansion = false;
		boolean useEllipses = false;
		boolean useERef = false;
		boolean addResEntropy = false;
		boolean addWtRots = false;
		ArrayList<String[]> moveableStrands = new ArrayList<String[]>();
		ArrayList<String[]> freeBBZones = new ArrayList<String[]>();
		SearchProblem search = new SearchProblem(
			"test", "test/1CC8/1CC8.ss.pdb", 
			flexResList, allowedAAs, addWt, doMinimize, useEpic, new EPICSettings(), useTupleExpansion,
			new DEEPerSettings(), moveableStrands, freeBBZones, useEllipses, useERef, addResEntropy, addWtRots
		);
		
		// compute the energy matrix
		File ematFile = new File(String.format("/tmp/emat.partcr.dat"));
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
			AStarScorer hscorer = new MPLPPairwiseHScorer(new NodeUpdater(), search.emat, 4, 0.0001);
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
			
			// on the first and last iterations, enumerate lower bounds to check against PartCR
			if (i == 0 || i == numPartCrIters - 1) {
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
			System.out.println("PartCR iteration " + (i+1));

			// sort all the positions by errors
			TreeMap<Double,Integer> positionsByError = new TreeMap<>();
			for (int pos=0; pos<search.confSpace.numPos; pos++) {
				int rot = search.confSpace.posFlex.get(pos).RCs.get(minBoundConf[pos]).rotNum;
				double posBoundEnergy = getPosBoundEnergy(search.emat, minBoundConf, pos);
				double posMinEnergy = getPosMinEnergy(search, pos);
				double posErr = posMinEnergy - posBoundEnergy;
				System.out.println(String.format("pos=%2d, rot=%3d, bound=%10f, energy=%10f, err=%10f",
					pos, rot, posBoundEnergy, posMinEnergy, posErr
				));
				positionsByError.put(posErr, pos);
			}
			List<Integer> positionsInOrder = new ArrayList<>(positionsByError.values());
			Collections.reverse(positionsInOrder);
			
			// the new A* implementation is super fast
			// which means we can bombard it with tons of new RCs and it won't slow down much
			// which means we should really aggressive with RC partitioning
			// since A* time is basically zero compared to energy calculation time
			// how about partition the top P% of worst errors instead of just the top 1?
			int numPosToSplit = search.confSpace.numPos/10; // 10%
			if (numPosToSplit <= 0) {
				numPosToSplit = 1;
			}
			
			// prep the mappings from new rcs to old rcs, to speed up energy matrix calculation
			List<List<Integer>> rcMaps = new ArrayList<>();
			for (int j=0; j<search.confSpace.numPos; j++) {
				rcMaps.add(null);
			}
			
			// split the positions with the worst bounds
			RCSplitter splitter = new NAryRCSplitter();
			for (int j=0; j<numPosToSplit; j++) {
			
				// get the rc at this pos
				int pos = positionsInOrder.get(j);
				PositionConfSpace posConfSpace = search.confSpace.posFlex.get(pos);
				RC rc = posConfSpace.RCs.get(minBoundConf[pos]);
				
				// TEMP
				System.out.println(String.format("splitting pos=%d, rot=%d", pos, rc.rotNum));
				
				// partition the rc and update the conf space
				List<RC> splitRCs = splitter.split(rc);
				List<Integer> rcMap = posConfSpace.replaceRC(rc, splitRCs);
				rcMaps.set(pos, rcMap);
			}
			
			// update energy matrix
			System.out.println("updating energy matrix...");
			search.emat = updateEnergyMatrix(search.confSpace, search.emat, search.shellResidues, rcMaps);
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
	
	private static EnergyMatrix updateEnergyMatrix(ConfSpace confSpace, EnergyMatrix oldEmat, ArrayList<Residue> shellResidues, List<List<Integer>> rcMaps) {
		
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
				if (rcMaps.get(pos1) != null) {
					oldRc1 = rcMaps.get(pos1).get(rc1);
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
						if (rcMaps.get(pos2) != null) {
							oldRc2 = rcMaps.get(pos2).get(rc2);
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
