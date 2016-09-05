package edu.duke.cs.osprey;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;

import edu.duke.cs.osprey.astar.conf.ConfAStarNode;
import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.astar.conf.order.AStarOrder;
import edu.duke.cs.osprey.astar.conf.order.StaticScoreHMeanAStarOrder;
import edu.duke.cs.osprey.astar.conf.scoring.AStarScorer;
import edu.duke.cs.osprey.astar.conf.scoring.MPLPPairwiseHScorer;
import edu.duke.cs.osprey.astar.conf.scoring.PairwiseGScorer;
import edu.duke.cs.osprey.astar.conf.scoring.mplp.NodeUpdater;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimpleEnergyCalculator;
import edu.duke.cs.osprey.ematrix.epic.EPICSettings;
import edu.duke.cs.osprey.energy.EnergyFunctionGenerator;
import edu.duke.cs.osprey.energy.MultiTermEnergyFunction;
import edu.duke.cs.osprey.pruning.PruningControl;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tools.ObjectIO;
import edu.duke.cs.osprey.tupexp.LUTESettings;

public class LowerBoundsPlayground extends TestBase {
	
	public static void main(String[] args)
	throws Exception {
		
		// config
		initDefaultEnvironment();
		
		MultiTermEnergyFunction.setNumThreads(4);
		
		// make the search problem
		//String aaNames = "ALA VAL LEU ILE PHE TYR TRP CYS MET SER THR LYS ARG HIE HID ASP GLU ASN GLN GLY";
		//String aaNames = "ALA VAL LEU ILE GLU ASN GLN GLY";
		String aaNames = "ALA VAL LEU ILE";
		//String aaNames = "ALA";
		//String flexRes = "38 39 40 41 42 43 44";
		String flexRes = "38 39 40 41";
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
			flexResList, allowedAAs, addWt, doMinimize, useEpic, new EPICSettings(), useTupleExpansion,  new LUTESettings(),
			new DEEPerSettings(), moveableStrands, freeBBZones, useEllipses, useERef, addResEntropy, addWtRots, null, false
		);
		
		EnergyFunctionGenerator egen = EnvironmentVars.curEFcnGenerator;
		//SimpleEnergyCalculator.ShellDistribution dist = SimpleEnergyCalculator.ShellDistribution.AllOnSingles;
		//SimpleEnergyCalculator.ShellDistribution dist = SimpleEnergyCalculator.ShellDistribution.HalfSinglesHalfPairs;
		SimpleEnergyCalculator.ShellDistribution dist = SimpleEnergyCalculator.ShellDistribution.Even;
		//SimpleEnergyCalculator.ShellDistribution dist = SimpleEnergyCalculator.ShellDistribution.AllOnPairs;
		SimpleEnergyCalculator ecalc = new SimpleEnergyCalculator(egen, search.confSpace, search.shellResidues, dist);
		
		File ematFile = new File(String.format("/tmp/emat.%s.lb.dat", dist.name()));
		if (!ematFile.exists()) {
			System.out.println("\nComputing energy matrix...");
			search.emat = ecalc.calcEnergyMatrix();
			ObjectIO.writeObject(search.emat, ematFile.getAbsolutePath());
		} else {
			System.out.println("\nReading energy matrix...");
			search.emat = (EnergyMatrix)ObjectIO.readObject(ematFile.getAbsolutePath(), true);
		}
		
		// do DEE
		double pruningInterval = 3;
		search.pruneMat = new PruningMatrix(search.confSpace, pruningInterval);
		boolean typeDep = false;
		double boundsThreshold = 100;
		int algoOption = 1;
		boolean useFlags = true;
		boolean useTriples = false;
		boolean useDacs = false;
		boolean useTupExp = false;
		double stericThreshold = 100;
		new PruningControl(
			search, search.pruneMat.getPruningInterval(), typeDep, boundsThreshold,
			algoOption, useFlags, useTriples, useDacs, useEpic, useTupExp, stericThreshold
		).prune();
		
		// get the min bound conformation
		System.out.println("Finding min bound GMEC among " + search.confSpace.getNumConformations().doubleValue() + " conformations...");
		RCs rcs = new RCs(search.pruneMat);
		AStarOrder order = new StaticScoreHMeanAStarOrder();
		AStarScorer hscorer = new MPLPPairwiseHScorer(new NodeUpdater(), search.emat, 4, 0.0001);
		ConfAStarTree tree = new ConfAStarTree(order, new PairwiseGScorer(search.emat), hscorer, rcs);
		tree.initProgress();
		ConfSearch.ScoredConf minBoundConf = tree.nextConf();
		System.out.println("conformation:     " + Arrays.toString(minBoundConf.getAssignments()));
		
		// get the min bound and minimized energies
		double minimizedEnergy = search.minimizedEnergy(minBoundConf.getAssignments());
		double boundError = minimizedEnergy - minBoundConf.getScore();
		System.out.println("min bound energy: " + minBoundConf.getScore());
		System.out.println("minimized energy: " + minimizedEnergy);
		System.out.println("Bound error:      " + boundError);
		
		// check how many bounds there are between the min bound and the min GMEC
		int numBounds = 0;
		while (true) {
			numBounds++;
			ConfAStarNode node = tree.nextLeafNode();
			if (node == null || node.getGScore() >= minimizedEnergy) {
				break;
			}
		}
		System.out.println("num min bounds:   " + numBounds);
		
		double minBoundCheck = 0;
		double minimizedCheck = 0;
		
		// check the minimized energy to make sure we didn't break anything
		// when we were messing around with the energy functions
		for (int pos1=0; pos1<search.emat.getNumPos(); pos1++) {
			int rc1 = minBoundConf.getAssignments()[pos1];
			
			// single term
			double singleBoundEnergy = search.emat.getOneBody(pos1, rc1);
			double singleMinimizedEnergy = ecalc.getSingleEfunc(pos1).getEnergy();
			
			//double singleErr = singleMinimizedEnergy - singleBoundEnergy;
			//System.out.println(String.format("%5d: %f, %f, %f", pos1, singleBoundEnergy, singleMinimizedEnergy, singleErr));
			
			minBoundCheck += singleBoundEnergy;
			minimizedCheck += singleMinimizedEnergy;
			
			// pairwise energies
			for (int pos2=0; pos2<pos1; pos2++) {
				int rc2 = minBoundConf.getAssignments()[pos2];
				
				double pairwiseBoundEnergy = search.emat.getPairwise(pos1, rc1, pos2, rc2);
				double pairwiseMinimizedEnergy = ecalc.getPairEfunc(pos1, pos2).getEnergy();
				
				//double pairwiseErr = pairwiseMinimizedEnergy - pairwiseBoundEnergy;
				//System.out.println(String.format("%2d,%2d: %f, %f, %f", pos1, pos2, pairwiseBoundEnergy, pairwiseMinimizedEnergy, pairwiseErr));
				
				minBoundCheck += pairwiseBoundEnergy;
				minimizedCheck += pairwiseMinimizedEnergy;
			}
		}
		
		final double Epsilon = 1e-12;
		assert (Math.abs(minBoundConf.getScore() - minBoundCheck) <= Epsilon);
		assert (Math.abs(minimizedEnergy - minimizedCheck) <= Epsilon);
		System.out.println("minimized energy checks passed, energies are correct");
	}
}
