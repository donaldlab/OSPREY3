package edu.duke.cs.osprey.astar;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import edu.duke.cs.osprey.TestBase;
import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.astar.conf.order.AStarOrder;
import edu.duke.cs.osprey.astar.conf.order.DynamicHMeanAStarOrder;
import edu.duke.cs.osprey.astar.conf.order.StaticScoreHMeanAStarOrder;
import edu.duke.cs.osprey.astar.conf.scoring.AStarScorer;
import edu.duke.cs.osprey.astar.conf.scoring.MPLPPairwiseHScorer;
import edu.duke.cs.osprey.astar.conf.scoring.PairwiseGScorer;
import edu.duke.cs.osprey.astar.conf.scoring.TraditionalPairwiseHScorer;
import edu.duke.cs.osprey.astar.conf.scoring.mplp.NodeUpdater;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.ConfSearch.EnergiedConf;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.confspace.SearchProblem.MatrixType;
import edu.duke.cs.osprey.control.ConfEnergyCalculator;
import edu.duke.cs.osprey.control.ConfSearchFactory;
import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.control.MinimizingEnergyCalculator;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.epic.EPICSettings;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.pruning.PruningControl;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tools.ObjectIO;
import edu.duke.cs.osprey.tools.Stopwatch;
import edu.duke.cs.osprey.tupexp.LUTESettings;

public class BenchmarkAStar extends TestBase {
	
	public static void main(String[] args)
	throws Exception {
		
		initDefaultEnvironment();
		
		// make a few positions with lots of options for a very wide A* tree
		ResidueFlexibility resFlex = new ResidueFlexibility();
		resFlex.addMutable("3", "VAL MET TYR CYS ALA THR HIE HID GLY SER GLN ARG LYS ASN GLU ASP");
		resFlex.addMutable("5", "PHE ILE TRP LEU VAL MET TYR CYS ALA THR HIE HID GLY SER GLN");
		resFlex.addMutable("6", "HID HIE TYR");
		resFlex.addMutable("7", "PHE ILE TRP LEU VAL MET TYR CYS ALA THR HIE HID GLY SER GLN");
		resFlex.addMutable("13", "THR");
		resFlex.addMutable("15", "PHE ILE TRP LEU VAL MET TYR CYS ALA THR HIE HID GLY SER GLN");
		resFlex.addMutable("17", "VAL MET TYR CYS ALA THR HIE HID GLY SER GLN ARG LYS ASN GLU ASP");
		
		// settings
		boolean doMinimize = true;
		boolean addWt = true;
		boolean useEpic = false;
		boolean useTupleExpansion = false;
		boolean useEllipses = false;
		boolean useERef = true; // is this the problem?
		boolean addResEntropy = true;
		boolean addWtRots = true;
		boolean typeDep = true;
		double boundsThresh = 100;
		int algOption = 1;
		boolean useFlags = true;
		boolean useTriples = false;
		boolean preDacs = false;
		boolean useTupExp = false;
		double stericThresh = 10;
		double eW = 1;
		boolean doIMinDEE = true;
		double I0 = 0;
		double pruningInterval = 0 + eW;
		boolean checkApproxE = true;
		boolean outputGMECStruct = false;
		boolean eFullConfOnly = false;
		String confFileName = "confs.txt";
		
		// make the search problem
		ArrayList<String[]> moveableStrands = new ArrayList<String[]>();
		ArrayList<String[]> freeBBZones = new ArrayList<String[]>();
		SearchProblem search = new SearchProblem(
			"test", "/home/jeff/donaldlab/osprey test cases/anna-memory/1gua_1.pdb", 
			resFlex.flexResList, resFlex.allowedAAs, addWt, doMinimize, useEpic, new EPICSettings(), useTupleExpansion, new LUTESettings(),
			new DEEPerSettings(), moveableStrands, freeBBZones, useEllipses, useERef, addResEntropy, addWtRots, null, false, new ArrayList<>()
		);
		search.numEmatThreads = 3;
		
		PruningControl pruner = new PruningControl(
			search, pruningInterval, typeDep, boundsThresh, algOption,
			useFlags, useTriples, preDacs, useEpic, useTupExp, stericThresh
		);
		
		ConfSearchFactory confSearchFactory = (EnergyMatrix emat, PruningMatrix pmat) -> {
			boolean useMPLP = false;
			RCs rcs = new RCs(pmat);
			AStarOrder order;
			AStarScorer hscorer;
			if (useMPLP) {
				hscorer = new MPLPPairwiseHScorer(new NodeUpdater(), emat, 1, 0.0001);
				order = new StaticScoreHMeanAStarOrder();
			} else {
				hscorer = new TraditionalPairwiseHScorer(emat, rcs);
				order = new DynamicHMeanAStarOrder();
			}
			ConfAStarTree tree = new ConfAStarTree(order, new PairwiseGScorer(emat), hscorer, rcs);
			tree.initProgress();
			return tree;
		};
		
		ForcefieldParams ffparams = EnvironmentVars.curEFcnGenerator.ffParams;
		ConfEnergyCalculator.Async ecalc = MinimizingEnergyCalculator.make(ffparams, search, Parallelism.makeCpu(1), false);
		
		/* TEMP
		// use GMEC finder
		GMECFinder finder = new GMECFinder();
		finder.init(
			search, pruner, confSearchFactory, ecalc, eW, doIMinDEE, I0, doMinimize, useTupExp, useEpic,
			checkApproxE, outputGMECStruct, eFullConfOnly, confFileName, stericThresh
		);
		finder.calcGMEC();
		*/
	
		// compute the energy matrix
		String ematPath = "/tmp/BenchmarkAStar.emat.dat";
		// TEMP
		//new java.io.File(ematPath).delete();
		search.emat = (EnergyMatrix)ObjectIO.readObject(ematPath, true);
		if (search.emat == null) {
			search.emat = (EnergyMatrix)search.calcMatrix(MatrixType.EMAT);
			ObjectIO.writeObject(search.emat, ematPath);
		}
		
		// make an identity pruning matrix
		search.pruneMat = new PruningMatrix(search.confSpace, 0);
		pruner.prune();
			
		// build A* tree
		ConfSearch tree = confSearchFactory.make(search.emat, search.pruneMat);
		
		System.out.println(String.format("searching A* tree (" + tree.getNumConformations().floatValue() + " confs)..."));
		Stopwatch stopwatch = new Stopwatch().start();
		ScoredConf conf = tree.nextConf();
		System.out.println("first leaf node: " + stopwatch.stop().getTime(2));
		
		EnergiedConf econf = ecalc.calcEnergy(conf);
		System.out.println("\nMIN SCORE CONFORMATION");
		System.out.println("\tRCs      " + Arrays.toString(conf.getAssignments()));
		System.out.println(String.format("\tEnergy   %.6f", econf.getEnergy()));
		System.out.println(String.format("\tScore    %.6f (gap: %.6f)", conf.getScore(), econf.getEnergy() - conf.getScore()));
		
		stopwatch.reset().start();
		List<ScoredConf> confs = tree.nextConfs(econf.getEnergy() + 1);
		System.out.println(String.format("energy window (%d confs) %s", confs.size(), stopwatch.stop().getTime(2)));
		
		/* TODO: check the result
		assertThat(conf.getAssignments(), is(new int[] { 0, 0, 3, 0, 22, 4, 1 }));
		assertThat(conf.getScore(), isAbsolutely(-93.1587662665223));
		*/
	}
}
