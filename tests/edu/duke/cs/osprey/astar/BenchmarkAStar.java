package edu.duke.cs.osprey.astar;

import java.util.ArrayList;
import java.util.Arrays;

import edu.duke.cs.osprey.TestBase;
import edu.duke.cs.osprey.astar.conf.ConfAStarNode;
import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.astar.conf.order.AStarOrder;
import edu.duke.cs.osprey.astar.conf.order.DynamicHMeanAStarOrder;
import edu.duke.cs.osprey.astar.conf.order.StaticScoreHMeanAStarOrder;
import edu.duke.cs.osprey.astar.conf.scoring.AStarScorer;
import edu.duke.cs.osprey.astar.conf.scoring.MPLPPairwiseHScorer;
import edu.duke.cs.osprey.astar.conf.scoring.PairwiseGScorer;
import edu.duke.cs.osprey.astar.conf.scoring.TraditionalPairwiseHScorer;
import edu.duke.cs.osprey.astar.conf.scoring.mplp.NodeUpdater;
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
import edu.duke.cs.osprey.ematrix.SimpleEnergyCalculator;
import edu.duke.cs.osprey.ematrix.TermECalculator;
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
		
		// turning this on apparently hurts the A* heuristic
		boolean useERef = true;
		
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
		String ematPath = String.format("/tmp/BenchmarkAStar.emat.%s.dat", useERef ? "eref" : "base");
		search.emat = (EnergyMatrix)ObjectIO.readObject(ematPath, true);
		if (search.emat == null) {
			search.emat = (EnergyMatrix)search.calcMatrix(MatrixType.EMAT);
			ObjectIO.writeObject(search.emat, ematPath);
		}
		
		// make an identity pruning matrix
		search.pruneMat = new PruningMatrix(search.confSpace, 0);
		//pruner.prune();
			
		// build A* tree
		ConfAStarTree tree = (ConfAStarTree)confSearchFactory.make(search.emat, search.pruneMat);
		
		System.out.println(String.format("searching A* tree (" + tree.getNumConformations().floatValue() + " confs)..."));
		Stopwatch stopwatch = new Stopwatch().start();
		ConfAStarNode node = tree.nextLeafNode();
		System.out.println("first leaf node: " + stopwatch.stop().getTime(2));
		
		ScoredConf conf = new ScoredConf(node.makeConf(tree.rcs.getNumPos()), node.getScore());
		EnergiedConf econf = ecalc.calcEnergy(conf);
		
		System.out.println("\nMIN SCORE CONFORMATION");
		System.out.println("\tRCs      " + Arrays.toString(econf.getAssignments()));
		System.out.println(String.format("\tEnergy   %.6f", econf.getEnergy()));
		System.out.println(String.format("\tScore    %.6f (gap: %.6f)", econf.getScore(), econf.getEnergy() - econf.getScore()));
		
		//ScoredConf aconf = new ScoredConf(new int[] { 103, 108, 13, 113, 0, 111, 122 }, 0);
		ScoredConf aconf = conf;
		System.out.println("\nANALYZE CONF: " + Arrays.toString(aconf.getAssignments()));
		analyzeAstar(search, tree, aconf);
		
		/* enumerate energy window
		stopwatch.reset().start();
		List<ScoredConf> confs = tree.nextConfs(econf.getEnergy() + 1);
		System.out.println(String.format("energy window (%d confs) %s", confs.size(), stopwatch.stop().getTime(2)));
		*/
		
		/* TODO: check the result
		assertThat(conf.getAssignments(), is(new int[] { 0, 0, 3, 0, 22, 4, 1 }));
		assertThat(conf.getScore(), isAbsolutely(-93.1587662665223));
		*/
	}
	
	private static void analyzeAstar(SearchProblem search, ConfAStarTree tree, ScoredConf conf) {
		
		int[] a = conf.getAssignments();
		int[][] sequence = {
			{},
			{ a[0] },
			{ a[0], a[1] },
			{ a[0], a[1], a[2] },
			{ a[0], a[1], a[2], a[3] },
			{ a[0], a[1], a[2], a[3], a[4] },
			{ a[0], a[1], a[2], a[3], a[4], a[5] },
			{ a[0], a[1], a[2], a[3], a[4], a[5], a[6] }
		};
		for (int i=0; i<sequence.length; i++) {
			analyzeConf(search, tree, conf, sequence[i]);
			breakdownHScore(search, tree, conf, sequence[i]);
		}
		breakdownGScore(search, tree, conf);
	}
	
	private static void breakdownGScore(SearchProblem search, ConfAStarTree tree, ScoredConf conf) {
		
		ForcefieldParams ffparams = EnvironmentVars.curEFcnGenerator.ffParams;
		SimpleEnergyCalculator ecalcShell = new SimpleEnergyCalculator.Cpu(ffparams, search.confSpace, search.shellResidues);
		SimpleEnergyCalculator ecalcNoShell = new SimpleEnergyCalculator.Cpu(ffparams, search.confSpace, new ArrayList<>());
		
		int numPos = tree.rcs.getNumPos();
		
		// break down single energies
		System.out.println("G singles breakdown");
		for (int pos1=0; pos1<numPos; pos1++) {
			int rc1 = conf.getAssignments()[pos1];
			
			double newIntraEnergy = ecalcNoShell.calcSingle(pos1, rc1).energy;
			double newIntraAndShellEnergy = ecalcShell.calcSingle(pos1, rc1).energy;
			double newShellEnergy = newIntraAndShellEnergy - newIntraEnergy;
			
			double refEnergy = 0;
			if (search.emat.geteRefMat() != null) {
				refEnergy = search.emat.geteRefMat().posERef(pos1, rc1);
			}
			double entropyEnergy = search.confSpace.getRCResEntropy(pos1, rc1);
			double singleEnergy = search.emat.getOneBody(pos1, rc1);
			
			System.out.println(String.format("%1d   %10.6f + %10.6f + %10.6f + %10.6f = %10.6f",
				pos1, newIntraEnergy, newShellEnergy, entropyEnergy, -refEnergy, singleEnergy
			));
		}
		
		// break down pair energies
		System.out.println("G breakdown");
		for (int pos1=0; pos1<numPos; pos1++) {
			int rc1 = conf.getAssignments()[pos1];
			for (int pos2=0; pos2<pos1; pos2++) {
				int rc2 = conf.getAssignments()[pos2];
				System.out.print(String.format(" %6.2f", search.emat.getPairwise(pos1, rc1, pos2, rc2)));
			}
			System.out.println(String.format(" %6.2f", search.emat.getOneBody(pos1, rc1)));
		}
	}
	
	private static ConfIndex makeConfIndex(ScoredConf conf, int ... definedRCs) {
		ConfIndex index = new ConfIndex(conf.getAssignments().length);
		index.setNumDefined(definedRCs.length);
		for (int i=0; i<index.getNumDefined(); i++) {
			index.getDefinedPos()[i] = i;
			index.getDefinedRCs()[i] = definedRCs[i];
		}
		index.setNumUndefined(conf.getAssignments().length - definedRCs.length);
		for (int i=0; i<index.getNumUndefined(); i++) {
			index.getUndefinedPos()[i] = definedRCs.length + i;
		}
		return index;
	}
	
	private static void analyzeConf(SearchProblem search, ConfAStarTree tree, ScoredConf conf, int ... definedRCs) {

		ConfIndex index = makeConfIndex(conf, definedRCs);
		double gscore = tree.gscorer.calc(index, tree.rcs);
		double hscore = tree.hscorer.calc(index, tree.rcs);
		
		double truehscore = calcTrueHScore(search.emat, index, conf);
		//double naivehscore = new NaiveTraditionalPairwiseHScorer(emat).calc(index, tree.rcs);
		
		System.out.println(String.format("%2d   f %12.6f   g %12.6f   h %12.6f   th %12.6f   dh %12.6f   %% %6.1f",
			definedRCs.length,
			gscore + hscore,
			gscore,
			hscore,
			truehscore,
			Math.abs(truehscore - hscore),
			100f*(hscore - truehscore)/truehscore
		));
	}
	
	private static double calcTrueHScore(EnergyMatrix emat, ConfIndex index, ScoredConf conf) {
	
		double hscore = 0;
		
		// get the score for each undefined position
		for (int i=0; i<index.getNumUndefined(); i++) {
			int pos1 = index.getUndefinedPos()[i];
			int rc1 = conf.getAssignments()[pos1];
			
			// assignment to pos1
			hscore += emat.getOneBody(pos1, rc1);

			// interactions with defined residues
			for (int j=0; j<index.getNumDefined(); j++) {
				int pos2 = index.getDefinedPos()[j];
				int rc2 = index.getDefinedRCs()[j];
				hscore += emat.getPairwise(pos1, rc1, pos2, rc2);
			}

			// interactions with undefined residues
			for (int j=0; j<index.getNumUndefined(); j++) {
				int pos2 = index.getUndefinedPos()[j];
				if (pos2 >= pos1) {
					break;
				}

				// assignment to pos2
				int rc2 = conf.getAssignments()[pos2];
				hscore += emat.getPairwise(pos1, rc1, pos2, rc2);
			}
		}
		
		return hscore;
	}
	
	private static void breakdownHScore(SearchProblem search, ConfAStarTree tree, ScoredConf conf, int ... definedRCs) {
		
		ConfIndex index = makeConfIndex(conf, definedRCs);
		int numPos = tree.rcs.getNumPos();
		
		System.out.println("H breakdown");
		for (int pos1=0; pos1<numPos; pos1++) {
			
			RCPicks rcPicks = null;
			if (index.isUndefined(pos1)) {
				rcPicks = pickRCs(search, tree, index, pos1);
			}
			
			// pairs
			double pairEnergy = 0;
			for (int pos2=0; pos2<pos1; pos2++) {
				
				if (index.isDefined(pos1) && index.isDefined(pos2)) {
					
					// both defined
					int rc1 = conf.getAssignments()[pos1];
					int rc2 = conf.getAssignments()[pos2];
					pairEnergy = search.emat.getPairwise(pos1, rc1, pos2, rc2);
					
				} else if (index.isUndefined(pos1)) {
					
					// one defined
					pairEnergy = rcPicks.energies[pos2];
					
				} else {
					
					// both undefined
					assert (false);
				}
				
				System.out.print(String.format(" %6.2f", pairEnergy));
			}
			
			// single
			double singleEnergy = 0;
			if (index.isDefined(pos1)) {
				
				int rc1 = conf.getAssignments()[pos1];
				singleEnergy = search.emat.getOneBody(pos1, rc1);
				
			} else {
				
				singleEnergy = rcPicks.energies[pos1];
			}
			
			System.out.println(String.format(" %6.2f", singleEnergy));
		}
	}
	
	private static class RCPicks {
		
		public int[] rcs;
		public double[] energies;
		public double energy;
		
		public RCPicks(int numPos) {
			rcs = new int[numPos];
			Arrays.fill(rcs, -1);
			energies = new double[numPos];
			Arrays.fill(energies, 0);
			energy = 0;
		}
	}
	
	private static RCPicks pickRCs(SearchProblem search, ConfAStarTree tree, ConfIndex index, int pos1) {
		
		RCPicks minRCPicks = null;
		
		for (int rc1 : tree.rcs.get(pos1)) {
			
			RCPicks rcPicks = new RCPicks(tree.rcs.getNumPos());
			
			rcPicks.rcs[pos1] = rc1;
			rcPicks.energies[pos1] = search.emat.getOneBody(pos1, rc1);

			// interactions with defined residues
			for (int j=0; j<index.getNumDefined(); j++) {
				int pos2 = index.getDefinedPos()[j];
				int rc2 = index.getDefinedRCs()[j];
				rcPicks.rcs[pos2] = rc2;
				rcPicks.energies[pos2] = search.emat.getPairwise(pos1, rc1, pos2, rc2);
			}

			// interactions with undefined residues
			for (int j=0; j<index.getNumUndefined(); j++) {
				int pos2 = index.getUndefinedPos()[j];
				if (pos2 >= pos1) {
					break;
				}

				int rc2 = pickMinRC2(search, tree, pos1, rc1, pos2);
				rcPicks.rcs[pos2] = rc2;
				rcPicks.energies[pos2] = search.emat.getPairwise(pos1, rc1, pos2, rc2);
			}
			
			rcPicks.energy = 0;
			for (int i=0; i<tree.rcs.getNumPos(); i++) {
				rcPicks.energy += rcPicks.energies[i];
			}

			if (minRCPicks == null || rcPicks.energy < minRCPicks.energy) {
				minRCPicks = rcPicks;
			}
		}
		
		return minRCPicks;
	}
	
	private static int pickMinRC2(SearchProblem search, ConfAStarTree tree, int pos1, int rc1, int pos2) {
		double minEnergy = Double.POSITIVE_INFINITY;
		int minRC = -1;
		for (int rc2 : tree.rcs.get(pos2)) {
			double pairwiseEnergy = search.emat.getPairwise(pos1, rc1, pos2, rc2);
			if (pairwiseEnergy < minEnergy) {
				minEnergy = pairwiseEnergy;
				minRC = rc2;
			}
		}
		return minRC;
	}
}
