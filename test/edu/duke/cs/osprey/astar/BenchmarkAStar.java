/*
** This file is part of OSPREY 3.0
** 
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
** 
** OSPREY relies on grants for its development, and since visibility
** in the scientific literature is essential for our success, we
** ask that users of OSPREY cite our papers. See the CITING_OSPREY
** document in this distribution for more information.
** 
** Contact Info:
**    Bruce Donald
**    Duke University
**    Department of Computer Science
**    Levine Science Research Center (LSRC)
**    Durham
**    NC 27708-0129
**    USA
**    e-mail: www.cs.duke.edu/brd/
** 
** <signature of Bruce Donald>, Mar 1, 2018
** Bruce Donald, Professor of Computer Science
*/

package edu.duke.cs.osprey.astar;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

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
import edu.duke.cs.osprey.astar.conf.scoring.mplp.EdgeUpdater;
import edu.duke.cs.osprey.confspace.ConfSearch.EnergiedConf;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.confspace.SearchProblem.MatrixType;
import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimpleEnergyCalculator;
import edu.duke.cs.osprey.ematrix.epic.EPICSettings;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.gmec.GMECConfEnergyCalculator;
import edu.duke.cs.osprey.gmec.EnergyRange;
import edu.duke.cs.osprey.gmec.MinimizingConfEnergyCalculator;
import edu.duke.cs.osprey.parallelism.Parallelism;
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
		boolean useVoxelG = false;
		ArrayList<String[]> moveableStrands = new ArrayList<String[]>();
		ArrayList<String[]> freeBBZones = new ArrayList<String[]>();
		
		// turning this on apparently hurts the A* heuristic
		boolean useERef = true;
		
		boolean useMPLP = true;
		int astarThreads = 4;
		double Ew = 1;
		
		// make the search problem
		SearchProblem search = new SearchProblem(
			"test", "/home/jeff/donaldlab/osprey test cases/anna-memory/1gua_1.pdb", 
			resFlex.flexResList, resFlex.allowedAAs, addWt, doMinimize, useEpic, new EPICSettings(), useTupleExpansion, new LUTESettings(),
			new DEEPerSettings(), moveableStrands, freeBBZones, useEllipses, useERef, addResEntropy, addWtRots, null, useVoxelG, new ArrayList<>()
		);
		
		ForcefieldParams ffparams = EnvironmentVars.curEFcnGenerator.ffParams;
		GMECConfEnergyCalculator.Async ecalc = MinimizingConfEnergyCalculator.make(ffparams, search, Parallelism.makeCpu(1));
		
		// compute the energy matrix
		String ematPath = String.format("/tmp/BenchmarkAStar.emat.%s.dat", useERef ? "eref" : "base");
		search.emat = (EnergyMatrix)ObjectIO.readObject(ematPath, true);
		if (search.emat == null) {
			search.numEmatThreads = 3;
			search.emat = (EnergyMatrix)search.calcMatrix(MatrixType.EMAT);
			ObjectIO.writeObject(search.emat, ematPath);
		}
		
		// make an identity pruning matrix
		search.pruneMat = new PruningMatrix(search.confSpace, 0);
			
		// build A* tree
		RCs rcs = new RCs(search.pruneMat);
		AStarOrder order;
		AStarScorer hscorer;
		if (useMPLP) {
			//hscorer = new MPLPPairwiseHScorer(new NodeUpdater(), search.emat, 1, 0.0001);
			hscorer = new MPLPPairwiseHScorer(new EdgeUpdater(), search.emat, 5, 0.0001);
			order = new StaticScoreHMeanAStarOrder();
		} else {
			hscorer = new TraditionalPairwiseHScorer(search.emat, rcs);
			order = new DynamicHMeanAStarOrder();
		}
		ConfAStarTree tree = new ConfAStarTree.Builder(search.emat, rcs)
			.setCustom(order, new PairwiseGScorer(search.emat), hscorer)
			.build();
		tree.initProgress();
		tree.setParallelism(Parallelism.makeCpu(astarThreads));
		
		System.out.println(String.format("searching A* tree (" + tree.getNumConformations().floatValue() + " confs)..."));
		ScoredConf conf = tree.nextConf();
		EnergiedConf econf = ecalc.calcEnergy(conf);
		
		System.out.println("\nMIN SCORE CONFORMATION");
		System.out.println("\tRCs      " + Arrays.toString(econf.getAssignments()));
		System.out.println(String.format("\tEnergy   %.6f", econf.getEnergy()));
		System.out.println(String.format("\tScore    %.6f (gap: %.6f)", econf.getScore(), econf.getEnergy() - econf.getScore()));
		
		/* DEBUG: heuristic analysis
		//ScoredConf aconf = new ScoredConf(new int[] { 103, 108, 13, 113, 0, 111, 122 }, 0);
		ScoredConf aconf = conf;
		System.out.println("\nANALYZE CONF: " + Arrays.toString(aconf.getAssignments()));
		analyzeAstar(search, tree, aconf);
		*/
		
		// TREE SIZE
		// Traditional:        expanded: 69k     queued:  9.4m
		// MPLP 1 node iter:   expanded: ~4.8k   queued: ~0.5m
		// MPLP 1 edge iter:   expanded: ~1.9k   queued: ~0.2m
		// MPLP 2 edge iters:  expanded:   130   queued:  ~11k
		// MPLP 3 edge iters:  expanded:    60   queued:   ~5k
		// MPLP 4 edge iters:  expanded:    42   queued:   ~4k
		// MPLP 5 edge iters:  expanded:    30   queued:   ~3k
		// MPLP 6 edge iters:  expanded:    21   queued:   ~2k
		// MPLP 7 edge iters:  expanded:    19   queued: ~1.8k
		// MPLP 8 edge iters:  expanded:    16   queued: ~1.6k
		
		// BENCHMARKING (on laptop)
		// traditional:     ~1.2 min    100k-160k scores/sec
		// MPLP 1 thread:   ~12.1 min   0.4k-0.9k scores/sec
		// MPLP 3 threads:  ~12.1 min   1.0k-1.5k scores/sec
		
		// BENCHMARKING (on desktop)
		// traditional, 1 thread:          ~1.1 min     100k-160k scores/sec
		// traditional, 2 threads:         ~1.7 min      60k-125k scores/sec
		
		// (desktop, with youtube playing)
		// traditional, 1 thread:          ~1.2 min      120k-170k scores/sec
		// traditional, 2 threads:         ~1.9 min       70k-110k scores/sec
		// traditional, 4 threads:         ~2.2 min        70k-90k scores/sec
		
		// with 1k queue
		// MPLP 1 node iter, 1 thread:      ??? min       500-700 scores/sec (for first 1k expansions)
		// MPLP 1 node iter, 2 threads:     ??? min      700-1000 scores/sec (for first 1k expansions)
		// MPLP 1 node iter, 4 threads:    ~6.4 min      900-1200 scores/sec (for first 1k expansions)
		
		// with no queue
		// MPLP 1 node iter, 4 threads:    ~5.3 min      1300-1700 scores/sec (for first 1k expansions)
		// MPLP 1 edge iter, 4 threads:    ~1.1 min      2600-3000 scores/sec (for first 1k expansions)
		// MPLP 2 edge iters, 4 threads:	 33 sec            ??? scores/sec
		// MPLP 3 edge iters, 4 threads:	  8 sec            ??? scores/sec
		// MPLP 4 edge iters, 4 threads:	  5 sec            ??? scores/sec
		// MPLP 5 edge iters, 1 thread: 	9.5 sec            ??? scores/sec
		// MPLP 5 edge iters, 4 threads:	  4 sec            ??? scores/sec
		// MPLP 6 edge iters, 4 threads:	2.5 sec            ??? scores/sec
		// MPLP 7 edge iters, 4 threads:	2.3 sec            ??? scores/sec
		// MPLP 8 edge iters, 4 threads:	2.3 sec            ??? scores/sec
		
		assertThat(conf.getAssignments(), is(new int[] { 165, 17, 26, 0, 0, 4, 122 }));
		assertThat(conf.getScore(), isAbsolutely(-50.792885, 1e-6));
		
		// enumerate energy window if needed
		if (Ew > 0) {
			Stopwatch stopwatch = new Stopwatch().start();
			
			// simple way
			//tree.nextConfs(econf.getEnergy() + eW);
			
			final EnergyRange window = new EnergyRange(econf.getEnergy(), Ew);
			tree.getProgress().setGoalScore(window.getMax());
			
			// enumerate all confs in order of the scores, up to the estimate of the top of the energy window
			System.out.println("Enumerating energy window...");
			List<ScoredConf> scoredConfs = new ArrayList<>();
			scoredConfs.add(conf);
			Stopwatch minStopwatch = new Stopwatch().start();
			int indexToMinimizeNext = 1;
			while (true) {
				
				ScoredConf nextConf = tree.nextConf();
				scoredConfs.add(nextConf);
				if (nextConf.getScore() >= window.getMax()) {
					break;
				}
				
				// if we've been enumerating confs for a while, try a minimization to see if we get a smaller window
				if (minStopwatch.getTimeS() >= 2) {
					minStopwatch.stop();
					minStopwatch.start();
					
					EnergiedConf nextEconf = ecalc.calcEnergy(scoredConfs.get(indexToMinimizeNext++));
					
					boolean changed = window.updateMin(nextEconf.getEnergy());
					if (changed) {
						System.out.println(String.format("Updated energy window! remaining: %14.8f", window.getMax() - nextConf.getScore()));
						tree.getProgress().setGoalScore(window.getMax());
					}
				}
			}
			
			System.out.println("Enumerated energy window in " + stopwatch.stop().getTime(2));
		}
	}
	
	@SuppressWarnings("unused")
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
		index.numDefined = definedRCs.length;
		for (int i=0; i<index.numDefined; i++) {
			index.definedPos[i] = i;
			index.definedRCs[i] = definedRCs[i];
		}
		index.numUndefined = conf.getAssignments().length - definedRCs.length;
		for (int i=0; i<index.numUndefined; i++) {
			index.undefinedPos[i] = definedRCs.length + i;
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
		for (int i=0; i<index.numUndefined; i++) {
			int pos1 = index.undefinedPos[i];
			int rc1 = conf.getAssignments()[pos1];
			
			// assignment to pos1
			hscore += emat.getOneBody(pos1, rc1);

			// interactions with defined residues
			for (int j=0; j<index.numDefined; j++) {
				int pos2 = index.definedPos[j];
				int rc2 = index.definedRCs[j];
				hscore += emat.getPairwise(pos1, rc1, pos2, rc2);
			}

			// interactions with undefined residues
			for (int j=0; j<index.numUndefined; j++) {
				int pos2 = index.undefinedPos[j];
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
			for (int j=0; j<index.numDefined; j++) {
				int pos2 = index.definedPos[j];
				int rc2 = index.definedRCs[j];
				rcPicks.rcs[pos2] = rc2;
				rcPicks.energies[pos2] = search.emat.getPairwise(pos1, rc1, pos2, rc2);
			}

			// interactions with undefined residues
			for (int j=0; j<index.numUndefined; j++) {
				int pos2 = index.undefinedPos[j];
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
