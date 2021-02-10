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

package edu.duke.cs.osprey;

import java.io.File;
import java.util.ArrayList;

import edu.duke.cs.osprey.astar.conf.ConfAStarNode;
import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.astar.conf.order.DynamicHMeanAStarOrder;
import edu.duke.cs.osprey.astar.conf.order.StaticScoreHMeanAStarOrder;
import edu.duke.cs.osprey.astar.conf.scoring.MPLPPairwiseHScorer;
import edu.duke.cs.osprey.astar.conf.scoring.PairwiseGScorer;
import edu.duke.cs.osprey.astar.conf.scoring.TraditionalPairwiseHScorer;
import edu.duke.cs.osprey.astar.conf.scoring.mplp.MPLPUpdater;
import edu.duke.cs.osprey.astar.conf.scoring.mplp.NodeUpdater;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.EnergyMatrixCalculator;
import edu.duke.cs.osprey.ematrix.epic.EPICSettings;
import edu.duke.cs.osprey.energy.MultiTermEnergyFunction;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tools.ObjectIO;
import edu.duke.cs.osprey.tupexp.LUTESettings;

public class MPLPPLayground {

	public static void main(String[] args)
	throws Exception {
		
		// check the cwd
		String path = new File("").getAbsolutePath();
		if (!path.endsWith("examples/DAGK")) {
			throw new Error("This profiler was designed to run in the examples/DAGK folder\n\tcwd: " + path);
		}

		// load configuration
		ConfigFileParser cfp = new ConfigFileParser();
		cfp.loadData();
		
		// multi-thread the energy function
		MultiTermEnergyFunction.setNumThreads(4);
		
		// init a conf space with lots of flexible residues, but no mutations
		final int NumFlexible = 16;
		//final int NumFlexible = 25;
		ArrayList<String> flexRes = new ArrayList<>();
		ArrayList<ArrayList<String>> allowedAAs = new ArrayList<>();
		for (int i=0; i<NumFlexible; i++) {
			flexRes.add(Integer.toString(i + 1));
			allowedAAs.add(new ArrayList<String>());
		}
		boolean addWt = true;
		boolean doMinimize = true;
		boolean useEpic = false;
		boolean useTupleExpansion = false;
		boolean useEllipses = false;
		boolean useERef = false;
		boolean addResEntropy = false;
		boolean addWtRots = true;
		ArrayList<String[]> moveableStrands = new ArrayList<String[]>();
		ArrayList<String[]> freeBBZones = new ArrayList<String[]>();
		SearchProblem search = new SearchProblem(
			"energyMatrixProfiling",
			"2KDC.P.forOsprey.pdb", 
			flexRes, allowedAAs, addWt, doMinimize, useEpic, new EPICSettings(), useTupleExpansion,  new LUTESettings(),
			new DEEPerSettings(), moveableStrands, freeBBZones, useEllipses, useERef, addResEntropy, addWtRots, null, 
                        false, new ArrayList<>()
		);
		
		// compute the energy matrix
		File ematFile = new File(String.format("emat.%s%d.dat", doMinimize ? "min." : "", NumFlexible));
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
		
		// don't bother with pruning, set all to unpruned
		search.pruneMat = new PruningMatrix(search.confSpace, search.emat.getPruningInterval());
		RCs rcs = new RCs(search.pruneMat);
		
		// define a partial conformation (set all flexible residues undefined for now)
		ConfIndex confIndex = new ConfIndex(NumFlexible);
		confIndex.numDefined = 0;
		confIndex.numUndefined = NumFlexible;
		for (int i=0; i<NumFlexible; i++) {
			confIndex.undefinedPos[i] = i;
		}
		
		/* TEMP
		// assign some positions
		confIndex = new ConfIndex(confIndex, NumFlexible/2, 0);
		confIndex = new ConfIndex(confIndex, NumFlexible/3, 0);
		confIndex = new ConfIndex(confIndex, NumFlexible*7/8, 0);
		*/
		
		// config the different heuristics
		TraditionalPairwiseHScorer tradHScorer = new TraditionalPairwiseHScorer(search.emat, rcs);
		//MPLPUpdater mplpUpdater = new EdgeUpdater();
		MPLPUpdater mplpUpdater = new NodeUpdater();
		//int numIterations = 0;
		int numIterations = 1;
		//int numIterations = 10;
		MPLPPairwiseHScorer mplpHScorer = new MPLPPairwiseHScorer(mplpUpdater, search.emat, numIterations, 0.000001);
		
		double tradHScore = tradHScorer.calc(confIndex, rcs);
		System.out.println(String.format("Trad H Score: %16.12f", tradHScore));
		
		double mplpHScore = mplpHScorer.calc(confIndex, rcs);
		System.out.println(String.format("MPLP H Score: %16.12f", mplpHScore));
		
		// get the real min bound conf using a trusted A* implementation
		ConfAStarTree tree = new ConfAStarTree.Builder(search.emat, rcs)
			.setCustom(
				new DynamicHMeanAStarOrder(),
				new PairwiseGScorer(search.emat),
				new TraditionalPairwiseHScorer(search.emat, rcs)
			).build();
		ConfSearch.ScoredConf minBoundConf = tree.nextConf();
		System.out.println(String.format("min bound e (trad):  %16.12f", minBoundConf.getScore()));
		
		// get the min bound conf using MPLP
		tree = new ConfAStarTree.Builder(search.emat, rcs)
			.setCustom(
				new StaticScoreHMeanAStarOrder(),
				new PairwiseGScorer(search.emat),
				mplpHScorer
			).build();
		ConfSearch.ScoredConf minBoundConfMplp = tree.nextConf();
		System.out.println(String.format("min bound e (MPLP):  %16.12f", minBoundConfMplp.getScore()));
		
		double minBoundMinimizedEnergy = search.minimizedEnergy(minBoundConf.getAssignments());
		System.out.println(String.format("min energy:   %16.12f", minBoundMinimizedEnergy));
	}
}
