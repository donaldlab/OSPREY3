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
import java.util.concurrent.TimeUnit;

import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.EnergyMatrixCalculator;
import edu.duke.cs.osprey.ematrix.epic.EPICSettings;
import edu.duke.cs.osprey.energy.MultiTermEnergyFunction;
import edu.duke.cs.osprey.pruning.PruningControl;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tools.ObjectIO;
import edu.duke.cs.osprey.tools.Stopwatch;
import edu.duke.cs.osprey.tupexp.LUTESettings;

public class PruningProfiling {
	
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
		// n=300 gives about 2.5e241 conformations to prune
		final int NumFlexible = 300;
		ArrayList<String> flexRes = new ArrayList<>();
		ArrayList<ArrayList<String>> allowedAAs = new ArrayList<>();
		for (int i=0; i<NumFlexible; i++) {
			flexRes.add(Integer.toString(i + 1));
			allowedAAs.add(new ArrayList<String>());
		}
		boolean addWt = true;
		boolean doMinimize = false;
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
			flexRes, allowedAAs, addWt, doMinimize, useEpic, new EPICSettings(), useTupleExpansion, new LUTESettings(),
			new DEEPerSettings(), moveableStrands, freeBBZones, useEllipses, useERef, addResEntropy, addWtRots, null, 
                        false, new ArrayList<>()
		);
		
		// compute/read the energy matrix
		File ematFile = new File(String.format("emat.%d.dat", NumFlexible));
		if (ematFile.exists()) {
			System.out.println("\nReading energy matrix...");
			search.emat = (EnergyMatrix)ObjectIO.readObject(ematFile.getAbsolutePath(), true);
		} else {
			System.out.println("\nComputing energy matrix...");
			EnergyMatrixCalculator emCalc = new EnergyMatrixCalculator(search.confSpace, search.shellResidues, useERef, addResEntropy);
			emCalc.calcPEM();
			search.emat = emCalc.getEMatrix();
			ObjectIO.writeObject(search.emat, ematFile.getAbsolutePath());
		}
		
		// init pruning
		double pruningInterval = 0;
		search.pruneMat = new PruningMatrix(search.confSpace, pruningInterval);
		boolean typeDep = false;
		double boundsThreshold = 100; // config default
		int algoOption = 3; // these two kind of do the same thing: turn on pair pruning
		boolean useFlags = true; //
		boolean useTriples = false;
		boolean preDacs = false;
		double stericThreshold = 100; // config default
		PruningControl pruner = new PruningControl(
			search, pruningInterval, typeDep, boundsThreshold, algoOption, 
			useFlags, useTriples, preDacs, useEpic, useTupleExpansion, stericThreshold
		);
		
		// notation below (trialN values in milliseconds):
		// numFlexPos: [trial1, trial2, trial2]
		
		// 2016-05-10
		// BEFORE OPTIMIZATIONS
		// 300: [18510, 18428, 18635]
		
		// optimize Pruner.canPrune(), PruningMatrix.isPruned()
		// 300: [10689, 10852, 10834] => 1.72x speedup
		
		// ugh... something happened to my computer after lunch, everything's slower now
		// re-benchmarking current code:
		// 300: [12973, 13275, 13109]
		
		// specialize TupleMatrixBoolean for PruningMatrix
		// 300: [12609, 12279, 12424] => 1.05x speedup over benchmark
		// not quite as dramatic an improvement as the EnergyMatrix double specialization, but I'll take it
		
		// more optimizations (particularly CPU cache optimizations)
		// 300: [11200, 11023, 10874] => 1.19x speedup over benchmark, 2.05x speedup over original
		
		System.out.println("\nPruning " + search.confSpace.getNumConformations().doubleValue() + " conformations...");
		Stopwatch stopwatch = new Stopwatch();
		stopwatch.start();
		pruner.prune();            
		stopwatch.stop();
		System.out.println("finished in " + stopwatch.getTime(TimeUnit.MILLISECONDS));
		
		/* TODO: check pruning accuracy automatically:
		
		Starting steric pruning.
		Pruned 82 in 1-body steric pruning
		Pruned 7745 in 2-body steric pruning
		Starting DEE cycle run: 0
		Starting pruning with GOLDSTEIN FOR 1 POSITION(S)
		Num pruned rot this run: 2275
		Num pruned pairs this run: 0

		Starting DEE cycle run: 1
		Starting pruning with GOLDSTEIN FOR 1 POSITION(S)
		Starting pruning with GOLDSTEIN FOR 2 POSITION(S)
		Num pruned rot this run: 0
		Num pruned pairs this run: 24861

		Starting DEE cycle run: 2
		Starting pruning with GOLDSTEIN FOR 1 POSITION(S)
		Num pruned rot this run: 73
		Num pruned pairs this run: 0

		Starting DEE cycle run: 3
		Starting pruning with GOLDSTEIN FOR 1 POSITION(S)
		Starting pruning with GOLDSTEIN FOR 2 POSITION(S)
		Num pruned rot this run: 0
		Num pruned pairs this run: 0
		*/
	}
}
