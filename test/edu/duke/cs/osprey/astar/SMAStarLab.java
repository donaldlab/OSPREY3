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

import static edu.duke.cs.osprey.TestBase.isAbsolutely;
import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.Conf;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tools.Stopwatch;

import java.util.*;

import static edu.duke.cs.osprey.tools.Log.log;


public class SMAStarLab {

	public static void main(String[] args) {

		Molecule mol = PDBIO.readResource("/1CC8.ss.pdb");
		Strand strand = new Strand.Builder(mol).build();
		//List<String> resNums = Arrays.asList("A2", "A3", "A4");
		//List<String> resNums = Arrays.asList("A2", "A3", "A4", "A5", "A6", "A7");
		List<String> resNums = Arrays.asList("A2", "A3", "A4", "A5", "A6", "A7", "A8", "A9", "A10", "A11", "A12");
		for (String resNum : resNums) {
			strand.flexibility.get(resNum).setLibraryRotamers("VAL");
		}

		SimpleConfSpace confSpace = new SimpleConfSpace.Builder()
			.addStrand(strand)
			.build();

		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpace, new ForcefieldParams())
			 .setParallelism(Parallelism.makeCpu(8))
			 .build()
		) {

			EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(confSpace, ecalc)
				.build()
				.calcEnergyMatrix();

			RCs rcs = new RCs(confSpace);

			// enumerate the confs using A*
			ConfAStarTree astar = new ConfAStarTree.Builder(emat, rcs)
				.setTraditional()
				.build();
			Stopwatch astarStopwatch = new Stopwatch().start();
			List<ConfSearch.ScoredConf> astarConfs = astar.nextConfs(Double.POSITIVE_INFINITY);
			astarStopwatch.stop();

			// enumerate the confs using SMA*
			//SMAStar smastar = new SMAStar(emat, rcs, rcs.getNumPos() + 1);
			ConfAStarTree smastar = new ConfAStarTree.Builder(emat, rcs)
				//.setMaxNumNodes(rcs.getNumPos() + 1)
				.setMaxNumNodes(1000000)
				.setTraditional()
				.build();
			Stopwatch smastarStopwatch = new Stopwatch().start();
			List<ConfSearch.ScoredConf> smastarConfs = smastar.nextConfs(Double.POSITIVE_INFINITY);
			smastarStopwatch.stop();

			dumpConfs("A*", astarConfs);
			dumpConfs("SMA*", smastarConfs);

			log("  A* finished in %s", astarStopwatch.getTime(2));
			log("SMA* finished in %s", smastarStopwatch.getTime(2));

			checkConfs(astarConfs, smastarConfs);
		}
	}

	private static void dumpConfs(String name, List<ConfSearch.ScoredConf> confs) {
		log("%s confs: %d", name, confs.size());
		int size = Math.min(100, confs.size());
		for (int i=0; i<size; i++) {
			ConfSearch.ScoredConf conf = confs.get(i);
			log("\t%14.4f   %s", conf.getScore(), Conf.toString(conf.getAssignments()));
		}
	}

	private static void checkConfs(List<ConfSearch.ScoredConf> expectedConfs, List<ConfSearch.ScoredConf> observedConfs) {

		assertThat(observedConfs.size(), is(expectedConfs.size()));

		for (int i=0; i<expectedConfs.size(); i++) {
			ConfSearch.ScoredConf exp = expectedConfs.get(i);
			ConfSearch.ScoredConf obs = observedConfs.get(i);
			assertThat(obs.getAssignments(), is(exp.getAssignments()));
			assertThat(obs.getScore(), isAbsolutely(exp.getScore(), 1e-10));
		}
	}
}
