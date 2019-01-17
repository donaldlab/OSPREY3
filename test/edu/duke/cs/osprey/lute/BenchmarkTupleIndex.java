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

package edu.duke.cs.osprey.lute;

import edu.duke.cs.osprey.confspace.*;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tools.Stopwatch;
import edu.duke.cs.osprey.tupexp.TupleIndexMatrix;
import org.apache.commons.collections4.ListUtils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

import static edu.duke.cs.osprey.tools.Log.log;

public class BenchmarkTupleIndex {

	public static void main(String[] args) {

		Strand strand = new Strand.Builder(PDBIO.readResource("/1CC8.ss.pdb")).build();
		List<String> mutableResNums = Arrays.asList("A38", "A39");
		//List<String> mutableResNums = Arrays.asList("A38", "A39", "A40"); // this works great on pairs
		//List<String> mutableResNums = Arrays.asList("A38", "A39", "A40", "A41"); // NOTE: need triples to get a good LUTE fit here
		for (String resNum : mutableResNums) {
			strand.flexibility.get(resNum)
				.setLibraryRotamers(Strand.WildType, "VAL", "LEU")
				.addWildTypeRotamers()
				.setContinuous();
		}

		// add some flexibility
		List<String> flexibleResNums = Arrays.asList("A40", "A41", "A42", "A43", "A44");
		for (String resNum : flexibleResNums) {
			strand.flexibility.get(resNum)
				.setLibraryRotamers(Strand.WildType)
				.addWildTypeRotamers()
				.setContinuous();
		}

		SimpleConfSpace confSpace = new SimpleConfSpace.Builder()
			.addStrand(strand)
			.build();

		// collect all the pair and triple tuples
		List<RCTuple> densePairs = new ArrayList<>();
		List<RCTuple> denseTriples = new ArrayList<>();
		for (SimpleConfSpace.Position pos1 : confSpace.positions) {

			for (SimpleConfSpace.Position pos2 : confSpace.positions) {
				if (pos2.index >= pos1.index) {
					continue;
				}

				for (SimpleConfSpace.ResidueConf rc1 : pos1.resConfs) {
					for (SimpleConfSpace.ResidueConf rc2 : pos2.resConfs) {

						densePairs.add(new RCTuple(
							pos2.index, rc2.index,
							pos1.index, rc1.index
						));
					}
				}

				for (SimpleConfSpace.Position pos3 : confSpace.positions) {
					if (pos3.index >= pos2.index) {
						continue;
					}

					for (SimpleConfSpace.ResidueConf rc1 : pos1.resConfs) {
						for (SimpleConfSpace.ResidueConf rc2 : pos2.resConfs) {
							for (SimpleConfSpace.ResidueConf rc3 : pos3.resConfs) {

								denseTriples.add(new RCTuple(
									pos3.index, rc3.index,
									pos2.index, rc2.index,
									pos1.index, rc1.index
								));
							}
						}
					}
				}
			}
		}

		Random rand = new Random(12345);

		// pick a random 1% of the dense triples
		int numSparseTriples = denseTriples.size()/100;
		List<RCTuple> sparseTriples = new ArrayList<>(numSparseTriples);
		while (sparseTriples.size() < numSparseTriples) {
			sparseTriples.add(denseTriples.get(rand.nextInt(denseTriples.size())));
		}

		// pick a bunch of random conformations
		final int numConfs = 2000000;
		List<int[]> confs = new ArrayList<>(numConfs);
		for (int i=0; i<numConfs; i++) {
			int[] conf = new int[confSpace.positions.size()];
			Arrays.fill(conf, Conf.Unassigned);
			for (SimpleConfSpace.Position pos : confSpace.positions) {
				conf[pos.index] = rand.nextInt(pos.resConfs.size());
			}
			confs.add(conf);
		}

		benchmark("dense", confSpace, ListUtils.union(densePairs, denseTriples), confs);
		benchmark("sparse", confSpace, ListUtils.union(densePairs, sparseTriples), confs);
	}

	private static void benchmark(String name, SimpleConfSpace confSpace, List<RCTuple> tuples, List<int[]> confs) {

		TupleIndexMatrix tmat = new TupleIndexMatrix(confSpace.getNumPos(), confSpace.getNumResConfsByPos(), 0.0);
		tmat.fill(-1);
		TupleTree<Integer> ttree = new TupleTree<>();
		TuplesIndex tindex = new TuplesIndex(confSpace);
		for (RCTuple tuple : tuples) {
			tindex.appendTuple(tuple);
		}

		// populate the indices with the tuples
		for (int t=0; t<tuples.size(); t++) {
			RCTuple tuple = tuples.get(t);

			tmat.setHigherOrder(tuple, t);
			ttree.put(tuple, t);
		}

		// benchmark
		log("%s:", name);

		Stopwatch tmatSw = new Stopwatch().start();
		for (int[] conf : confs) {
			tmat.calcSampleTuples(conf);
		}
		log("\ttmat      %s", tmatSw.stop().getTime(2));

		Stopwatch ttreeSw = new Stopwatch().start();
		for (int[] conf : confs) {
			ttree.forEachIn(conf, (tuple, index) -> {});
		}
		log("\tttree     %s    %.2fx speedup", ttreeSw.stop().getTime(2), tmatSw.getTimeS()/ttreeSw.getTimeS());

		Stopwatch tindexSw = new Stopwatch().start();
		for (int[] conf : confs) {
			tindex.forEachIn(conf, false, true, (index) -> {});
		}
		log("\ttindex    %s    %.2fx speedup", tindexSw.stop().getTime(2), tmatSw.getTimeS()/tindexSw.getTimeS());
	}
}
