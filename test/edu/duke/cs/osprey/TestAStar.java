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

import static edu.duke.cs.osprey.tools.Log.log;
import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import com.google.common.collect.Lists;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tools.MathTools;
import org.junit.BeforeClass;
import org.junit.Test;

import edu.duke.cs.osprey.astar.ConfTree;
import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.astar.conf.order.DynamicHMeanAStarOrder;
import edu.duke.cs.osprey.astar.conf.order.StaticScoreHMeanAStarOrder;
import edu.duke.cs.osprey.astar.conf.scoring.MPLPPairwiseHScorer;
import edu.duke.cs.osprey.astar.conf.scoring.PairwiseGScorer;
import edu.duke.cs.osprey.astar.conf.scoring.TraditionalPairwiseHScorer;
import edu.duke.cs.osprey.astar.conf.scoring.mplp.EdgeUpdater;
import edu.duke.cs.osprey.astar.conf.scoring.mplp.NodeUpdater;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.externalMemory.ExternalMemory;

import java.util.Arrays;
import java.util.List;

public class TestAStar extends TestBase {
	
	@BeforeClass
	public static void before() {
		initDefaultEnvironment();
	}
	
	
	// RIGID TESTS (without pruning)
	
	private SearchProblem makeSearchProblemDagkRigid() {
		EnergyMatrixConfig emConfig = new EnergyMatrixConfig();
		emConfig.pdbPath = "examples/DAGK/2KDC.P.forOsprey.pdb";
		emConfig.numFlexible = 8;
		emConfig.addWtRots = true;
		emConfig.doMinimize = false;
		return makeSearchProblem(emConfig);
	}
	
	private void checkDagkRigid(ConfSearch tree, SearchProblem search) {
		ConfSearch.ScoredConf conf = tree.nextConf();
		assertThat(conf.getScore(), isRelatively(-78.78903544331037));
		assertThat(conf.getAssignments(), is(new int[] { 0, 6, 7, 0, 16, 1, 0, 6 }));
	}
	
	@Test
	public void testDagkRigidConfTree() {
		SearchProblem search = makeSearchProblemDagkRigid();
		
		ConfSearch tree = ConfTree.makeFull(search);
		
		checkDagkRigid(tree, search);
	}
	
	@Test
	public void testDagkRigidDynamicScoreOrderTraditional() {
		SearchProblem search = makeSearchProblemDagkRigid();
		
		RCs rcs = new RCs(search.pruneMat);
		ConfAStarTree tree = new ConfAStarTree.Builder(search.emat, rcs)
			.setCustom(
				new DynamicHMeanAStarOrder(),
				new PairwiseGScorer(search.emat),
				new TraditionalPairwiseHScorer(search.emat, rcs)
			).build();
		
		checkDagkRigid(tree, search);
	}
	
	@Test
	public void testDagkRigidStaticScoreOrderTraditional() {
		SearchProblem search = makeSearchProblemDagkRigid();
		
		RCs rcs = new RCs(search.pruneMat);
		ConfAStarTree tree = new ConfAStarTree.Builder(search.emat, rcs)
			.setCustom(
				new StaticScoreHMeanAStarOrder(),
				new PairwiseGScorer(search.emat),
				new TraditionalPairwiseHScorer(search.emat, rcs)
			).build();
		
		checkDagkRigid(tree, search);
	}
	
	@Test
	public void testDagkRigidDynamicOrderMPLP0Iter() {
		SearchProblem search = makeSearchProblemDagkRigid();
		
		ConfAStarTree tree = new ConfAStarTree.Builder(search.emat, search.pruneMat)
			.setCustom(
				new DynamicHMeanAStarOrder(),
				new PairwiseGScorer(search.emat),
				new MPLPPairwiseHScorer(null, search.emat, 0, 0.0001)
			).build();
		
		checkDagkRigid(tree, search);
	}
	
	@Test
	public void testDagkRigidStaticScoreOrderMPLP0Iter() {
		SearchProblem search = makeSearchProblemDagkRigid();
		
		ConfAStarTree tree = new ConfAStarTree.Builder(search.emat, search.pruneMat)
			.setCustom(
				new StaticScoreHMeanAStarOrder(),
				new PairwiseGScorer(search.emat),
				new MPLPPairwiseHScorer(null, search.emat, 0, 0.0001)
			).build();
		
		checkDagkRigid(tree, search);
	}
	
	@Test
	public void testDagkRigidStaticScoreOrderMPLPNode1Iter() {
		SearchProblem search = makeSearchProblemDagkRigid();
		
		ConfAStarTree tree = new ConfAStarTree.Builder(search.emat, search.pruneMat)
			.setCustom(
				new StaticScoreHMeanAStarOrder(),
				new PairwiseGScorer(search.emat),
				new MPLPPairwiseHScorer(new NodeUpdater(), search.emat, 1, 0.0001)
			).build();
		
		checkDagkRigid(tree, search);
	}
	
	@Test
	public void testDagkRigidStaticScoreOrderMPLPNode5Iter() {
		SearchProblem search = makeSearchProblemDagkRigid();
		
		ConfAStarTree tree = new ConfAStarTree.Builder(search.emat, search.pruneMat)
			.setCustom(
				new StaticScoreHMeanAStarOrder(),
				new PairwiseGScorer(search.emat),
				new MPLPPairwiseHScorer(new NodeUpdater(), search.emat, 5, 0.0001)
			).build();
		
		checkDagkRigid(tree, search);
	}
	
	@Test
	public void testDagkRigidStaticScoreOrderMPLPEdge1Iter() {
		SearchProblem search = makeSearchProblemDagkRigid();
		
		ConfAStarTree tree = new ConfAStarTree.Builder(search.emat, search.pruneMat)
			.setCustom(
				new StaticScoreHMeanAStarOrder(),
				new PairwiseGScorer(search.emat),
				new MPLPPairwiseHScorer(new EdgeUpdater(), search.emat, 1, 0.0001)
			).build();
		
		checkDagkRigid(tree, search);
	}
	
	@Test
	public void testDagkRigidStaticScoreOrderMPLPEdge20Iter() {
		SearchProblem search = makeSearchProblemDagkRigid();
		
		ConfAStarTree tree = new ConfAStarTree.Builder(search.emat, search.pruneMat)
			.setCustom(
				new StaticScoreHMeanAStarOrder(),
				new PairwiseGScorer(search.emat),
				new MPLPPairwiseHScorer(new EdgeUpdater(), search.emat, 20, 0.0001)
			).build();
		
		checkDagkRigid(tree, search);
	}
	
	@Test
	public void testDagkRigidDynamicOrderMPLPNode1Iter() {
		SearchProblem search = makeSearchProblemDagkRigid();
		
		ConfAStarTree tree = new ConfAStarTree.Builder(search.emat, search.pruneMat)
			.setCustom(
				new DynamicHMeanAStarOrder(),
				new PairwiseGScorer(search.emat),
				new MPLPPairwiseHScorer(new NodeUpdater(), search.emat, 1, 0.0001)
			).build();
		
		checkDagkRigid(tree, search);
	}
	
	
	// RIGID TESTS (with pruning)
	
	private void setPruning(SearchProblem search, int[] partialConf) {
		
		int numChosen = 0;
		for (int pos=0; pos<search.confSpace.numPos; pos++) {
			int chosenRc = partialConf[pos];
			if (chosenRc > 0) {
				numChosen++;
				
				// prune all but this rc
				for (int rc=0; rc<search.confSpace.posFlex.get(pos).RCs.size(); rc++) {
					if (rc != chosenRc) {
						search.pruneMat.setOneBody(pos, rc, true);
					}
				}
			}
		}
		
		assert (new RCs(search.pruneMat).getNumTrivialPos() == numChosen);
	}
	
	private SearchProblem makePrunedSearchProblemDagkRigid() {
		SearchProblem search = makeSearchProblemDagkRigid();
		setPruning(search, new int[] { -1, -1, 7, -1, 16, -1, 0, 6 });
		return search;
	}
	
	@Test
	public void testDagkRigidConfTreePruned() {
		SearchProblem search = makePrunedSearchProblemDagkRigid();
		
		ConfSearch tree = ConfTree.makeFull(search);
		
		checkDagkRigid(tree, search);
	}

	@Test
	public void testDagkRigidStaticScoreOrderMPLPNode1IterPruned() {
		SearchProblem search = makePrunedSearchProblemDagkRigid();
		
		ConfAStarTree tree = new ConfAStarTree.Builder(search.emat, search.pruneMat)
			.setCustom(
				new StaticScoreHMeanAStarOrder(),
				new PairwiseGScorer(search.emat),
				new MPLPPairwiseHScorer(new NodeUpdater(), search.emat, 1, 0.0001)
			).build();
		
		checkDagkRigid(tree, search);
	}
	
	
	// RIGID TESTS (with infinite energies)
	
	private SearchProblem makeSearchProblemDagkRigidInf() {
		SearchProblem search = makeSearchProblemDagkRigid();
		
		// add some infinite energies to see if any of the A* algorithms break
		search.emat = new EnergyMatrix(search.emat);
		for (int pos1=0; pos1<search.emat.getNumPos(); pos1++) {
			for (int pos2=0; pos2<pos1; pos2++) {
				search.emat.setPairwise(pos1, 0, pos2, 0, Double.POSITIVE_INFINITY);
			}
		}
		
		return search;
	}
	
	private void checkDagkRigidInf(ConfSearch tree, SearchProblem search) {
		ConfSearch.ScoredConf conf = tree.nextConf();
		assertThat(conf.getScore(), isRelatively(-73.66240935481423));
		assertThat(conf.getAssignments(), is(new int[] { 0, 6, 7, 16, 16, 1, 2, 6 }));
	}
	
	@Test
	public void testDagkRigidInfConfTree() {
		SearchProblem search = makeSearchProblemDagkRigidInf();
		
		ConfSearch tree = ConfTree.makeFull(search);
		
		checkDagkRigidInf(tree, search);
	}
	
	@Test
	public void testDagkRigidInfDynamicScoreOrderTraditional() {
		SearchProblem search = makeSearchProblemDagkRigidInf();
		
		RCs rcs = new RCs(search.pruneMat);
		ConfAStarTree tree = new ConfAStarTree.Builder(search.emat, rcs)
			.setCustom(
				new DynamicHMeanAStarOrder(),
				new PairwiseGScorer(search.emat),
				new TraditionalPairwiseHScorer(search.emat, rcs)
			).build();
		
		checkDagkRigidInf(tree, search);
	}
	
	@Test
	public void testDagkRigidInfStaticScoreOrderMPLP0Iter() {
		SearchProblem search = makeSearchProblemDagkRigidInf();
		
		ConfAStarTree tree = new ConfAStarTree.Builder(search.emat, search.pruneMat)
			.setCustom(
				new StaticScoreHMeanAStarOrder(),
				new PairwiseGScorer(search.emat),
				new MPLPPairwiseHScorer(null, search.emat, 0, 0.0001)
			).build();
		
		checkDagkRigidInf(tree, search);
	}
	
	@Test
	public void testDagkRigidInfStaticScoreOrderMPLPNode1Iter() {
		SearchProblem search = makeSearchProblemDagkRigidInf();
		
		ConfAStarTree tree = new ConfAStarTree.Builder(search.emat, search.pruneMat)
			.setCustom(
				new StaticScoreHMeanAStarOrder(),
				new PairwiseGScorer(search.emat),
				new MPLPPairwiseHScorer(new NodeUpdater(), search.emat, 1, 0.0001)
			).build();
		
		checkDagkRigidInf(tree, search);
	}
	
	@Test
	public void testDagkRigidInfStaticScoreOrderMPLPEdge1Iter() {
		SearchProblem search = makeSearchProblemDagkRigidInf();
		
		ConfAStarTree tree = new ConfAStarTree.Builder(search.emat, search.pruneMat)
			.setCustom(
				new StaticScoreHMeanAStarOrder(),
				new PairwiseGScorer(search.emat),
				new MPLPPairwiseHScorer(new EdgeUpdater(), search.emat, 1, 0.0001)
			).build();
		
		checkDagkRigidInf(tree, search);
	}
	
	
	// CONTINUOUS TESTS
	
	private SearchProblem makeSearchProblemDagkContinuous() {
		EnergyMatrixConfig emConfig = new EnergyMatrixConfig();
		emConfig.pdbPath = "examples/DAGK/2KDC.P.forOsprey.pdb";
		emConfig.numFlexible = 8;
		emConfig.addWtRots = true;
		emConfig.doMinimize = true;
		return makeSearchProblem(emConfig);
	}
	
	private void checkDagkContinuous(ConfSearch tree, SearchProblem search) {
		ConfSearch.ScoredConf conf = tree.nextConf();
		assertThat(conf.getScore(), isRelatively(-84.85105597709511));
		assertThat(conf.getAssignments(), is(new int[] { 0, 0, 7, 16, 16, 1, 1, 17 }));
	}
	
	@Test
	public void testDagkContinuousConfTree() {
		SearchProblem search = makeSearchProblemDagkContinuous();
		
		ConfSearch tree = ConfTree.makeFull(search);
		
		checkDagkContinuous(tree, search);
	}
	
	@Test
	public void testDagkContinuousDynamicScoreOrderTraditional() {
		SearchProblem search = makeSearchProblemDagkContinuous();
		
		RCs rcs = new RCs(search.pruneMat);
		ConfAStarTree tree = new ConfAStarTree.Builder(search.emat, rcs)
			.setCustom(
				new DynamicHMeanAStarOrder(),
				new PairwiseGScorer(search.emat),
				new TraditionalPairwiseHScorer(search.emat, rcs)
			).build();
		
		checkDagkContinuous(tree, search);
	}
	
	@Test
	public void testDagkContinuousStaticScoreOrderTraditional() {
		SearchProblem search = makeSearchProblemDagkContinuous();
		
		RCs rcs = new RCs(search.pruneMat);
		ConfAStarTree tree = new ConfAStarTree.Builder(search.emat, rcs)
			.setCustom(
				new StaticScoreHMeanAStarOrder(),
				new PairwiseGScorer(search.emat),
				new TraditionalPairwiseHScorer(search.emat, rcs)
			).build();
		
		checkDagkContinuous(tree, search);
	}
	
	@Test
	public void testDagkContinuousDynamicOrderMPLP0Iter() {
		SearchProblem search = makeSearchProblemDagkContinuous();
		
		ConfAStarTree tree = new ConfAStarTree.Builder(search.emat, search.pruneMat)
			.setCustom(
				new DynamicHMeanAStarOrder(),
				new PairwiseGScorer(search.emat),
				new MPLPPairwiseHScorer(null, search.emat, 0, 0.0001)
			).build();
		
		checkDagkContinuous(tree, search);
	}
	
	@Test
	public void testDagkContinuousStaticScoreOrderMPLP0Iter() {
		SearchProblem search = makeSearchProblemDagkContinuous();
		
		ConfAStarTree tree = new ConfAStarTree.Builder(search.emat, search.pruneMat)
			.setCustom(
				new StaticScoreHMeanAStarOrder(),
				new PairwiseGScorer(search.emat),
				new MPLPPairwiseHScorer(null, search.emat, 0, 0.0001)
			).build();
		
		checkDagkContinuous(tree, search);
	}
	
	@Test
	public void testDagkContinuousStaticScoreOrderMPLPNode1Iter() {
		SearchProblem search = makeSearchProblemDagkContinuous();
		
		ConfAStarTree tree = new ConfAStarTree.Builder(search.emat, search.pruneMat)
			.setCustom(
				new StaticScoreHMeanAStarOrder(),
				new PairwiseGScorer(search.emat),
				new MPLPPairwiseHScorer(new NodeUpdater(), search.emat, 1, 0.0001)
			).build();
		
		checkDagkContinuous(tree, search);
	}
	
	@Test
	public void testDagkContinuousStaticScoreOrderMPLPNode5Iter() {
		SearchProblem search = makeSearchProblemDagkContinuous();
		
		ConfAStarTree tree = new ConfAStarTree.Builder(search.emat, search.pruneMat)
			.setCustom(
				new StaticScoreHMeanAStarOrder(),
				new PairwiseGScorer(search.emat),
				new MPLPPairwiseHScorer(new NodeUpdater(), search.emat, 5, 0.0001)
			).build();
		
		checkDagkContinuous(tree, search);
	}
	
	@Test
	public void testDagkContinuousStaticScoreOrderMPLPEdge1Iter() {
		SearchProblem search = makeSearchProblemDagkContinuous();
		
		ConfAStarTree tree = new ConfAStarTree.Builder(search.emat, search.pruneMat)
			.setCustom(
				new StaticScoreHMeanAStarOrder(),
				new PairwiseGScorer(search.emat),
				new MPLPPairwiseHScorer(new EdgeUpdater(), search.emat, 1, 0.0001)
			).build();
		
		checkDagkContinuous(tree, search);
	}
	
	@Test
	public void testDagkContinuousStaticScoreOrderMPLPEdge20Iter() {
		SearchProblem search = makeSearchProblemDagkContinuous();
		
		ConfAStarTree tree = new ConfAStarTree.Builder(search.emat, search.pruneMat)
			.setCustom(
				new StaticScoreHMeanAStarOrder(),
				new PairwiseGScorer(search.emat),
				new MPLPPairwiseHScorer(new EdgeUpdater(), search.emat, 20, 0.0001)
			).build();
		
		checkDagkContinuous(tree, search);
	}
	
	@Test
	public void testDagkContinuousDynamicOrderMPLPNode1Iter() {
		SearchProblem search = makeSearchProblemDagkContinuous();
		
		ConfAStarTree tree = new ConfAStarTree.Builder(search.emat, search.pruneMat)
			.setCustom(
				new DynamicHMeanAStarOrder(),
				new PairwiseGScorer(search.emat),
				new MPLPPairwiseHScorer(new NodeUpdater(), search.emat, 1, 0.0001)
			).build();
		
		checkDagkContinuous(tree, search);
	}
	
	
	// EXTERNAL MEMORY TESTS
	
	@Test
	public void testExternalMemory() {
		SearchProblem search = makeSearchProblemDagkContinuous();
		
		ExternalMemory.use(16, () -> {
			ConfAStarTree tree = new ConfAStarTree.Builder(search.emat, search.pruneMat)
				.setTraditional()
				.useExternalMemory()
				.build();
			
			checkDagkContinuous(tree, search);
		});
	}

	@Test
	public void optimization() {

		Strand strand = new Strand.Builder(PDBIO.readResource("/1CC8.ss.pdb")).build();
		for (String resNum : Arrays.asList("A5", "A6", "A7")) {
			strand.flexibility.get(resNum)
				.setLibraryRotamers("ALA")
				.addWildTypeRotamers();
		}

		SimpleConfSpace confSpace = new SimpleConfSpace.Builder()
			.addStrand(strand)
			.build();

		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpace, new ForcefieldParams()).build()) {

			EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(confSpace, ecalc)
				.build()
				.calcEnergyMatrix();

			// get all the confs in ascending order
			List<ConfSearch.ScoredConf> ascendingConfs = new ConfAStarTree.Builder(emat, confSpace)
				.setTraditionalOpt(MathTools.Optimizer.Minimize)
				.build()
				.nextConfs(Double.POSITIVE_INFINITY);

			// make sure the sequences are in weakly ascending order
			for (int i=1; i<ascendingConfs.size(); i++) {
				assertThat(ascendingConfs.get(i-1).getScore(), lessThanOrEqualTo(ascendingConfs.get(i).getScore()));
			}

			// try maximization instead of minimization
			List<ConfSearch.ScoredConf> descendingConfs = new ConfAStarTree.Builder(emat, confSpace)
				.setTraditionalOpt(MathTools.Optimizer.Maximize)
				.build()
				.nextConfs(Double.NEGATIVE_INFINITY);

			// make sure the sequences are in weakly descending order
			for (int i=1; i<descendingConfs.size(); i++) {
				assertThat(descendingConfs.get(i-1).getScore(), greaterThanOrEqualTo(descendingConfs.get(i).getScore()));
			}

			// the two lists should match if we reverse the ascending confs
			assertThat(descendingConfs.size(), is(ascendingConfs.size()));
			for (int i=0; i<descendingConfs.size(); i++) {
				assertThat(descendingConfs.get(i).getScore(), isAbsolutely(ascendingConfs.get(ascendingConfs.size() - 1 - i).getScore()));
			}
		}
	}
}
