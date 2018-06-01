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

package edu.duke.cs.osprey.mplp;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import org.junit.BeforeClass;
import org.junit.Test;

import edu.duke.cs.osprey.TestBase;
import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.astar.conf.linked.LinkedConfAStarNode;
import edu.duke.cs.osprey.astar.conf.scoring.TraditionalPairwiseHScorer;
import edu.duke.cs.osprey.astar.conf.scoring.mplp.MessageVars;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;

public class TestMessageVars extends TestBase {
	
	@BeforeClass
	public static void before() {
		initDefaultEnvironment();
	}
	
	private SearchProblem makeSearchProblemDagkRigid() {
		EnergyMatrixConfig emConfig = new EnergyMatrixConfig();
		emConfig.pdbPath = "examples/DAGK/2KDC.P.forOsprey.pdb";
		emConfig.numFlexible = 8;
		emConfig.addWtRots = true;
		emConfig.doMinimize = false;
		return makeSearchProblem(emConfig);
	}
	
	private void checkTotalEnergy(EnergyMatrix emat, RCs rcs, ConfIndex confIndex, MessageVars vars) {
		
		// make sure the vars total energy matches the traditional A* heuristic
		double tradEnergy = new TraditionalPairwiseHScorer(emat, rcs).calc(confIndex, rcs);
		assertThat(vars.getTotalEnergy(), isRelatively(tradEnergy));
	}
	
	@Test
	public void testTotalEnergyRootNodeWithPrecompute() {
		
		// get the RCs and a conf index
		SearchProblem search = makeSearchProblemDagkRigid();
		RCs rcs = new RCs(search.pruneMat);
		ConfIndex confIndex = new ConfIndex(rcs.getNumPos());
		
		// initialize message vars with traditional A* heuristic
		MessageVars vars = new MessageVars(rcs, confIndex, true);
		vars.initTraditionalAStar(search.emat);
		
		checkTotalEnergy(search.emat, rcs, confIndex, vars);
	}
	
	@Test
	public void testTotalEnergyRootNodeNoPrecompute() {
		
		// get the RCs and a conf index
		SearchProblem search = makeSearchProblemDagkRigid();
		RCs rcs = new RCs(search.pruneMat);
		ConfIndex confIndex = new ConfIndex(rcs.getNumPos());
		
		// initialize message vars with traditional A* heuristic
		MessageVars vars = new MessageVars(rcs, confIndex, false);
		vars.initTraditionalAStar(search.emat);
		
		checkTotalEnergy(search.emat, rcs, confIndex, vars);
	}
	
	@Test
	public void testTotalEnergyInternalNodeWithPrecompute() {
		
		// get the RCs and a conf index
		SearchProblem search = makeSearchProblemDagkRigid();
		RCs rcs = new RCs(search.pruneMat);
		ConfIndex confIndex = new ConfIndex(rcs.getNumPos());
		new LinkedConfAStarNode()
			.assign(0, 0)
			.assign(1, 0)
			.assign(2, 0)
		.index(confIndex);
		
		// initialize message vars with traditional A* heuristic
		MessageVars vars = new MessageVars(rcs, confIndex, true);
		vars.initTraditionalAStar(search.emat);
		
		checkTotalEnergy(search.emat, rcs, confIndex, vars);
	}
	
	@Test
	public void testTotalEnergyInternalNodeNoPrecompute() {
		
		// get the RCs and a conf index
		SearchProblem search = makeSearchProblemDagkRigid();
		RCs rcs = new RCs(search.pruneMat);
		ConfIndex confIndex = new ConfIndex(rcs.getNumPos());
		new LinkedConfAStarNode()
			.assign(0, 0)
			.assign(1, 0)
			.assign(2, 0)
		.index(confIndex);
		
		// initialize message vars with traditional A* heuristic
		MessageVars vars = new MessageVars(rcs, confIndex, false);
		vars.initTraditionalAStar(search.emat);
		
		checkTotalEnergy(search.emat, rcs, confIndex, vars);
	}
	
	@Test
	public void testInf() {
		
		// get the RCs and a conf index
		SearchProblem search = makeSearchProblemDagkRigid();
		RCs rcs = new RCs(search.pruneMat);
		ConfIndex confIndex = new ConfIndex(rcs.getNumPos());
		new LinkedConfAStarNode()
			.assign(0, 0)
			.assign(1, 0)
			.assign(2, 0)
		.index(confIndex);
		
		// make some zero messages
		MessageVars vars = new MessageVars(rcs, confIndex);
		
		assertThat(vars.getEnergy(0, 0), is(0.0));
		assertThat(vars.getEnergyWithout(0, 0, 1), is(0.0));
		assertThat(vars.getEnergyWithout(0, 0, 2), is(0.0));
		
		// set some inf
		vars.set(1, 0, 0, Double.POSITIVE_INFINITY);
		
		assertThat(vars.getEnergy(0, 0), is(Double.POSITIVE_INFINITY));
		assertThat(vars.getEnergyWithout(0, 0, 1), is(0.0));
		assertThat(vars.getEnergyWithout(0, 0, 2), is(Double.POSITIVE_INFINITY));
		
		// set some more inf
		vars.set(2, 0, 0, Double.POSITIVE_INFINITY);
		
		assertThat(vars.getEnergy(0, 0), is(Double.POSITIVE_INFINITY));
		assertThat(vars.getEnergyWithout(0, 0, 1), is(Double.POSITIVE_INFINITY));
		assertThat(vars.getEnergyWithout(0, 0, 2), is(Double.POSITIVE_INFINITY));
		
		// go back one
		vars.set(1, 0, 0, 0);
		assertThat(vars.getEnergy(0, 0), is(Double.POSITIVE_INFINITY));
		assertThat(vars.getEnergyWithout(0, 0, 1), is(Double.POSITIVE_INFINITY));
		assertThat(vars.getEnergyWithout(0, 0, 2), is(0.0));
		
		// go back to all zero
		vars.set(2, 0, 0, 0);
		assertThat(vars.getEnergy(0, 0), is(0.0));
		assertThat(vars.getEnergyWithout(0, 0, 1), is(0.0));
		assertThat(vars.getEnergyWithout(0, 0, 2), is(0.0));
	}
}
