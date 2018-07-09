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

package edu.duke.cs.osprey.confspace;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import org.junit.BeforeClass;
import org.junit.Test;

import edu.duke.cs.osprey.TestBase;
import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;

public class TestConfSearchMultiSplitter extends TestBase {
	
	private static SearchProblem search;
	
	@BeforeClass
	public static void before() {
		initDefaultEnvironment();
		
		// make any search problem, doesn't matter
		EnergyMatrixConfig emConfig = new EnergyMatrixConfig();
		emConfig.pdbPath = "examples/DAGK/2KDC.P.forOsprey.pdb";
		emConfig.numFlexible = 2; // total confs in tree is 16
		emConfig.addWtRots = true;
		emConfig.doMinimize = false;
		search = makeSearchProblem(emConfig);
	}
	
	private ConfSearch makeTree() {
		
		// make any search tree, doesn't matter
		return new ConfAStarTree.Builder(search.emat, search.pruneMat).build();
	}
	
	private void assertConfs(ConfSearch.MultiSplitter.Stream stream, ConfSearch expectedTree, int numConfs) {
		
		for (int i=0; i<numConfs; i++) {
			ScoredConf expectedConf = expectedTree.nextConf();
			ScoredConf observedConf = stream.nextConf();
			if (observedConf == null) {
				assertThat(observedConf, is(nullValue()));
			} else {
				assertThat(observedConf.getAssignments(), is(expectedConf.getAssignments()));
			}
		}
	}
	
	@Test
	public void testOne() {
		
		ConfSearch tree = makeTree();
		ConfSearch.MultiSplitter splitter = new ConfSearch.MultiSplitter(tree);
		ConfSearch.MultiSplitter.Stream stream = splitter.makeStream();
		ConfSearch check1 = makeTree();

		assertThat(splitter.getBufferSize(), is(0));
		assertConfs(stream, check1, 16);
		assertThat(splitter.getBufferSize(), is(0));
	}
	
	@Test
	public void testTwo() {
		
		ConfSearch tree = makeTree();
		ConfSearch.MultiSplitter splitter = new ConfSearch.MultiSplitter(tree);
		ConfSearch.MultiSplitter.Stream stream1 = splitter.makeStream();
		ConfSearch.MultiSplitter.Stream stream2 = splitter.makeStream();
		ConfSearch check1 = makeTree();
		ConfSearch check2 = makeTree();

		assertThat(splitter.getBufferSize(), is(0));
		assertConfs(stream1, check1, 16);
		assertThat(splitter.getBufferSize(), is(16));
		assertConfs(stream2, check2, 16);
		assertThat(splitter.getBufferSize(), is(0));
	}
	
	@Test
	public void testOneEnd() {
		
		ConfSearch tree = makeTree();
		ConfSearch.MultiSplitter splitter = new ConfSearch.MultiSplitter(tree);
		ConfSearch.MultiSplitter.Stream stream = splitter.makeStream();
		ConfSearch check1 = makeTree();

		assertThat(splitter.getBufferSize(), is(0));
		assertConfs(stream, check1, 17);
		assertThat(splitter.getBufferSize(), is(0));
	}
	
	@Test
	public void testTwoEnd() {
		
		ConfSearch tree = makeTree();
		ConfSearch.MultiSplitter splitter = new ConfSearch.MultiSplitter(tree);
		ConfSearch.MultiSplitter.Stream stream1 = splitter.makeStream();
		ConfSearch.MultiSplitter.Stream stream2 = splitter.makeStream();
		ConfSearch check1 = makeTree();
		ConfSearch check2 = makeTree();

		assertThat(splitter.getBufferSize(), is(0));
		assertConfs(stream1, check1, 30);
		assertThat(splitter.getBufferSize(), is(16));
		assertConfs(stream2, check2, 30);
		assertThat(splitter.getBufferSize(), is(0));
	}
	
	@Test
	public void testTwoTrade() {
		
		ConfSearch tree = makeTree();
		ConfSearch.MultiSplitter splitter = new ConfSearch.MultiSplitter(tree);
		ConfSearch.MultiSplitter.Stream stream1 = splitter.makeStream();
		ConfSearch.MultiSplitter.Stream stream2 = splitter.makeStream();
		ConfSearch check1 = makeTree();
		ConfSearch check2 = makeTree();
		
		assertThat(splitter.getBufferSize(), is(0));
		
		assertConfs(stream1, check1, 4);
		assertThat(splitter.getBufferSize(), is(4));
		assertConfs(stream2, check2, 4);
		assertThat(splitter.getBufferSize(), is(0));
		
		assertConfs(stream1, check1, 4);
		assertThat(splitter.getBufferSize(), is(4));
		assertConfs(stream2, check2, 4);
		assertThat(splitter.getBufferSize(), is(0));
		
		assertConfs(stream1, check1, 4);
		assertThat(splitter.getBufferSize(), is(4));
		assertConfs(stream2, check2, 4);
		assertThat(splitter.getBufferSize(), is(0));
		
		assertConfs(stream1, check1, 4);
		assertThat(splitter.getBufferSize(), is(4));
		assertConfs(stream2, check2, 4);
		assertThat(splitter.getBufferSize(), is(0));
	}
	
	@Test
	public void testTwoLeapfrog() {
		
		ConfSearch tree = makeTree();
		ConfSearch.MultiSplitter splitter = new ConfSearch.MultiSplitter(tree);
		ConfSearch.MultiSplitter.Stream stream1 = splitter.makeStream();
		ConfSearch.MultiSplitter.Stream stream2 = splitter.makeStream();
		ConfSearch check1 = makeTree();
		ConfSearch check2 = makeTree();
		
		assertThat(splitter.getBufferSize(), is(0));
		
		assertConfs(stream1, check1, 4);
		assertThat(splitter.getBufferSize(), is(4));
		assertConfs(stream2, check2, 8);
		assertThat(splitter.getBufferSize(), is(4));
		
		assertConfs(stream1, check1, 8);
		assertThat(splitter.getBufferSize(), is(4));
		assertConfs(stream2, check2, 4);
		assertThat(splitter.getBufferSize(), is(0));
	}
	
	@Test
	public void testTwoLeapfrogEnd() {
		
		ConfSearch tree = makeTree();
		ConfSearch.MultiSplitter splitter = new ConfSearch.MultiSplitter(tree);
		ConfSearch.MultiSplitter.Stream stream1 = splitter.makeStream();
		ConfSearch.MultiSplitter.Stream stream2 = splitter.makeStream();
		ConfSearch check1 = makeTree();
		ConfSearch check2 = makeTree();
		
		assertThat(splitter.getBufferSize(), is(0));
		
		assertConfs(stream1, check1, 4);
		assertThat(splitter.getBufferSize(), is(4));
		assertConfs(stream2, check2, 20);
		assertThat(splitter.getBufferSize(), is(12));
	}
	
	@Test
	public void testTwoRemoveOne() {
		
		ConfSearch tree = makeTree();
		ConfSearch.MultiSplitter splitter = new ConfSearch.MultiSplitter(tree);
		ConfSearch.MultiSplitter.Stream stream1 = splitter.makeStream();
		ConfSearch.MultiSplitter.Stream stream2 = splitter.makeStream();
		ConfSearch check1 = makeTree();
		ConfSearch check2 = makeTree();
	
		assertThat(splitter.getBufferSize(), is(0));
		assertConfs(stream1, check1, 8);
		assertThat(splitter.getBufferSize(), is(8));
		assertConfs(stream2, check2, 4);
		stream2.close();
		assertThat(splitter.getBufferSize(), is(0));
		assertConfs(stream1, check1, 8);
		assertThat(splitter.getBufferSize(), is(0));
	}
}
