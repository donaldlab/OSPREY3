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

package edu.duke.cs.osprey.partcr;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import java.util.List;

import org.junit.BeforeClass;
import org.junit.Test;

import edu.duke.cs.osprey.TestBase;
import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.confspace.RC;
import edu.duke.cs.osprey.partcr.RCSplits.RCInfo;
import edu.duke.cs.osprey.partcr.splitters.BinaryRCSplitter;

public class TestRCSplits extends TestBase {
	
	@BeforeClass
	public static void before() {
		initDefaultEnvironment();
	}
	
	private ConfSpace makeConfSpace() {
		
		// need to get RCs from anywhere, doesn't matter where
		
		EnergyMatrixConfig emConfig = new EnergyMatrixConfig();
		emConfig.pdbPath = "examples/DAGK/2KDC.P.forOsprey.pdb";
		emConfig.numFlexible = 2;
		emConfig.addWtRots = false;
		emConfig.doMinimize = true;
		return makeSearchProblem(emConfig).confSpace;
	}
	
	@Test
	public void testBase() {
		ConfSpace confSpace = makeConfSpace();
		RCSplits splits = new RCSplits(confSpace);
		
		RCInfo info;
	
		// just spot check a few RCs
		info = splits.getRCInfo(0, 0);
		assertThat(info.isSplit(), is(false));
		assertThat(info.getRCs(), contains(0));
		assertThat(splits.getRCInfo(confSpace.posFlex.get(0).RCs.get(0)), is(info));
		
		info = splits.getRCInfo(1, 5);
		assertThat(info.isSplit(), is(false));
		assertThat(info.getRCs(), contains(5));
		assertThat(splits.getRCInfo(confSpace.posFlex.get(1).RCs.get(5)), is(info));
	}
	
	@Test
	public void testOneSplit() {
		ConfSpace confSpace = makeConfSpace();
		RCSplits splits = new RCSplits(confSpace);
		
		// make the split
		RC parentRC = confSpace.posFlex.get(1).RCs.get(5);
		List<RC> childRCs = new BinaryRCSplitter().split(1, parentRC);
		assertThat(childRCs, is(not(nullValue())));
		splits.split(parentRC, childRCs);
		
		// renumber the split RCs
		childRCs.get(0).RCIndex = parentRC.RCIndex;
		childRCs.get(1).RCIndex = confSpace.posFlex.get(1).RCs.size() - 1;
		
		// check the results
		RCInfo info;
		
		// nothing should have changed at 0,0
		info = splits.getRCInfo(0, 0);
		assertThat(info.isSplit(), is(false));
		assertThat(info.getRCs(), containsInAnyOrder(0));
		assertThat(splits.getRCInfo(confSpace.posFlex.get(0).RCs.get(0)), is(info));
		
		// check changes at 1,5
		info = splits.getRCInfo(1, 5);
		assertThat(info.isSplit(), is(true));
		assertThat(info.getRCs(), containsInAnyOrder(
			childRCs.get(0).RCIndex,
			childRCs.get(1).RCIndex
		));
		assertThat(splits.getRCInfo(confSpace.posFlex.get(1).RCs.get(5)), is(nullValue()));
		assertThat(splits.getRCInfo(childRCs.get(0)), is(info));
		assertThat(splits.getRCInfo(childRCs.get(1)), is(info));
	}
	
	@Test
	public void testTwoSplits() {
		ConfSpace confSpace = makeConfSpace();
		RCSplits splits = new RCSplits(confSpace);
		
		// make first split
		RC parentRC = confSpace.posFlex.get(1).RCs.get(5);
		List<RC> childRCs = new BinaryRCSplitter().split(1, parentRC);
		assertThat(childRCs, is(not(nullValue())));
		splits.split(parentRC, childRCs);
		
		// renumber the split RCs
		childRCs.get(0).RCIndex = parentRC.RCIndex;
		childRCs.get(1).RCIndex = confSpace.posFlex.get(1).RCs.size();
		
		// make second split
		List<RC> childChildRCs = new BinaryRCSplitter().split(1, childRCs.get(0));
		assertThat(childChildRCs, is(not(nullValue())));
		splits.split(childRCs.get(0), childChildRCs);
		
		// renumber the split RCs
		childChildRCs.get(0).RCIndex = childRCs.get(0).RCIndex;
		childChildRCs.get(1).RCIndex = confSpace.posFlex.get(1).RCs.size() - 1;
		
		// check the results
		RCInfo info;
		
		// nothing should have changed at 0,0
		info = splits.getRCInfo(0, 0);
		assertThat(info.isSplit(), is(false));
		assertThat(info.getRCs(), containsInAnyOrder(0));
		assertThat(splits.getRCInfo(confSpace.posFlex.get(0).RCs.get(0)), is(info));
		
		// check changes at 1,5
		info = splits.getRCInfo(1, 5);
		assertThat(info.isSplit(), is(true));
		assertThat(info.getRCs(), containsInAnyOrder(
			childRCs.get(1).RCIndex,
			childChildRCs.get(0).RCIndex,
			childChildRCs.get(1).RCIndex
		));
		assertThat(splits.getRCInfo(confSpace.posFlex.get(1).RCs.get(5)), is(nullValue()));
		assertThat(splits.getRCInfo(childRCs.get(0)), is(nullValue()));
		assertThat(splits.getRCInfo(childRCs.get(1)), is(info));
		assertThat(splits.getRCInfo(childChildRCs.get(0)), is(info));
		assertThat(splits.getRCInfo(childChildRCs.get(1)), is(info));
	}
}
