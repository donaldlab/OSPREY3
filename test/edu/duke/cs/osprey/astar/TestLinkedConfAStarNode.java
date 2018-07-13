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

import static edu.duke.cs.osprey.astar.Matchers.*;
import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import org.junit.Test;

import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.linked.LinkedConfAStarNode;

public class TestLinkedConfAStarNode {
	
	@Test
	public void indexRoot() {
		
		LinkedConfAStarNode node = new LinkedConfAStarNode();
		
		ConfIndex confIndex = new ConfIndex(5);
		node.index(confIndex);
		
		assertThat(confIndex.node, is(node));
		assertThat(confIndex.numPos, is(5));
		assertThat(confIndex.numDefined, is(0));
		assertThat(confIndex.numUndefined, is(5));
		assertThat(confIndex.undefinedPos, startsWith(new int[] { 0, 1, 2, 3, 4 }));
	}
	
	@Test
	public void indexChild0() {
		
		LinkedConfAStarNode node = new LinkedConfAStarNode()
			.assign(0, 5);
		
		ConfIndex confIndex = new ConfIndex(5);
		node.index(confIndex);
		
		assertThat(confIndex.node, is(node));
		assertThat(confIndex.numPos, is(5));
		assertThat(confIndex.numDefined, is(1));
		assertThat(confIndex.definedPos, startsWith(0));
		assertThat(confIndex.definedRCs, startsWith(5));
		assertThat(confIndex.numUndefined, is(4));
		assertThat(confIndex.undefinedPos, startsWith(1, 2, 3, 4));
	}
	
	@Test
	public void indexChild3() {
		
		LinkedConfAStarNode node = new LinkedConfAStarNode()
			.assign(3, 6);
		
		ConfIndex confIndex = new ConfIndex(5);
		node.index(confIndex);
		
		assertThat(confIndex.node, is(node));
		assertThat(confIndex.numPos, is(5));
		assertThat(confIndex.numDefined, is(1));
		assertThat(confIndex.definedPos, startsWith(3));
		assertThat(confIndex.definedRCs, startsWith(6));
		assertThat(confIndex.numUndefined, is(4));
		assertThat(confIndex.undefinedPos, startsWith(0, 1, 2, 4));
	}
	
	@Test
	public void indexChild03() {
		
		LinkedConfAStarNode node = new LinkedConfAStarNode()
			.assign(0, 5)
			.assign(3, 6);
		
		ConfIndex confIndex = new ConfIndex(5);
		node.index(confIndex);
		
		assertThat(confIndex.node, is(node));
		assertThat(confIndex.numPos, is(5));
		assertThat(confIndex.numDefined, is(2));
		assertThat(confIndex.definedPos, startsWith(0, 3));
		assertThat(confIndex.definedRCs, startsWith(5, 6));
		assertThat(confIndex.numUndefined, is(3));
		assertThat(confIndex.undefinedPos, startsWith(1, 2, 4));
	}
	
	@Test
	public void indexChild30() {
		
		LinkedConfAStarNode node = new LinkedConfAStarNode()
			.assign(3, 6)
			.assign(0, 5);
		
		ConfIndex confIndex = new ConfIndex(5);
		node.index(confIndex);
		
		assertThat(confIndex.node, is(node));
		assertThat(confIndex.numPos, is(5));
		assertThat(confIndex.numDefined, is(2));
		assertThat(confIndex.definedPos, startsWith(0, 3));
		assertThat(confIndex.definedRCs, startsWith(5, 6));
		assertThat(confIndex.numUndefined, is(3));
		assertThat(confIndex.undefinedPos, startsWith(1, 2, 4));
	}
	
	@Test
	public void indexChild01234() {
		
		LinkedConfAStarNode node = new LinkedConfAStarNode()
			.assign(0, 5)
			.assign(1, 7)
			.assign(2, 9)
			.assign(3, 6)
			.assign(4, 2);
		
		ConfIndex confIndex = new ConfIndex(5);
		node.index(confIndex);
		
		assertThat(confIndex.node, is(node));
		assertThat(confIndex.numPos, is(5));
		assertThat(confIndex.numDefined, is(5));
		assertThat(confIndex.definedPos, startsWith(0, 1, 2, 3, 4));
		assertThat(confIndex.definedRCs, startsWith(5, 7, 9, 6, 2));
		assertThat(confIndex.numUndefined, is(0));
	}
	
	@Test
	public void indexChild43210() {
		
		LinkedConfAStarNode node = new LinkedConfAStarNode()
			.assign(4, 2)
			.assign(3, 6)
			.assign(2, 9)
			.assign(1, 7)
			.assign(0, 5);
		
		ConfIndex confIndex = new ConfIndex(5);
		node.index(confIndex);
		
		assertThat(confIndex.node, is(node));
		assertThat(confIndex.numPos, is(5));
		assertThat(confIndex.numDefined, is(5));
		assertThat(confIndex.definedPos, startsWith(0, 1, 2, 3, 4));
		assertThat(confIndex.definedRCs, startsWith(5, 7, 9, 6, 2));
		assertThat(confIndex.numUndefined, is(0));
	}
}
