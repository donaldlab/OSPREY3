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
