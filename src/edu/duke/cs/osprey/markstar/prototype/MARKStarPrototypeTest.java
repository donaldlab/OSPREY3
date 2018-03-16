package edu.duke.cs.osprey.markstar.prototype;

import org.junit.Test;

public class MARKStarPrototypeTest {
	
	@Test
	public void testMARKStarPrototype()
	{
		SimpleAStarSearch baseTree = new SimpleAStarSearch();
		RecursiveKStarBound bound = new RecursiveKStarBound(baseTree,0.99);
		bound.computeBound();
		System.out.println("Test!");
	}

}
