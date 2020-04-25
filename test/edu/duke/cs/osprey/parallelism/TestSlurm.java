package edu.duke.cs.osprey.parallelism;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import org.junit.Test;

import java.util.Arrays;


public class TestSlurm {

	@Test
	public void parseNodes() {

		check(1, "apple", 1, "1", "apple");
		check(1, "apple", 2, "2", "apple", "apple");
		check(1, "apple", 3, "3", "apple", "apple", "apple");

		check(2, "apple,pear", 2, "1,1", "apple", "pear");
		check(2, "apple,pear", 2, "1(x2)", "apple", "pear");
		check(2, "apple,pear", 4, "2,2", "apple", "apple", "pear", "pear");

		check(1, "beef[steak]", 1, "1", "beefsteak");
		check(2, "beef[steak,pie]", 2, "1(x2)", "beefsteak", "beefpie");

		check(1, "node[42]", 1, "1", "node42");
		check(3, "node[1,2,3]", 3, "1(x3)", "node1", "node2", "node3");
		check(3, "node[1-3]", 3, "1(x3)", "node1", "node2", "node3");

		check(4, "node[1-2,5-6]", 4, "1(x4)", "node1", "node2", "node5", "node6");
		check(4, "node[1-2,5-6]", 5, "1(x2),2,1", "node1", "node2", "node5", "node5", "node6");
		check(4, "node[1-2,5-6]", 6, "1(x2),2(x2)", "node1", "node2", "node5", "node5", "node6", "node6");

		check(4, "foo[1-2],bar[5-6]", 4, "1(x4)", "foo1", "foo2", "bar5", "bar6");
		check(4, "foo[1,2],bar[5,6]", 4, "1,1,1,1", "foo1", "foo2", "bar5", "bar6");

		check(10, "grisman-[37-40],jerry[1-6]", 10, "1(x10)",
			"grisman-37", "grisman-38", "grisman-39", "grisman-40",
			"jerry1", "jerry2", "jerry3", "jerry4", "jerry5", "jerry6"
		);
	}

	private static void check(int numNodes, String nodelist, int numTasks, String tasksPerNode, String ... nodes) {
		assertThat(numTasks, is(nodes.length));
		assertThat(Slurm.parseNodes(numNodes, nodelist, numTasks, tasksPerNode), is(Arrays.asList(nodes)));
	}
}
