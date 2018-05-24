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

package edu.duke.cs.osprey.externalMemory;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import edu.duke.cs.osprey.astar.conf.ConfAStarNode;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.tpie.EntrySize;
import edu.duke.cs.tpie.serialization.SerializingDoublePriorityQueue;
import org.junit.Test;

import java.lang.reflect.Field;
import java.util.Arrays;
import java.util.stream.Collectors;
import java.util.stream.IntStream;


public class TestAStarSerialization {

	@Test
	public void astarNode1x10() {
		ExternalMemory.use(16, () -> {

			Queue<EMConfAStarNode> q = makeQueue(10); // 21 bytes
			assertQueue(q, EntrySize.Bytes32, AssignmentsSerializer.Encoding.Byte);

			q.push(makeNode(7, 4.2, 7.9, 1));

			assertNode(q.poll(), 7, 4.2, 7.9, 1);
		});
	}

	@Test
	public void astarNode2x10() {
		ExternalMemory.use(16, () -> {

			Queue<EMConfAStarNode> q = makeQueue(10, 10); // 22 bytes
			assertQueue(q, EntrySize.Bytes32, AssignmentsSerializer.Encoding.Byte);

			q.push(makeNode(7, 4.2, 7.9, 7, 2));

			assertNode(q.poll(), 7, 4.2, 7.9, 7, 2);
		});
	}

	@Test
	public void astarNode12x120() {
		ExternalMemory.use(16, () -> {

			Queue<EMConfAStarNode> q = makeQueue(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120); // 32 bytes
			assertQueue(q, EntrySize.Bytes32, AssignmentsSerializer.Encoding.Byte);

			q.push(makeNode(7, 4.2, 7.9, 9, 19, 29, 39, 49, 59, 69, 79, 89, 99, 109, 119));

			assertNode(q.poll(), 7, 4.2, 7.9, 9, 19, 29, 39, 49, 59, 69, 79, 89, 99, 109, 119);
		});
	}

	@Test
	public void astarNode13x10() {
		ExternalMemory.use(16, () -> {

			Queue<EMConfAStarNode> q = makeQueue(10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10); // 33 bytes
			assertQueue(q, EntrySize.Bytes64, AssignmentsSerializer.Encoding.Byte);

			q.push(makeNode(8, 0.3, 4.0, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3));
			q.push(makeNode(2, 3.8, 2.4, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6));
			q.push(makeNode(7, 4.2, 7.9, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, 9, 8, 7));
			q.push(makeNode(4, 9.2, 5.7, 0, 0, 3, 0, 0, 5, 0, 0, 7, 0, 0, 9, 0));

			assertNode(q.poll(), 8, 0.3, 4.0, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3); // 4.3
			assertNode(q.poll(), 2, 3.8, 2.4, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6); // 6.2
			assertNode(q.poll(), 7, 4.2, 7.9, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, 9, 8, 7); // 12.1
			assertNode(q.poll(), 4, 9.2, 5.7, 0, 0, 3, 0, 0, 5, 0, 0, 7, 0, 0, 9, 0); // 14.9
		});
	}

	@Test
	public void astarNode1x128() {
		ExternalMemory.use(16, () -> {

			Queue<EMConfAStarNode> q = makeQueue(128); // 21 bytes
			assertQueue(q, EntrySize.Bytes32, AssignmentsSerializer.Encoding.Byte);

			q.push(makeNode(7, 4.2, 7.9, 0));
			q.push(makeNode(2, 3.8, 2.4, 1));
			q.push(makeNode(4, 9.2, 5.7, 126));
			q.push(makeNode(8, 0.3, 4.0, 127));

			assertNode(q.poll(), 8, 0.3, 4.0, 127); // 4.3
			assertNode(q.poll(), 2, 3.8, 2.4, 1); // 6.2
			assertNode(q.poll(), 7, 4.2, 7.9, 0); // 12.1
			assertNode(q.poll(), 4, 9.2, 5.7, 126); // 14.9
		});
	}

	@Test(expected = RuntimeException.class)
	public void astarNode1x128TooBig() {
		ExternalMemory.use(16, () -> {

			Queue<EMConfAStarNode> q = makeQueue(128); // 21 bytes
			assertQueue(q, EntrySize.Bytes32, AssignmentsSerializer.Encoding.Byte);

			q.push(makeNode(7, 4.2, 7.9, 128));
		});
	}

	@Test
	public void astarNode1x129() {
		ExternalMemory.use(16, () -> {
			Queue<EMConfAStarNode> q = makeQueue(129); // 22 bytes
			assertQueue(q, EntrySize.Bytes32, AssignmentsSerializer.Encoding.Short);
		});
	}

	@Test
	public void astarNode1x32768() {
		ExternalMemory.use(16, () -> {

			Queue<EMConfAStarNode> q = makeQueue(32768); // 22 bytes
			assertQueue(q, EntrySize.Bytes32, AssignmentsSerializer.Encoding.Short);

			q.push(makeNode(7, 4.2, 7.9, 0));
			q.push(makeNode(2, 3.8, 2.4, 1));
			q.push(makeNode(4, 9.2, 5.7, 32766));
			q.push(makeNode(8, 0.3, 4.0, 32767));

			assertNode(q.poll(), 8, 0.3, 4.0, 32767); // 4.3
			assertNode(q.poll(), 2, 3.8, 2.4, 1); // 6.2
			assertNode(q.poll(), 7, 4.2, 7.9, 0); // 12.1
			assertNode(q.poll(), 4, 9.2, 5.7, 32766); // 14.9
		});
	}

	@Test(expected = RuntimeException.class)
	public void astarNode1x32768TooBig() {
		ExternalMemory.use(16, () -> {

			Queue<EMConfAStarNode> q = makeQueue(32768); // 22 bytes
			assertQueue(q, EntrySize.Bytes32, AssignmentsSerializer.Encoding.Short);

			q.push(makeNode(7, 4.2, 7.9, 32768));
		});
	}

	@Test
	public void astarNode1x32769() {
		ExternalMemory.use(16, () -> {
			Queue<EMConfAStarNode> q = makeQueue(32769); // 24 bytes
			assertQueue(q, EntrySize.Bytes32, AssignmentsSerializer.Encoding.Int);
		});
	}

	@Test
	public void astarNode9x200() {
		ExternalMemory.use(16, () -> {

			Queue<EMConfAStarNode> q = makeQueue(200, 200, 200, 200, 200, 200, 200, 200, 200); // 38 bytes
			assertQueue(q, EntrySize.Bytes64, AssignmentsSerializer.Encoding.Short);

			q.push(makeNode(8, 0.3, 4.0,  10,  20,  30,  40,  50,  60,  70,  80,  90));
			q.push(makeNode(2, 3.8, 2.4,   0,   0,   0,   0,   0,   0,   0,   0,   0));
			q.push(makeNode(7, 4.2, 7.9, 199, 199, 199, 199, 199, 199, 199, 199, 199));
			q.push(makeNode(4, 9.2, 5.7,   0,   9,  10,  19,  99, 100, 109, 159, 199));

			assertNode(q.poll(), 8, 0.3, 4.0,  10,  20,  30,  40,  50,  60,  70,  80,  90); // 4.3
			assertNode(q.poll(), 2, 3.8, 2.4,   0,   0,   0,   0,   0,   0,   0,   0,   0); // 6.2
			assertNode(q.poll(), 7, 4.2, 7.9, 199, 199, 199, 199, 199, 199, 199, 199, 199); // 12.1
			assertNode(q.poll(), 4, 9.2, 5.7,   0,   9,  10,  19,  99, 100, 109, 159, 199); // 14.9
		});
	}

	@Test
	public void astarNode22x200() {
		ExternalMemory.use(16, () -> {

			Queue<EMConfAStarNode> q = makeQueue(
				200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200,
				200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200
			); // 64 bytes
			assertQueue(q, EntrySize.Bytes64, AssignmentsSerializer.Encoding.Short);

			q.push(makeNode(8, 0.3, 4.0,  10,  20,  30,  40,  50,  60,  70,  80,  90,   4,   2,  10,  20,  30,  40,  50,  60,  70,  80,  90,   5,   5));
			q.push(makeNode(2, 3.8, 2.4,   0,   0,   0,   0,   0,   0,   0,   0,   0,   4,   2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   5,   5));
			q.push(makeNode(7, 4.2, 7.9, 199, 199, 199, 199, 199, 199, 199, 199, 199, 199, 199, 199, 199, 199, 199, 199, 199, 199, 199, 199, 199, 199));
			q.push(makeNode(4, 9.2, 5.7,   0,   9,  10,  19,  99, 100, 109, 159, 199,   4,   2,   0,   9,  10,  19,  99, 100, 109, 159, 199,   5,   5));

			assertNode(q.poll(), 8, 0.3, 4.0,  10,  20,  30,  40,  50,  60,  70,  80,  90,   4,   2,  10,  20,  30,  40,  50,  60,  70,  80,  90,   5,   5); // 4.3
			assertNode(q.poll(), 2, 3.8, 2.4,   0,   0,   0,   0,   0,   0,   0,   0,   0,   4,   2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   5,   5); // 6.2
			assertNode(q.poll(), 7, 4.2, 7.9, 199, 199, 199, 199, 199, 199, 199, 199, 199, 199, 199, 199, 199, 199, 199, 199, 199, 199, 199, 199, 199, 199); // 12.1
			assertNode(q.poll(), 4, 9.2, 5.7,   0,   9,  10,  19,  99, 100, 109, 159, 199,   4,   2,   0,   9,  10,  19,  99, 100, 109, 159, 199,   5,   5); // 14.9
		});
	}

	@Test
	public void astarNode1x200NonContiguousRCs() {
		ExternalMemory.use(16, () -> {

			Queue<EMConfAStarNode> q = makeQueue(new int[][] {
				{ 1, 2, 3, 4, 5 },
				{ 201, 202, 203, 204, 205 }
			}); // 24 bytes
			assertQueue(q, EntrySize.Bytes32, AssignmentsSerializer.Encoding.Short);

			q.push(makeNode(8, 0.3, 4.0, 1, 205));
			q.push(makeNode(2, 3.8, 2.4, 3, 202));
			q.push(makeNode(7, 4.2, 7.9, 5, 201));
			q.push(makeNode(4, 9.2, 5.7, 4, 204));

			assertNode(q.poll(), 8, 0.3, 4.0, 1, 205); // 4.3
			assertNode(q.poll(), 2, 3.8, 2.4, 3, 202); // 6.2
			assertNode(q.poll(), 7, 4.2, 7.9, 5, 201); // 12.1
			assertNode(q.poll(), 4, 9.2, 5.7, 4, 204); // 14.9
		});
	}

	@SuppressWarnings("unchecked")
	private static Queue<EMConfAStarNode> makeQueue(int ... rcSizes) {

		// make RCs where each pos has a N RCs in the range [0,N)
		RCs rcs = new RCs(IntStream.of(rcSizes)
			.mapToObj((rcSize) -> {
				return IntStream.range(0, rcSize)
					.boxed()
					.collect(Collectors.toList());
			})
			.collect(Collectors.toList())
		);

		// make the queue (and do ridiculous casting, because Java's type system is dumb)
		Queue<? extends ConfAStarNode> q = new EMConfAStarFactory().makeQueue(rcs);
		return (Queue<EMConfAStarNode>)q;
	}

	@SuppressWarnings("unchecked")
	private static Queue<EMConfAStarNode> makeQueue(int[][] rcVals) {

		// make RCs where each pos has a N RCs in the range [0,N)
		RCs rcs = new RCs(Arrays.stream(rcVals)
			.map((rcsAtPos) -> {
				return IntStream.of(rcsAtPos)
					.boxed()
					.collect(Collectors.toList());
			})
			.collect(Collectors.toList())
		);

		// make the queue (and do ridiculous casting, because Java's type system is dumb)
		Queue<? extends ConfAStarNode> q = new EMConfAStarFactory().makeQueue(rcs);
		return (Queue<EMConfAStarNode>)q;
	}


	private static EMConfAStarNode makeNode(int level, double gscore, double hscore, int ... assignments) {
		EMConfAStarNode node = new EMConfAStarNode(assignments.length);
		node.setLevel(level);
		node.setGScore(gscore);
		node.setHScore(hscore);
		System.arraycopy(assignments, 0, node.getConf(), 0, assignments.length);
		return node;
	}

	private static void assertNode(EMConfAStarNode node, int level, double gscore, double hscore, int ... assignments) {
		assertThat(node.getLevel(), is(level));
		assertThat(node.getGScore(), is(gscore));
		assertThat(node.getHScore(), is(hscore));
		assertThat(node.getConf(), is(assignments));
	}

	private static void assertQueue(Queue<EMConfAStarNode> q, EntrySize entrySize, AssignmentsSerializer.Encoding encoding) {

		SerializingDoublePriorityQueue<?> pq;
		try {
			// hack out the wrapped queue using reflection
			Field field = q.getClass().getDeclaredField("q");
			field.setAccessible(true);
			pq = (SerializingDoublePriorityQueue<?>)field.get(q);
		} catch (Exception ex) {
			throw new Error(ex);
		}

		AssignmentsSerializer serializer = (AssignmentsSerializer)pq.serializer;
		assertThat(serializer.entrySize, is(entrySize));
		assertThat(serializer.encoding, is(encoding));
	}
}
