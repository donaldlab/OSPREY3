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

import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tools.Streams;
import org.junit.BeforeClass;
import org.junit.Test;

import java.util.Arrays;
import java.util.Iterator;
import java.util.stream.Collectors;

public class TestSequence {

	private static Molecule mol;
	private static SimpleConfSpace confSpace;

	@BeforeClass
	public static void beforeClass() {

		mol = PDBIO.readFile("examples/1CC8/1CC8.ss.pdb");
		confSpace = makeConfSpace("A5", "A6");

		assertThat(confSpace.positions.get(0).resNum, is("A5"));
		assertThat(confSpace.positions.get(0).resFlex.wildType, is("LYS"));

		assertThat(confSpace.positions.get(1).resNum, is("A6"));
		assertThat(confSpace.positions.get(1).resFlex.wildType, is("HIE"));
	}

	private static SimpleConfSpace makeConfSpace(String ... resNums) {

		Strand strand = new Strand.Builder(mol).build();
		for (String resNum : resNums) {
			strand.flexibility.get(resNum).setLibraryRotamers(Strand.WildType, "ALA");
		}

		return new SimpleConfSpace.Builder()
			.addStrand(strand)
			.build();
	}

	@Test
	public void makeUnassigned() {

		Sequence s = confSpace.makeUnassignedSequence();

		assertThat(s.get("A5"), is(nullValue()));
		assertThat(s.get("A6"), is(nullValue()));

		assertThat(s.isWildType(), is(false));
		assertThat(s.countAssignments(), is(0));
		assertThat(s.isFullyAssigned(), is(false));
		assertThat(s.countMutations(), is(0));
	}

	@Test
	public void makeWildType() {

		Sequence s = confSpace.makeWildTypeSequence();

		assertThat(s.get("A5").name, is("LYS"));
		assertThat(s.get("A6").name, is("HIE"));

		assertThat(s.isWildType(), is(true));
		assertThat(s.countAssignments(), is(2));
		assertThat(s.isFullyAssigned(), is(true));
		assertThat(s.countMutations(), is(0));
	}

	@Test
	public void assignments() {

		Sequence s = confSpace.makeWildTypeSequence();

		Iterator<Sequence.Assignment> iter = s.assignments().iterator();

		Sequence.Assignment a = iter.next();
		assertThat(a.pos, is(confSpace.seqSpace.getPositionOrThrow("A5")));
		assertThat(a.getResNum(), is("A5"));
		assertThat(a.getResType().name, is("LYS"));
		assertThat(a.getWildType().name, is("LYS"));
		assertThat(a.isWildType(), is(true));
		assertThat(a.isAssigned(), is(true));
		assertThat(a.isMutated(), is(false));

		a = iter.next();
		assertThat(a.pos, is(confSpace.seqSpace.getPositionOrThrow("A6")));
		assertThat(a.getResNum(), is("A6"));
		assertThat(a.getResType().name, is("HIE"));
		assertThat(a.getWildType().name, is("HIE"));
		assertThat(a.isWildType(), is(true));
		assertThat(a.isAssigned(), is(true));
		assertThat(a.isMutated(), is(false));

		assertThat(iter.hasNext(), is(false));
	}

	@Test
	public void fillWildType() {

		Sequence s = confSpace.makeUnassignedSequence();
		s.set("A6", "HIE");
		s.fillWildType();

		assertThat(s, is(confSpace.makeWildTypeSequence()));
	}

	@Test
	public void mutation() {

		Sequence s = confSpace.makeWildTypeSequence();

		Sequence m = s.copy().set("A6", "ALA");

		assertThat(m, is(not(s)));
		assertThat(m.get("A5").name, is("LYS"));
		assertThat(m.get("A6").name, is("ALA"));
		assertThat(m.isMutated("A5"), is(false));
		assertThat(m.isMutated("A6"), is(true));
		assertThat(m.isWildType(), is(false));
		assertThat(m.countAssignments(), is(2));
		assertThat(m.isFullyAssigned(), is(true));
		assertThat(m.countMutations(), is(1));
	}

	@Test
	public void equals() {

		assertThat(confSpace.makeUnassignedSequence(), is(confSpace.makeUnassignedSequence()));
		assertThat(confSpace.makeUnassignedSequence().set("A5", "ALA"), is(confSpace.makeUnassignedSequence().set("A5", "ALA")));

		assertThat(confSpace.makeUnassignedSequence(), is(not(confSpace.makeUnassignedSequence().set("A5", "ALA"))));
		assertThat(confSpace.makeUnassignedSequence().set("A5", "LYS"), is(not(confSpace.makeUnassignedSequence().set("A5", "ALA"))));
	}

	@Test
	public void equalsDifferentConfSpace() {

		assertThat(
			makeConfSpace("A4", "A6").makeWildTypeSequence()
				.set("A4", "ILE")
				.set("A6", "HIE"),
			is(not(confSpace.makeWildTypeSequence()))
		);
	}

	@Test
	public void filter() {

		Sequence s = confSpace.makeWildTypeSequence();

		Sequence f = s.filter(makeConfSpace("A6").seqSpace);

		assertThat(
			Streams.of(f.assignments())
				.map((a) -> a.getResNum())
				.collect(Collectors.toList()),
			is(Arrays.asList("A6"))
		);
	}

	@Test
	public void testToString() {

		Sequence s = confSpace.makeWildTypeSequence();

		assertThat(s.toString(Sequence.Renderer.ResNum), is("A5 A6"));
		assertThat(s.toString(Sequence.Renderer.ResType), is("LYS HIE"));
		assertThat(s.toString(Sequence.Renderer.ResTypeMutations), is("lys hie"));
		assertThat(s.toString(Sequence.Renderer.Assignment), is("A5=LYS A6=HIE"));
		assertThat(s.toString(Sequence.Renderer.AssignmentMutations), is("A5=lys A6=hie"));

		assertThat(
			s.toString(Sequence.Renderer.ResNum, null, Arrays.asList("A6", "A5")),
			is("A6 A5")
		);

		s = makeConfSpace("A7", "A42").makeWildTypeSequence();

		assertThat(s.toString(Sequence.Renderer.ResNum), is("A7 A42"));
		assertThat(s.toString(Sequence.Renderer.ResNum, 5), is("A7    A42  "));

		assertThat(s.calcCellSize(), is(3));
	}
}
