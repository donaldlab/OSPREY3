package edu.duke.cs.osprey.confspace;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
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
			strand.flexibility.get(resNum).setLibraryRotamers(Strand.WildType);
		}

		return new SimpleConfSpace.Builder()
			.addStrand(strand)
			.build();
	}

	@Test
	public void makeUnassigned() {

		Sequence s = Sequence.makeUnassigned(confSpace);

		assertThat(s.get("A5"), is(nullValue()));
		assertThat(s.get("A6"), is(nullValue()));

		assertThat(s.isWildType(), is(false));
		assertThat(s.countAssignments(), is(0));
		assertThat(s.isFullyAssigned(), is(false));
		assertThat(s.countMutations(), is(0));
	}

	@Test
	public void makeWildType() {

		Sequence s = Sequence.makeWildType(confSpace);

		assertThat(s.get("A5"), is("LYS"));
		assertThat(s.get("A6"), is("HIE"));

		assertThat(s.isWildType(), is(true));
		assertThat(s.countAssignments(), is(2));
		assertThat(s.isFullyAssigned(), is(true));
		assertThat(s.countMutations(), is(0));
	}

	@Test
	public void assignments() {

		Sequence s = Sequence.makeWildType(confSpace);

		Iterator<Sequence.Assignment> iter = s.iterator();

		Sequence.Assignment a = iter.next();
		assertThat(a.getPos(), is(confSpace.getPositionOrThrow("A5")));
		assertThat(a.getResNum(), is("A5"));
		assertThat(a.getResType(), is("LYS"));
		assertThat(a.getWildType(), is("LYS"));
		assertThat(a.isWildType(), is(true));
		assertThat(a.isAssigned(), is(true));
		assertThat(a.isMutated(), is(false));

		a = iter.next();
		assertThat(a.getPos(), is(confSpace.getPositionOrThrow("A6")));
		assertThat(a.getResNum(), is("A6"));
		assertThat(a.getResType(), is("HIE"));
		assertThat(a.getWildType(), is("HIE"));
		assertThat(a.isWildType(), is(true));
		assertThat(a.isAssigned(), is(true));
		assertThat(a.isMutated(), is(false));

		assertThat(iter.hasNext(), is(false));
	}

	@Test
	public void iteratorSubset() {

		Sequence s = Sequence.makeWildType(confSpace);

		Iterator<Sequence.Assignment> iter = s.iterator(Arrays.asList(confSpace.getPositionOrThrow("A6")));

		assertThat(iter.next().getResNum(), is("A6"));
		assertThat(iter.hasNext(), is(false));
	}

	@Test
	public void fillWildType() {

		Sequence s = Sequence.makeUnassigned(confSpace)
			.set("A6", "HIE")
			.fillWildType();

		assertThat(s, is(Sequence.makeWildType(confSpace)));
	}

	@Test
	public void mutation() {

		Sequence s = Sequence.makeWildType(confSpace);

		Sequence m = s.makeMutatedSequence("A6", "ALA");

		assertThat(m, is(not(s)));
		assertThat(m.get("A5"), is("LYS"));
		assertThat(m.get("A6"), is("ALA"));
		assertThat(m.isMutated("A5"), is(false));
		assertThat(m.isMutated("A6"), is(true));
		assertThat(m.isWildType(), is(false));
		assertThat(m.countAssignments(), is(2));
		assertThat(m.isFullyAssigned(), is(true));
		assertThat(m.countMutations(), is(1));
	}

	@Test
	public void equals() {

		assertThat(Sequence.makeUnassigned(confSpace), is(Sequence.makeUnassigned(confSpace)));
		assertThat(Sequence.makeUnassigned(confSpace).set("A5", "ALA"), is(Sequence.makeUnassigned(confSpace).set("A5", "ALA")));

		assertThat(Sequence.makeUnassigned(confSpace), is(not(Sequence.makeUnassigned(confSpace).set("A5", "ALA"))));
		assertThat(Sequence.makeUnassigned(confSpace).set("A5", "VAL"), is(not(Sequence.makeUnassigned(confSpace).set("A5", "ALA"))));
	}

	@Test
	public void equalsDifferentConfSpace() {

		assertThat(
			Sequence.makeWildType(makeConfSpace("A4", "A6"))
				.set("A4", "LYS")
				.set("A6", "HIE"),
			is(not(Sequence.makeWildType(confSpace)))
		);
	}

	@Test
	public void filter() {

		Sequence s = Sequence.makeWildType(confSpace);

		Sequence f = s.filter(makeConfSpace("A6"));

		assertThat(
			f.stream()
				.map((a) -> a.getResNum())
				.collect(Collectors.toList()),
			is(Arrays.asList("A6"))
		);
	}

	@Test
	public void testToString() {

		Sequence s = Sequence.makeWildType(confSpace);

		assertThat(s.toString(Sequence.Renderer.ResNum), is("A5 A6"));
		assertThat(s.toString(Sequence.Renderer.ResType), is("LYS HIE"));
		assertThat(s.toString(Sequence.Renderer.ResTypeMutations), is("lys hie"));
		assertThat(s.toString(Sequence.Renderer.Assignment), is("A5=LYS A6=HIE"));
		assertThat(s.toString(Sequence.Renderer.AssignmentMutations), is("A5=lys A6=hie"));

		assertThat(
			s.toString(Sequence.Renderer.ResNum, null, Arrays.asList(
				confSpace.getPositionOrThrow("A6"),
				confSpace.getPositionOrThrow("A5")
			)),
			is("A6 A5")
		);

		s = Sequence.makeWildType(makeConfSpace("A7", "A42"));

		assertThat(s.toString(Sequence.Renderer.ResNum), is("A7 A42"));
		assertThat(s.toString(Sequence.Renderer.ResNum, 5), is("A7    A42  "));

		assertThat(s.getMaxResNumLength(), is(3));
	}
}
