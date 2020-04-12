package edu.duke.cs.osprey.confspace.compiled;

import edu.duke.cs.osprey.confspace.SeqSpace;
import edu.duke.cs.osprey.tools.FileTools;
import org.junit.Test;

import java.util.Arrays;
import java.util.List;
import java.util.function.Consumer;
import java.util.stream.Collectors;

import static org.hamcrest.Matchers.containsInAnyOrder;
import static org.hamcrest.Matchers.is;
import static org.hamcrest.core.IsNull.notNullValue;
import static org.junit.Assert.assertThat;


public class TestConfSpace {

	@Test
	public void check2RL0() {

		ConfSpace complex = ConfSpace.fromBytes(FileTools.readResourceBytes("/confSpaces/2RL0.complex.ccsx"));
		ConfSpace chainA = ConfSpace.fromBytes(FileTools.readResourceBytes("/confSpaces/2RL0.A.ccsx"));
		ConfSpace chainG = ConfSpace.fromBytes(FileTools.readResourceBytes("/confSpaces/2RL0.G.ccsx"));

		Consumer<ConfSpace> assertG = (confSpace) -> {
			assertPosition(confSpace, "649 PHE", "PHE", 29, Arrays.asList("PHE:5",  "TYR:8",  "ALA:1", "VAL:3", "ILE:7", "LEU:5"));
			assertPosition(confSpace, "650 ASP", "ASP", 14, Arrays.asList("ASP:6",  "GLU:8"));
			assertPosition(confSpace, "651 GLU", "GLU", 14, Arrays.asList("GLU:9",  "ASP:5"));
		};

		Consumer<ConfSpace> assertA = (confSpace) -> {
			assertPosition(confSpace, "156 PHE", "PHE", 29, Arrays.asList("PHE:5",  "TYR:8",  "ALA:1", "VAL:3", "ILE:7", "LEU:5"));
			assertPosition(confSpace, "172 LYS", "LYS", 41, Arrays.asList("LYS:28", "ASP:5",  "GLU:8"));
			assertPosition(confSpace, "192 ILE", "ILE", 29, Arrays.asList("ILE:8",  "ALA:1",  "VAL:3", "LEU:5", "PHE:4", "TYR:8"));
			assertPosition(confSpace, "193 THR", "THR", 44, Arrays.asList("THR:19", "SER:18", "ASN:7"));
		};

		// check chain G
		assertThat(chainG.positions.length, is(3));
		assertG.accept(chainG);
		assertThat(chainG.countSingles(), is(57));
		assertThat(chainG.countPairs(), is(1008));

		// check chain A
		assertThat(chainA.positions.length, is(4));
		assertA.accept(chainA);
		assertThat(chainA.countSingles(), is(143));
		assertThat(chainA.countPairs(), is(7575));

		// check the complex
		assertThat(complex.positions.length, is(7));
		assertG.accept(complex);
		assertA.accept(complex);
		assertThat(complex.countSingles(), is(200));
		assertThat(complex.countPairs(), is(16734));

		// TODO: check atom counts and positions are the same?
	}

	private static void assertPosition(ConfSpace confSpace, String name, String wildType, int numConfs, List<String> descriptions) {

		ConfSpace.Pos confPos = confSpace.findPos(name);
		assertThat("No position with name " + name, confPos, is(notNullValue()));
		SeqSpace.Position seqPos = confSpace.seqSpace.getPositionOrThrow(name);

		// check the name
		assertThat(confPos.name, is(name));
		assertThat(seqPos.resNum, is(name));

		// check the wild type
		assertThat(confPos.wildType, is(wildType));
		assertThat(seqPos.wildType.name, is(wildType));

		// check the types
		String[] resTypes = descriptions.stream()
			.map(desc -> desc.split(":")[0])
			.toArray(String[]::new);
		assertThat(
			Arrays.stream(confPos.confs)
				.map(conf -> conf.type)
				.collect(Collectors.toSet()), // collapse duplicates
			containsInAnyOrder(resTypes)
		);
		assertThat(
			seqPos.resTypes.stream()
				.map(resType -> resType.name)
				.collect(Collectors.toList()),
			containsInAnyOrder(resTypes)
		);

		// count the conformations for each type
		int[] numConfsByDesc = descriptions.stream()
			.mapToInt(desc -> Integer.parseInt(desc.split(":")[1]))
			.toArray();
		for (int i=0; i<descriptions.size(); i++) {
			int fi = i; // javac is dumb ...

			List<ConfSpace.Conf> resConfs = Arrays.stream(confPos.confs)
				.filter(conf -> conf.type.equals(resTypes[fi]))
				.collect(Collectors.toList());

			assertThat(resConfs.size(), is(numConfsByDesc[i]));
		}

		// check the total number of conformations
		assertThat(confPos.confs.length, is(numConfs));
	}
}
