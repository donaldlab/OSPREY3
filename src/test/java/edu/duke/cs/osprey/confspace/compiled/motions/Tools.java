package edu.duke.cs.osprey.confspace.compiled.motions;

import edu.duke.cs.osprey.confspace.compiled.AssignedCoords;
import edu.duke.cs.osprey.confspace.compiled.CoordsList;
import org.joml.Vector3d;
import org.joml.Vector3dc;

import java.util.Set;
import java.util.function.BiConsumer;
import java.util.stream.Collectors;

import static edu.duke.cs.osprey.TestBase.isAbsolutely;
import static org.hamcrest.Matchers.is;
import static org.hamcrest.MatcherAssert.assertThat;


public class Tools {

	public static void forEachAtom(AssignedCoords assignedCoords, CoordsList initialCoords, BiConsumer<Vector3dc,Vector3dc> func) {
		forEachAtom(
			assignedCoords.atoms().stream()
				.map(i -> i.coordsi)
				.collect(Collectors.toSet()),
			assignedCoords,
			initialCoords,
			func
		);
	}

	public static void forEachAtom(Set<Integer> atomIndices, AssignedCoords assignedCoords, CoordsList initialCoords, BiConsumer<Vector3dc,Vector3dc> func) {

		assertThat(assignedCoords.coords.size, is(initialCoords.size));

		// call the function for every matched pair of atoms
		var pos = new Vector3d();
		var ipos = new Vector3d();
		for (var i : atomIndices) {
			if (atomIndices.contains(i)) {
				assignedCoords.coords.get(i, pos);
				initialCoords.get(i, ipos);
				func.accept(pos, ipos);
			}
		}
	}

	public static void assertPos(Vector3dc obs, Vector3dc exp, double epsilon) {
		assertThat(obs.x(), isAbsolutely(exp.x(), epsilon));
		assertThat(obs.y(), isAbsolutely(exp.y(), epsilon));
		assertThat(obs.z(), isAbsolutely(exp.z(), epsilon));
	}
}
