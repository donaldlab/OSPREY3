package edu.duke.cs.osprey.energy.compiled;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import edu.duke.cs.osprey.confspace.compiled.AssignedCoords;
import edu.duke.cs.osprey.confspace.compiled.ConfSpace;
import edu.duke.cs.osprey.confspace.compiled.TestConfSpace;
import org.junit.Test;

import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;


public class TestNativeConfEnergyCalculator {

	private void assign(NativeConfEnergyCalculator.Precision precision) {

		ConfSpace confSpace = TestConfSpace.Design2RL0Interface7Mut.makeCompiled().complex;
		var confEcalc = new NativeConfEnergyCalculator(confSpace, precision);

		assertThat(confSpace.positions.length, is(7));

		// compare the coords for a few different conformations
		int[][] confs = {
			{ 0, 0, 0, 0, 0, 0, 0 },
			{ 1, 0, 0, 0, 0, 0, 0 },
			{ 0, 1, 0, 0, 0, 0, 0 },
			{ 0, 0, 1, 0, 0, 0, 0 },
			{ 0, 0, 0, 1, 0, 0, 0 },
			{ 0, 0, 0, 0, 1, 0, 0 },
			{ 0, 0, 0, 0, 0, 1, 0 },
			{ 0, 0, 0, 0, 0, 0, 1 }
		};
		for (int[] conf : confs) {

			var exp = new AssignedCoords(confEcalc.confSpace, conf).coords;

			// convert the expected coords to the needed precision
			var coords = new AssignedCoords(confEcalc.confSpace, conf).coords;
			for (int i=0; i<coords.size; i++) {
				coords.setX(i, precision.toDouble(precision.fromDouble(coords.x(i))));
				coords.setY(i, precision.toDouble(precision.fromDouble(coords.y(i))));
				coords.setZ(i, precision.toDouble(precision.fromDouble(coords.z(i))));
			}

			var obs = confEcalc.assign(conf);

			// diff the two coord lists
			List<String> diffs = IntStream.range(0, confSpace.maxNumConfAtoms)
				.filter(i ->
					exp.x(i) != obs.x(i)
					|| exp.y(i) != obs.y(i)
					|| exp.z(i) != obs.z(i)
				)
				.mapToObj(i ->
					String.format("%4d   EXP: %10.6f, %10.6f, %10.6f    OBS: %10.6f, %10.6f, %10.6f",
						i,
						exp.x(i), exp.y(i), exp.z(i),
						obs.x(i), obs.y(i), obs.z(i)
					)
				)
				.collect(Collectors.toList());
			String diffMsg = String.format("Coords are different at %d/%d positions:\n%s",
				diffs.size(),
				confSpace.maxNumConfAtoms,
				String.join("\n", diffs)
			);

			// TODO: position coords being written too early in the list, even in f64

			assertThat(diffMsg, diffs.size(), is(0));
		}
	}
	@Test public void assign_f32() { assign(NativeConfEnergyCalculator.Precision.Single); }
	@Test public void assign_f64() { assign(NativeConfEnergyCalculator.Precision.Double); }
}
