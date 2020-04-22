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
			for (int i=0; i<exp.size; i++) {
				exp.setX(i, precision.toDouble(precision.fromDouble(exp.x(i))));
				exp.setY(i, precision.toDouble(precision.fromDouble(exp.y(i))));
				exp.setZ(i, precision.toDouble(precision.fromDouble(exp.z(i))));
			}

			var obs = confEcalc.assign(conf);

			// diff the two coord lists
			int[] indices = IntStream.range(0, confSpace.maxNumConfAtoms)
				.filter(i ->
					exp.x(i) != obs.x(i)
					|| exp.y(i) != obs.y(i)
					|| exp.z(i) != obs.z(i)
				)
				.toArray();
			List<String> diffs = IntStream.range(0, indices.length)
				.mapToObj(i -> {
					int prevAtomi = i > 0 ? indices[i - 1] : -1;
					int atomi = indices[i];
					return String.format("%s%4d   EXP: %10.6f, %10.6f, %10.6f    OBS: %10.6f, %10.6f, %10.6f",
						atomi - prevAtomi > 1 ? "...\n" : "",
						atomi,
						exp.x(atomi), exp.y(atomi), exp.z(atomi),
						obs.x(atomi), obs.y(atomi), obs.z(atomi)
					);
				})
				.collect(Collectors.toList());
			String diffMsg = String.format("Coords are different at %d/%d positions:\n%s",
				diffs.size(),
				confSpace.maxNumConfAtoms,
				String.join("\n", diffs)
			);

			assertThat(diffMsg, diffs.size(), is(0));
		}
	}
	@Test public void assign_f32() { assign(NativeConfEnergyCalculator.Precision.Single); }
	@Test public void assign_f64() { assign(NativeConfEnergyCalculator.Precision.Double); }
}
