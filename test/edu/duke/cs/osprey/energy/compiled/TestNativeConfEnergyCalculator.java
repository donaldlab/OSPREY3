package edu.duke.cs.osprey.energy.compiled;

import static edu.duke.cs.osprey.TestBase.isAbsolutely;
import static edu.duke.cs.osprey.TestBase.isRelatively;
import static edu.duke.cs.osprey.tools.Log.log;
import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import edu.duke.cs.osprey.confspace.compiled.*;
import org.junit.Test;

import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;


public class TestNativeConfEnergyCalculator {

	private static final ConfSpace confSpace_2RL0 = TestConfSpace.Design2RL0Interface7Mut.makeCompiled().complex;
	private static final int[][] confs_2RL0 = {
		{ 0, 0, 0, 0, 0, 0, 0 },
		{ 1, 0, 0, 0, 0, 0, 0 },
		{ 0, 1, 0, 0, 0, 0, 0 },
		{ 0, 0, 1, 0, 0, 0, 0 },
		{ 0, 0, 0, 1, 0, 0, 0 },
		{ 0, 0, 0, 0, 1, 0, 0 },
		{ 0, 0, 0, 0, 0, 1, 0 },
		{ 0, 0, 0, 0, 0, 0, 1 }
	};

	private void assign(ConfSpace confSpace, int[][] confs, NativeConfEnergyCalculator.Precision precision) {

		var confEcalc = new NativeConfEnergyCalculator(confSpace, precision);

		assertThat(confSpace.positions.length, is(7));

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
	@Test public void assign_2RL0_f32() { assign(confSpace_2RL0, confs_2RL0, NativeConfEnergyCalculator.Precision.Single); }
	@Test public void assign_2RL0_f64() { assign(confSpace_2RL0, confs_2RL0, NativeConfEnergyCalculator.Precision.Double); }

	private void calc_all(ConfSpace confSpace, int[][] confs, NativeConfEnergyCalculator.Precision precision) {

		var confEcalc = new NativeConfEnergyCalculator(confSpace, precision);

		// make complete position interactions
		List<PosInter> inters = PosInterDist.all(confSpace);

		switch (precision) {
			case Single -> {
				final double epsilon = 1e-5;
				assertThat(confEcalc.calcEnergy(confs[0], inters), isRelatively(   2199.44093411, epsilon));
				assertThat(confEcalc.calcEnergy(confs[1], inters), isRelatively(   2205.48998686, epsilon));
				assertThat(confEcalc.calcEnergy(confs[2], inters), isRelatively(   2607.45981769, epsilon));
				assertThat(confEcalc.calcEnergy(confs[3], inters), isRelatively(   2307.90672767, epsilon));
				assertThat(confEcalc.calcEnergy(confs[4], inters), isRelatively( 749133.92904943, epsilon));
				assertThat(confEcalc.calcEnergy(confs[5], inters), isRelatively(   2241.54003600, epsilon));
				assertThat(confEcalc.calcEnergy(confs[6], inters), isRelatively(   2179.54796288, epsilon));
				assertThat(confEcalc.calcEnergy(confs[7], inters), isRelatively(   2171.14773794, epsilon));
			}
			case Double -> {
				final double epsilon = 1e-8;
				assertThat(confEcalc.calcEnergy(confs[0], inters), isAbsolutely(   2199.44093411, epsilon));
				assertThat(confEcalc.calcEnergy(confs[1], inters), isAbsolutely(   2205.48998686, epsilon));
				assertThat(confEcalc.calcEnergy(confs[2], inters), isAbsolutely(   2607.45981769, epsilon));
				assertThat(confEcalc.calcEnergy(confs[3], inters), isAbsolutely(   2307.90672767, epsilon));
				assertThat(confEcalc.calcEnergy(confs[4], inters), isAbsolutely( 749133.92904943, epsilon));
				assertThat(confEcalc.calcEnergy(confs[5], inters), isAbsolutely(   2241.54003600, epsilon));
				assertThat(confEcalc.calcEnergy(confs[6], inters), isAbsolutely(   2179.54796288, epsilon));
				assertThat(confEcalc.calcEnergy(confs[7], inters), isAbsolutely(   2171.14773794, epsilon));
			}
		}
	}
	@Test public void calc_complete_2RL0_f32() { calc_all(confSpace_2RL0, confs_2RL0, NativeConfEnergyCalculator.Precision.Single); }
	@Test public void calc_complete_2RL0_f64() { calc_all(confSpace_2RL0, confs_2RL0, NativeConfEnergyCalculator.Precision.Double); }


	public static void main(String[] args) {

		// generate the expected values
		dumpEnergiesAll(confSpace_2RL0, confs_2RL0);
	}

	private static void dumpEnergiesAll(ConfSpace confSpace, int[][] confs) {
		dumpEnergies(confSpace, confs, PosInterDist.all(confSpace));
	}

	private static void dumpEnergies(ConfSpace confSpace, int[][] confs, List<PosInter> inters) {

		var confEcalc = new CPUConfEnergyCalculator(confSpace);

		for (int i=0; i<confs.length; i++) {
			double energy = confEcalc.calcEnergy(confs[i], inters);
			log("assertThat(confEcalc.calcEnergy(confs[%d], inters), isAbsolutely(%16.8f, epsilon));", i, energy);
		}
	}
}
