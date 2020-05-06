package edu.duke.cs.osprey.energy.compiled;

import static edu.duke.cs.osprey.TestBase.isAbsolutely;
import static edu.duke.cs.osprey.TestBase.isRelatively;
import static edu.duke.cs.osprey.tools.Log.log;
import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import edu.duke.cs.osprey.confspace.Conf;
import edu.duke.cs.osprey.confspace.compiled.*;
import edu.duke.cs.osprey.gpu.Structs;
import org.junit.Test;

import java.util.List;
import java.util.function.Function;
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
	private static final double[] calcEnergy_all_2RL0 = {
		2199.44093411,
		2205.48998686,
		2607.45981769,
		2307.90672767,
		749133.92904943,
		2241.54003600,
		2179.54796288,
		2171.14773794
	};
	private static final double[] minimize_all_2RL0 = {
		-1359.27208010,
		-1357.74512549,
		-1110.74689221,
		-1328.74084045,
		14471.08373607,
		-1348.55609879,
		-1364.70178141,
		-1360.70959143
	};

	public static void assertCoords(int[] conf, ConfSpace confSpace, CoordsList exp, CoordsList obs) {

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
		String diffMsg = String.format("Coords are different at %d/%d positions:\n%s\nfor conformation: %s",
			diffs.size(),
			confSpace.maxNumConfAtoms,
			String.join("\n", diffs),
			Conf.toString(conf)
		);

		assertThat(diffMsg, diffs.size(), is(0));
	}

	public static void setPrecision(CoordsList coords, Structs.Precision precision) {
		for (int i=0; i<coords.size; i++) {
			coords.setX(i, precision.cast(coords.x(i)));
			coords.setY(i, precision.cast(coords.y(i)));
			coords.setZ(i, precision.cast(coords.z(i)));
		}
	}

	private void assign(ConfEnergyCalculator confEcalc, int[][] confs, Function<int[],CoordsList> f) {
		for (int[] conf : confs) {

			var exp = confEcalc.confSpace().makeCoords(conf).coords;
			setPrecision(exp, confEcalc.precision());

			var obs = f.apply(conf);
			assertCoords(conf, confEcalc.confSpace(), exp, obs);
		}
	}

	private void assign_native(ConfSpace confSpace, int[][] confs, Structs.Precision precision) {
		try (var confEcalc = new NativeConfEnergyCalculator(confSpace, precision)) {
			assign(confEcalc, confs, conf -> confEcalc.assign(conf).coords);
		}
	}
	@Test public void assign_native_2RL0_f32() { assign_native(confSpace_2RL0, confs_2RL0, Structs.Precision.Float32); }
	@Test public void assign_native_2RL0_f64() { assign_native(confSpace_2RL0, confs_2RL0, Structs.Precision.Float64); }

	private void assign_cuda(ConfSpace confSpace, int[][] confs, Structs.Precision precision) {
		try (var confEcalc = new CudaConfEnergyCalculator(confSpace, precision)) {
			assign(confEcalc, confs, conf -> confEcalc.assign(conf).coords);
		}
	}
	@Test public void assign_cuda_2RL0_f32() { assign_cuda(confSpace_2RL0, confs_2RL0, Structs.Precision.Float32); }
	@Test public void assign_cuda_2RL0_f64() { assign_cuda(confSpace_2RL0, confs_2RL0, Structs.Precision.Float64); }

	private void calcEnergy_all(ConfEnergyCalculator confEcalc, int[][] confs, double[] energies, double epsilon) {

		// make complete position interactions
		List<PosInter> inters = PosInterDist.all(confEcalc.confSpace());

		assertThat(energies.length, is(confs.length));

		for (int i=0; i<confs.length; i++) {
			double energy = confEcalc.calcEnergy(confs[i], inters);
			switch (confEcalc.precision()) {
				case Float32 -> assertThat("conf " + i, energy, isRelatively(energies[i], epsilon));
				case Float64 -> assertThat("conf " + i, energy, isAbsolutely(energies[i], epsilon));
			}
		}
	}

	private void calcEnergy_native_all(ConfSpace confSpace, int[][] confs, double[] energies, Structs.Precision precision, double epsilon) {
		try (var confEcalc = new NativeConfEnergyCalculator(confSpace, precision)) {
			calcEnergy_all(confEcalc, confs, energies, epsilon);
		}
	}
	@Test public void calcEnergy_native_all_2RL0_f32() { calcEnergy_native_all(confSpace_2RL0, confs_2RL0, calcEnergy_all_2RL0, Structs.Precision.Float32, 1e-5); }
	@Test public void calcEnergy_native_all_2RL0_f64() { calcEnergy_native_all(confSpace_2RL0, confs_2RL0, calcEnergy_all_2RL0, Structs.Precision.Float64, 1e-8); }

	private void calcEnergy_cuda_all(ConfSpace confSpace, int[][] confs, double[] energies, Structs.Precision precision, double epsilon) {
		try (var confEcalc = new CudaConfEnergyCalculator(confSpace, precision)) {
			calcEnergy_all(confEcalc, confs, energies, epsilon);
		}
	}
	@Test public void calcEnergy_cuda_all_2RL0_f32() { calcEnergy_cuda_all(confSpace_2RL0, confs_2RL0, calcEnergy_all_2RL0, Structs.Precision.Float32, 1e-5); }
	@Test public void calcEnergy_cuda_all_2RL0_f64() { calcEnergy_cuda_all(confSpace_2RL0, confs_2RL0, calcEnergy_all_2RL0, Structs.Precision.Float64, 1e-8); }

	private void calc_all(ConfEnergyCalculator confEcalc, int[][] confs, double[] energies, double epsilon) {

		// make complete position interactions
		List<PosInter> inters = PosInterDist.all(confEcalc.confSpace());

		assertThat(energies.length, is(confs.length));

		for (int i=0; i<confs.length; i++) {
			ConfEnergyCalculator.EnergiedCoords energiedCoords = confEcalc.calc(confs[i], inters);

			switch (confEcalc.precision()) {
				case Float32 -> assertThat("conf " + i, energiedCoords.energy, isRelatively(energies[i], epsilon));
				case Float64 -> assertThat("conf " + i, energiedCoords.energy, isAbsolutely(energies[i], epsilon));
			}

			var expCoords = confEcalc.confSpace().makeCoords(confs[i]).coords;
			setPrecision(expCoords, confEcalc.precision());
			assertCoords(confs[i], confEcalc.confSpace(), expCoords, energiedCoords.coords.coords);
		}
	}

	private void calc_native_all(ConfSpace confSpace, int[][] confs, double[] energies, Structs.Precision precision, double epsilon) {
		try (var confEcalc = new NativeConfEnergyCalculator(confSpace, precision)) {
			calc_all(confEcalc, confs, energies, epsilon);
		}
	}
	@Test public void calc_native_all_2RL0_f32() { calc_native_all(confSpace_2RL0, confs_2RL0, calcEnergy_all_2RL0, Structs.Precision.Float32, 1e-5); }
	@Test public void calc_native_all_2RL0_f64() { calc_native_all(confSpace_2RL0, confs_2RL0, calcEnergy_all_2RL0, Structs.Precision.Float64, 1e-8); }

	private void calc_cuda_all(ConfSpace confSpace, int[][] confs, double[] energies, Structs.Precision precision, double epsilon) {
		try (var confEcalc = new CudaConfEnergyCalculator(confSpace, precision)) {
			calc_all(confEcalc, confs, energies, epsilon);
		}
	}
	@Test public void calc_cuda_all_2RL0_f32() { calc_cuda_all(confSpace_2RL0, confs_2RL0, calcEnergy_all_2RL0, Structs.Precision.Float32, 1e-5); }
	@Test public void calc_cuda_all_2RL0_f64() { calc_cuda_all(confSpace_2RL0, confs_2RL0, calcEnergy_all_2RL0, Structs.Precision.Float64, 1e-8); }


	private void minimizeEnergy_all(ConfEnergyCalculator confEcalc, int[][] confs, double[] energies, double epsilon) {

		// make complete position interactions
		List<PosInter> inters = PosInterDist.all(confEcalc.confSpace());

		assertThat(energies.length, is(confs.length));

		for (int i=0; i<confs.length; i++) {
			double energy = confEcalc.minimizeEnergy(confs[i], inters);
			switch (confEcalc.precision()) {
				case Float32 -> assertThat("conf " + i, energy, isRelatively(energies[i], epsilon));
				case Float64 -> assertThat("conf " + i, energy, isAbsolutely(energies[i], epsilon));
			}
		}
	}

	private void minimizeEnergy_native_all(ConfSpace confSpace, int[][] confs, double[] energies, Structs.Precision precision, double epsilon) {
		try (var confEcalc = new NativeConfEnergyCalculator(confSpace, precision)) {
			minimizeEnergy_all(confEcalc, confs, energies, epsilon);
		}
	}
	@Test public void confDihedrals_native_2RL0_f32() { minimizeEnergy_native_all(confSpace_2RL0, confs_2RL0, minimize_all_2RL0, Structs.Precision.Float32, 1e-4); }
	@Test public void confDihedrals_native_2RL0_f64() { minimizeEnergy_native_all(confSpace_2RL0, confs_2RL0, minimize_all_2RL0, Structs.Precision.Float64, 1e-8); }

	private void minimizeEnergy_cuda_all(ConfSpace confSpace, int[][] confs, double[] energies, Structs.Precision precision, double epsilon) {
		try (var confEcalc = new CudaConfEnergyCalculator(confSpace, precision)) {
			minimizeEnergy_all(confEcalc, confs, energies, epsilon);
		}
	}
	@Test public void confDihedrals_cuda_2RL0_f32() { minimizeEnergy_cuda_all(confSpace_2RL0, confs_2RL0, minimize_all_2RL0, Structs.Precision.Float32, 1e-4); }
	@Test public void confDihedrals_cuda_2RL0_f64() { minimizeEnergy_cuda_all(confSpace_2RL0, confs_2RL0, minimize_all_2RL0, Structs.Precision.Float64, 1e-8); }

	// TODO: test minimize(), ie not just minimizeEnergy()


	public static void main(String[] args) {

		// generate the expected values
		dumpEnergiesAll(confSpace_2RL0, confs_2RL0);
	}

	private static void dumpEnergiesAll(ConfSpace confSpace, int[][] confs) {
		dumpEnergies(confSpace, confs, PosInterDist.all(confSpace));
	}

	private static void dumpEnergies(ConfSpace confSpace, int[][] confs, List<PosInter> inters) {

		var confEcalc = new CPUConfEnergyCalculator(confSpace);

		log("calcEnergy");
		for (int[] conf : confs) {
			double energy = confEcalc.calcEnergy(conf, inters);
			log("%16.8f", energy);
		}

		log("minimizeEnergy");
		for (int[] conf : confs) {
			double energy = confEcalc.minimizeEnergy(conf, inters);
			log("%16.8f", energy);
		}
	}
}
