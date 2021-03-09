package edu.duke.cs.osprey.energy.compiled;

import static edu.duke.cs.osprey.TestBase.*;
import static edu.duke.cs.osprey.tools.Log.log;
import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;
import static org.junit.Assume.assumeTrue;

import edu.duke.cs.osprey.confspace.Conf;
import edu.duke.cs.osprey.confspace.compiled.*;
import edu.duke.cs.osprey.gpu.Structs;
import edu.duke.cs.osprey.tools.MathTools;
import org.joml.Vector3d;
import org.junit.Test;

import java.util.Arrays;
import java.util.List;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;


/**
 * NOTE: running these tests directly from your IDE may require extra JVM flags:
 * --add-modules=jdk.incubator.foreign
 */
public class TestNativeConfEnergyCalculator {

	private static final ConfSpace confSpace_2RL0 = TestConfSpace.Design2RL0Interface7Mut.makeCompiled().complex;
	private static final int[][] confs_2RL0 = {
		// full conformations
		{ 0, 0, 0, 0, 0, 0, 0 },
		{ 1, 0, 0, 0, 0, 0, 0 },
		{ 0, 1, 0, 0, 0, 0, 0 },
		{ 0, 0, 1, 0, 0, 0, 0 },
		{ 0, 0, 0, 1, 0, 0, 0 },
		{ 0, 0, 0, 0, 1, 0, 0 },
		{ 0, 0, 0, 0, 0, 1, 0 },
		{ 0, 0, 0, 0, 0, 0, 1 },
		// no assignments
		{ -1, -1, -1, -1, -1, -1, -1 },
		// partial assignments
		{ 0, -1, -1, -1, -1, -1, -1 },
		{ 0, 0, -1, -1, -1, -1, -1 },
		{ 0, 0, 0, -1, -1, -1, -1 },
		{ 0, 0, 0, 0, -1, -1, -1 },
		{ 0, 0, 0, 0, 0, -1, -1 },
		{ 0, 0, 0, 0, 0, 0, -1 }
	};
	private static final double[] calcEnergy_all_2RL0 = {
		2199.44093411,
		2205.48998686,
		2607.45981769,
		2307.90672767,
		749133.92904943,
		2241.54003600,
		2179.54796288,
		2171.14773794,
		-1456.82332401,
		-1462.23330457,
		-1436.56711552,
		-1444.30456665,
		2166.86260702,
		2162.91402230,
		2179.35366839
	};
	private static final double[] minimize_all_2RL0 = {
		-1359.27208010,
		-1357.74512549,
		-1110.74689221,
		-1328.74084045,
		14471.08373607,
		-1348.55609879,
		-1364.70178141,
		-1360.70959143,
		-1456.82332401,
		-1462.23330457,
		-1447.62657956,
		-1455.35974041,
		-1344.76116367,
		-1348.71185049,
		-1352.20687426
	};

	private static final ConfSpace confSpace_1DG9_6f = TestConfSpace.DesignSmallMolAffinity6f.makeCompiled().complex;
	private static final int[][] confs_1DG9_6f = {
		// pick a few arbitrary full conformations
		{ 5, 7, 34, 8, 8 }, // wild-type
		{ 4, 7, 34, 8, 8 },
		{ 5, 6, 34, 8, 8 },
		{ 5, 7,  5, 8, 8 },
		{ 5, 7, 34, 7, 8 },
		{ 5, 7, 34, 8, 7 },
		{ 4, 6,  5, 7, 7 }
	};
	private static final double[] calcEnergy_all_1DG9_6f = {
		-2410.40260979,
		-2409.43285354,
		-2401.20983442,
		-2308.15152838,
		-2410.56187477,
		-2405.16517485,
		-2292.86345610
	};
	private static final double[] minimize_all_1DG9_6f = {
		-2413.13885410,
		-2412.70804489,
		-2405.26589687,
		-2399.24997771,
		-2413.15282630,
		-2409.89634506,
		-2386.73740292
	};

	public static void assertCoords(int[] conf, AssignedCoords exp, CoordsList obs, double maxDist) {

		// diff the two coord lists
		int[] indices = IntStream.range(0, exp.confSpace.maxNumConfAtoms)
			.filter(i -> {
				Vector3d expa = new Vector3d();
				Vector3d obsa = new Vector3d();
				exp.coords.get(i, expa);
				obs.get(i, obsa);
				return expa.distance(obsa) > maxDist;
			})
			.toArray();
		List<String> diffs = IntStream.range(0, indices.length)
			.mapToObj(i -> {
				int prevAtomi = i > 0 ? indices[i - 1] : -1;
				int atomi = indices[i];
				return String.format("%s%4d   EXP: %10.6f, %10.6f, %10.6f    OBS: %10.6f, %10.6f, %10.6f   %s",
					atomi - prevAtomi > 1 ? "...\n" : "",
					atomi,
					exp.coords.x(atomi), exp.coords.y(atomi), exp.coords.z(atomi),
					obs.x(atomi), obs.y(atomi), obs.z(atomi),
					exp.getAtomName(i)
				);
			})
			.collect(Collectors.toList());
		String diffMsg = String.format("Coords are different at %d/%d positions:\n%s\nfor conformation: %s",
			diffs.size(),
			exp.confSpace.maxNumConfAtoms,
			String.join("\n", diffs),
			Conf.toString(conf)
		);

		assertThat(diffMsg, diffs.size(), is(0));
	}

	public static void assertCoords(int[] conf, AssignedCoords exp, CoordsList obs) {
		assertCoords(conf, exp, obs, 0.0);
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

			var exp = confEcalc.confSpace().makeCoords(conf);
			setPrecision(exp.coords, confEcalc.precision());

			var obs = f.apply(conf);
			assertCoords(conf, exp, obs);
		}
	}

	private void assign_native(ConfSpace confSpace, int[][] confs, Structs.Precision precision) {
		try (var confEcalc = new NativeConfEnergyCalculator(confSpace, precision)) {
			assign(confEcalc, confs, conf -> confEcalc.assign(conf).coords);
		}
	}
	@Test public void assign_native_2RL0_f32() { assign_native(confSpace_2RL0, confs_2RL0, Structs.Precision.Float32); }
	@Test public void assign_native_2RL0_f64() { assign_native(confSpace_2RL0, confs_2RL0, Structs.Precision.Float64); }
	@Test public void assign_native_1DG9_6f_f32() { assign_native(confSpace_1DG9_6f, confs_1DG9_6f, Structs.Precision.Float32); }
	@Test public void assign_native_1DG9_6f_f64() { assign_native(confSpace_1DG9_6f, confs_1DG9_6f, Structs.Precision.Float64); }

	private void assign_cuda(ConfSpace confSpace, int[][] confs, Structs.Precision precision) {
		assumeTrue(CudaConfEnergyCalculator.isSupported());
		try (var confEcalc = new CudaConfEnergyCalculator(confSpace, precision)) {
			assign(confEcalc, confs, conf -> confEcalc.assign(conf).coords);
		}
	}
	@Test public void assign_cuda_2RL0_f32() { assign_cuda(confSpace_2RL0, confs_2RL0, Structs.Precision.Float32); }
	@Test public void assign_cuda_2RL0_f64() { assign_cuda(confSpace_2RL0, confs_2RL0, Structs.Precision.Float64); }
	@Test public void assign_cuda_1DG9_6f_f32() { assign_cuda(confSpace_1DG9_6f, confs_1DG9_6f, Structs.Precision.Float32); }
	@Test public void assign_cuda_1DG9_6f_f64() { assign_cuda(confSpace_1DG9_6f, confs_1DG9_6f, Structs.Precision.Float64); }

	private void calcEnergy_all(ConfEnergyCalculator confEcalc, int[][] confs, double[] energies, double epsilon) {

		assertThat(energies.length, is(confs.length));

		for (int i=0; i<confs.length; i++) {
			var inters = PosInterDist.all(confEcalc.confSpace(), confs[i]);
			double energy = confEcalc.calcEnergy(confs[i], inters);
			assertThat("conf " + i, energy, isRelatively(energies[i], epsilon));
		}
	}

	private void calcEnergy_native_all(ConfSpace confSpace, int[][] confs, double[] energies, Structs.Precision precision, double epsilon) {
		try (var confEcalc = new NativeConfEnergyCalculator(confSpace, precision)) {
			calcEnergy_all(confEcalc, confs, energies, epsilon);
		}
	}
	@Test public void calcEnergy_native_all_2RL0_f32() { calcEnergy_native_all(confSpace_2RL0, confs_2RL0, calcEnergy_all_2RL0, Structs.Precision.Float32, 1e-5); }
	@Test public void calcEnergy_native_all_2RL0_f64() { calcEnergy_native_all(confSpace_2RL0, confs_2RL0, calcEnergy_all_2RL0, Structs.Precision.Float64, 1e-8); }
	@Test public void calcEnergy_native_all_1DG9_6f_f32() { calcEnergy_native_all(confSpace_1DG9_6f, confs_1DG9_6f, calcEnergy_all_1DG9_6f, Structs.Precision.Float32, 1e-5); }
	@Test public void calcEnergy_native_all_1DG9_6f_f64() { calcEnergy_native_all(confSpace_1DG9_6f, confs_1DG9_6f, calcEnergy_all_1DG9_6f, Structs.Precision.Float64, 1e-8); }

	private void calcEnergy_cuda_all(ConfSpace confSpace, int[][] confs, double[] energies, Structs.Precision precision, double epsilon) {
		assumeTrue(CudaConfEnergyCalculator.isSupported());
		try (var confEcalc = new CudaConfEnergyCalculator(confSpace, precision)) {
			calcEnergy_all(confEcalc, confs, energies, epsilon);
		}
	}
	@Test public void calcEnergy_cuda_all_2RL0_f32() { calcEnergy_cuda_all(confSpace_2RL0, confs_2RL0, calcEnergy_all_2RL0, Structs.Precision.Float32, 1e-5); }
	@Test public void calcEnergy_cuda_all_2RL0_f64() { calcEnergy_cuda_all(confSpace_2RL0, confs_2RL0, calcEnergy_all_2RL0, Structs.Precision.Float64, 1e-6); }
	@Test public void calcEnergy_cuda_all_1DG9_6f_f32() { calcEnergy_cuda_all(confSpace_1DG9_6f, confs_1DG9_6f, calcEnergy_all_1DG9_6f, Structs.Precision.Float32, 1e-5); }
	@Test public void calcEnergy_cuda_all_1DG9_6f_f64() { calcEnergy_cuda_all(confSpace_1DG9_6f, confs_1DG9_6f, calcEnergy_all_1DG9_6f, Structs.Precision.Float64, 1e-6); }

	private void calc_all(ConfEnergyCalculator confEcalc, int[][] confs, double[] energies, double epsilon) {

		assertThat(energies.length, is(confs.length));

		for (int i=0; i<confs.length; i++) {
			var inters = PosInterDist.all(confEcalc.confSpace(), confs[i]);
			ConfEnergyCalculator.EnergiedCoords energiedCoords = confEcalc.calc(confs[i], inters);
			assertThat("conf " + i, energiedCoords.energy, isRelatively(energies[i], epsilon));

			var expCoords = confEcalc.confSpace().makeCoords(confs[i]);
			setPrecision(expCoords.coords, confEcalc.precision());
			assertCoords(confs[i], expCoords, energiedCoords.coords.coords);
		}
	}

	private void calc_native_all(ConfSpace confSpace, int[][] confs, double[] energies, Structs.Precision precision, double epsilon) {
		try (var confEcalc = new NativeConfEnergyCalculator(confSpace, precision)) {
			calc_all(confEcalc, confs, energies, epsilon);
		}
	}
	@Test public void calc_native_all_2RL0_f32() { calc_native_all(confSpace_2RL0, confs_2RL0, calcEnergy_all_2RL0, Structs.Precision.Float32, 1e-5); }
	@Test public void calc_native_all_2RL0_f64() { calc_native_all(confSpace_2RL0, confs_2RL0, calcEnergy_all_2RL0, Structs.Precision.Float64, 1e-8); }
	@Test public void calc_native_all_1DG9_6f_f32() { calc_native_all(confSpace_1DG9_6f, confs_1DG9_6f, calcEnergy_all_1DG9_6f, Structs.Precision.Float32, 1e-5); }
	@Test public void calc_native_all_1DG9_6f_f64() { calc_native_all(confSpace_1DG9_6f, confs_1DG9_6f, calcEnergy_all_1DG9_6f, Structs.Precision.Float64, 1e-8); }

	private void calc_cuda_all(ConfSpace confSpace, int[][] confs, double[] energies, Structs.Precision precision, double epsilon) {
		assumeTrue(CudaConfEnergyCalculator.isSupported());
		try (var confEcalc = new CudaConfEnergyCalculator(confSpace, precision)) {
			calc_all(confEcalc, confs, energies, epsilon);
		}
	}
	@Test public void calc_cuda_all_2RL0_f32() { calc_cuda_all(confSpace_2RL0, confs_2RL0, calcEnergy_all_2RL0, Structs.Precision.Float32, 1e-5); }
	@Test public void calc_cuda_all_2RL0_f64() { calc_cuda_all(confSpace_2RL0, confs_2RL0, calcEnergy_all_2RL0, Structs.Precision.Float64, 1e-6); }
	@Test public void calc_cuda_all_1DG9_6f_f32() { calc_cuda_all(confSpace_1DG9_6f, confs_1DG9_6f, calcEnergy_all_1DG9_6f, Structs.Precision.Float32, 1e-5); }
	@Test public void calc_cuda_all_1DG9_6f_f64() { calc_cuda_all(confSpace_1DG9_6f, confs_1DG9_6f, calcEnergy_all_1DG9_6f, Structs.Precision.Float64, 1e-6); }


	private void minimizeEnergy_all(ConfEnergyCalculator confEcalc, int[][] confs, double[] energies, double epsilon) {

		assertThat(energies.length, is(confs.length));

		for (int i=0; i<confs.length; i++) {
			var inters = PosInterDist.all(confEcalc.confSpace(), confs[i]);
			double energy = confEcalc.minimizeEnergy(confs[i], inters);
			assertThat("conf " + i, energy, isRelatively(energies[i], epsilon));
		}
	}

	private void minimizeEnergy_native_all(ConfSpace confSpace, int[][] confs, double[] energies, Structs.Precision precision, double epsilon) {
		try (var confEcalc = new NativeConfEnergyCalculator(confSpace, precision)) {
			minimizeEnergy_all(confEcalc, confs, energies, epsilon);
		}
	}
	@Test public void minimizeEnergy_native_2RL0_f32() { minimizeEnergy_native_all(confSpace_2RL0, confs_2RL0, minimize_all_2RL0, Structs.Precision.Float32, 1e-4); }
	@Test public void minimizeEnergy_native_2RL0_f64() { minimizeEnergy_native_all(confSpace_2RL0, confs_2RL0, minimize_all_2RL0, Structs.Precision.Float64, 1e-8); }
	@Test public void minimizeEnergy_native_1DG9_6f_f32() { minimizeEnergy_native_all(confSpace_1DG9_6f, confs_1DG9_6f, minimize_all_1DG9_6f, Structs.Precision.Float32, 1e-3); }
	@Test public void minimizeEnergy_native_1DG9_6f_f64() { minimizeEnergy_native_all(confSpace_1DG9_6f, confs_1DG9_6f, minimize_all_1DG9_6f, Structs.Precision.Float64, 1e-7); }

	private void minimizeEnergy_cuda_all(ConfSpace confSpace, int[][] confs, double[] energies, Structs.Precision precision, double epsilon) {
		assumeTrue(CudaConfEnergyCalculator.isSupported());
		try (var confEcalc = new CudaConfEnergyCalculator(confSpace, precision)) {
			minimizeEnergy_all(confEcalc, confs, energies, epsilon);
		}
	}
	@Test public void minimizeEnergy_cuda_2RL0_f32() { minimizeEnergy_cuda_all(confSpace_2RL0, confs_2RL0, minimize_all_2RL0, Structs.Precision.Float32, 1e-4); }
	@Test public void minimizeEnergy_cuda_2RL0_f64() { minimizeEnergy_cuda_all(confSpace_2RL0, confs_2RL0, minimize_all_2RL0, Structs.Precision.Float64, 1e-7); }
	@Test public void minimizeEnergy_cuda_1DG9_6f_f32() { minimizeEnergy_cuda_all(confSpace_1DG9_6f, confs_1DG9_6f, minimize_all_1DG9_6f, Structs.Precision.Float32, 1e-3); }
	@Test public void minimizeEnergy_cuda_1DG9_6f_f64() { minimizeEnergy_cuda_all(confSpace_1DG9_6f, confs_1DG9_6f, minimize_all_1DG9_6f, Structs.Precision.Float64, 1e-4); }


	private void minimize_all(ConfEnergyCalculator confEcalc, int[][] confs, double[] energies, double epsilon) {

		assertThat(energies.length, is(confs.length));

		for (int i=0; i<confs.length; i++) {
			int[] conf = confs[i];

			AssignedCoords coords = confEcalc.confSpace().makeCoords(conf);

			var inters = PosInterDist.all(confEcalc.confSpace(), conf);
			var energiedCoords = confEcalc.minimize(conf, inters);

			// check the energy, obviously
			assertThat("conf " + i, energiedCoords.energy, isRelatively(energies[i], epsilon));

			// make sure we have coords that look reasonable
			assertCoords(conf, coords, energiedCoords.coords.coords, 10.0);

			// make sure the final dof values are still in range
			assertThat(energiedCoords.dofValues.size(), is(coords.dofs.size()));
			for (int d=0; d<coords.dofs.size(); d++) {
				assertThat(String.format("dof %d", d), energiedCoords.dofValues.get(d), isAbsolutelyBounded(new MathTools.DoubleBounds(
					coords.dofs.get(d).min(),
					coords.dofs.get(d).max()
				), 1e-4));
			}
		}
	}

	private void minimize_native_all(ConfSpace confSpace, int[][] confs, double[] energies, Structs.Precision precision, double epsilon) {
		try (var confEcalc = new NativeConfEnergyCalculator(confSpace, precision)) {
			minimize_all(confEcalc, confs, energies, epsilon);
		}
	}
	@Test public void minimize_native_2RL0_f32() { minimize_native_all(confSpace_2RL0, confs_2RL0, minimize_all_2RL0, Structs.Precision.Float32, 1e-4); }
	@Test public void minimize_native_2RL0_f64() { minimize_native_all(confSpace_2RL0, confs_2RL0, minimize_all_2RL0, Structs.Precision.Float64, 1e-8); }
	@Test public void minimize_native_1DG9_6f_f32() { minimize_native_all(confSpace_1DG9_6f, confs_1DG9_6f, minimize_all_1DG9_6f, Structs.Precision.Float32, 1e-3); }
	@Test public void minimize_native_1DG9_6f_f64() { minimize_native_all(confSpace_1DG9_6f, confs_1DG9_6f, minimize_all_1DG9_6f, Structs.Precision.Float64, 1e-7); }

	private void minimize_cuda_all(ConfSpace confSpace, int[][] confs, double[] energies, Structs.Precision precision, double epsilon) {
		assumeTrue(CudaConfEnergyCalculator.isSupported());
		try (var confEcalc = new CudaConfEnergyCalculator(confSpace, precision)) {
			minimize_all(confEcalc, confs, energies, epsilon);
		}
	}
	@Test public void minimize_cuda_2RL0_f32() { minimize_cuda_all(confSpace_2RL0, confs_2RL0, minimize_all_2RL0, Structs.Precision.Float32, 1e-4); }
	@Test public void minimize_cuda_2RL0_f64() { minimize_cuda_all(confSpace_2RL0, confs_2RL0, minimize_all_2RL0, Structs.Precision.Float64, 1e-7); }
	@Test public void minimize_cuda_1DG9_6f_f32() { minimize_cuda_all(confSpace_1DG9_6f, confs_1DG9_6f, minimize_all_1DG9_6f, Structs.Precision.Float32, 1e-3); }
	@Test public void minimize_cuda_1DG9_6f_f64() { minimize_cuda_all(confSpace_1DG9_6f, confs_1DG9_6f, minimize_all_1DG9_6f, Structs.Precision.Float64, 1e-4); }


	private void minimizeEnergies_all(ConfEnergyCalculator confEcalc, int[][] confs, double[] energies, double epsilon) {

		assertThat(energies.length, is(confs.length));

		// make the minimization jobs
		List<ConfEnergyCalculator.MinimizationJob> jobs = Arrays.stream(confs)
			.map(conf -> {
				var inters = PosInterDist.all(confEcalc.confSpace(), conf);
				return new ConfEnergyCalculator.MinimizationJob(conf, inters);
			})
			.collect(Collectors.toList());

		confEcalc.minimizeEnergies(jobs);
		for (int i=0; i<confs.length; i++) {
			assertThat("conf " + i, jobs.get(i).energy, isRelatively(energies[i], epsilon));
		}
	}

	private void minimizeEnergies_native_all(ConfSpace confSpace, int[][] confs, double[] energies, Structs.Precision precision, double epsilon) {
		try (var confEcalc = new NativeConfEnergyCalculator(confSpace, precision)) {
			minimizeEnergies_all(confEcalc, confs, energies, epsilon);
		}
	}
	@Test public void minimizeEnergies_native_2RL0_f32() { minimizeEnergies_native_all(confSpace_2RL0, confs_2RL0, minimize_all_2RL0, Structs.Precision.Float32, 1e-4); }
	@Test public void minimizeEnergies_native_2RL0_f64() { minimizeEnergies_native_all(confSpace_2RL0, confs_2RL0, minimize_all_2RL0, Structs.Precision.Float64, 1e-8); }
	@Test public void minimizeEnergies_native_1DG9_6f_f32() { minimizeEnergies_native_all(confSpace_1DG9_6f, confs_1DG9_6f, minimize_all_1DG9_6f, Structs.Precision.Float32, 1e-3); }
	@Test public void minimizeEnergies_native_1DG9_6f_f64() { minimizeEnergies_native_all(confSpace_1DG9_6f, confs_1DG9_6f, minimize_all_1DG9_6f, Structs.Precision.Float64, 1e-7); }

	private void minimizeEnergies_cuda_all(ConfSpace confSpace, int[][] confs, double[] energies, Structs.Precision precision, double epsilon) {
		assumeTrue(CudaConfEnergyCalculator.isSupported());
		try (var confEcalc = new CudaConfEnergyCalculator(confSpace, precision)) {
			minimizeEnergies_all(confEcalc, confs, energies, epsilon);
		}
	}
	@Test public void minimizeEnergies_cuda_2RL0_f32() { minimizeEnergies_cuda_all(confSpace_2RL0, confs_2RL0, minimize_all_2RL0, Structs.Precision.Float32, 1e-2); }
	@Test public void minimizeEnergies_cuda_2RL0_f64() { minimizeEnergies_cuda_all(confSpace_2RL0, confs_2RL0, minimize_all_2RL0, Structs.Precision.Float64, 1e-2); }
	@Test public void minimizeEnergies_cuda_1DG9_6f_f32() { minimizeEnergies_cuda_all(confSpace_1DG9_6f, confs_1DG9_6f, minimize_all_1DG9_6f, Structs.Precision.Float32, 1e-3); }
	@Test public void minimizeEnergies_cuda_1DG9_6f_f64() { minimizeEnergies_cuda_all(confSpace_1DG9_6f, confs_1DG9_6f, minimize_all_1DG9_6f, Structs.Precision.Float64, 1e-4); }


	public static void main(String[] args) {

		// generate the expected values
		//dumpEnergiesAll(confSpace_2RL0, confs_2RL0);
		dumpEnergiesAll(confSpace_1DG9_6f, confs_1DG9_6f);
	}

	private static void dumpEnergiesAll(ConfSpace confSpace, int[][] confs) {

		var confEcalc = new CPUConfEnergyCalculator(confSpace);

		log("calcEnergy");
		for (int[] conf : confs) {
			var inters = PosInterDist.all(confSpace, conf);
			double energy = confEcalc.calcEnergy(conf, inters);
			log("%16.8f", energy);
		}

		log("minimizeEnergy");
		for (int[] conf : confs) {
			var inters = PosInterDist.all(confSpace, conf);
			double energy = confEcalc.minimizeEnergy(conf, inters);
			log("%16.8f", energy);
		}
	}
}
