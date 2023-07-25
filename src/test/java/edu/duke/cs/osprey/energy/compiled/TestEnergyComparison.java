package edu.duke.cs.osprey.energy.compiled;

import static org.hamcrest.Matchers.*;
import static org.hamcrest.MatcherAssert.*;
import static edu.duke.cs.osprey.TestBase.isAbsolutely;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.compiled.ConfSpace;
import edu.duke.cs.osprey.confspace.compiled.PosInterDist;
import edu.duke.cs.osprey.confspace.compiled.TestConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.ematrix.compiled.EmatCalculator;
import edu.duke.cs.osprey.energy.*;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.kstar.pfunc.GradientDescentPfunc;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.parallelism.TaskExecutor;
import edu.duke.cs.osprey.tools.MathTools.DoubleBounds;
import org.junit.jupiter.api.Test;

import java.util.List;
import java.util.stream.Collectors;


/**
 * Compare energies computed by the old and new conf space tools
 * on conf spaces that have been prepared to be as similar to each other as possible
 */
public class TestEnergyComparison {

	// NOTE: comparisons with anything less than full conformation energies won't work,
	// since in compiled conf spaces, the design positions generally won't contain all the atoms from the residue

	// also, use rigid energies for comparisons, since that will avoid noise from minimization
	private static final boolean Minimize = false;

	// NOTE: don't use reference energies here either, since they break comparability between the old and new code
	// the atoms at each design position are totally different, and reference energies are based on the "single" energies

	private static double calcEnergy(SimpleConfSpace confSpace, int[] conf) {
		var frag = new RCTuple(conf);
		var inters = ResInterGen.of(confSpace)
			.addAll(frag)
			.make();
		var pmol = confSpace.makeMolecule(conf);
		try (var ecalc = new EnergyCalculator.Builder(new ForcefieldParams())
			.setIsMinimizing(Minimize)
			.build()) {
			return ecalc.calcEnergy(pmol, inters).energy;
		}
	}

	private static double calcEnergy(ConfSpace confSpace, int[] conf) {
		var inters = PosInterDist.all(confSpace, null, conf);
		var coords = confSpace.makeCoords(conf);
		CPUConfEnergyCalculator confEcalc = new CPUConfEnergyCalculator(confSpace);
		return confEcalc.calcEnergy(coords, inters);
	}

	private static void assertEnergy(double classic, double compiled, double delta) {
		assertThat(compiled, isAbsolutely(classic, delta));
	}


	@Test
	public void wtconf_2RL0() {

		TestConfSpace.AffinityClassic classic = TestConfSpace.Design2RL0Interface7Mut.makeClassic();
		TestConfSpace.AffinityCompiled compiled = TestConfSpace.Design2RL0Interface7Mut.makeCompiled();

		// the energies won't match exactly, there are way too many tiny differences for that
		// but the energies should be within a few kcal/mol
		assertEnergy(
			calcEnergy(classic.chainA, classic.makeConfWt(classic.chainA)),
			calcEnergy(compiled.chainA, compiled.makeConfWt(compiled.chainA)),
			3.0
		);
		assertEnergy(
			calcEnergy(classic.chainB, classic.makeConfWt(classic.chainB)),
			calcEnergy(compiled.chainB, compiled.makeConfWt(compiled.chainB)),
			3.0
		);
		assertEnergy(
			calcEnergy(classic.complex, classic.makeConfWt(classic.complex)),
			calcEnergy(compiled.complex, compiled.makeConfWt(compiled.complex)),
			9.0
		);
	}


	private static List<ConfSearch.EnergiedConf> astarConfs(SimpleConfSpace confSpace, int num) {
		try (var ecalc = new EnergyCalculator.Builder(new ForcefieldParams())
			.setIsMinimizing(Minimize)
			.build()) {

			var confEcalc = new edu.duke.cs.osprey.energy.ConfEnergyCalculator.Builder(confSpace, ecalc)
				.setEnergyPartition(EnergyPartition.Traditional)
				.setAddShellInters(true)
				.build();
			var emat = new SimplerEnergyMatrixCalculator.Builder(confEcalc)
				.setCalcConstantTerm(true)
				.build()
				.calcEnergyMatrix();

			return new ConfAStarTree.Builder(emat, confSpace)
				.setTraditional()
				.build()
				.nextConfs(num)
				.stream()
				.map(c -> confEcalc.calcEnergy(c))
				.collect(Collectors.toList());
		}
	}

	private static List<ConfSearch.EnergiedConf> astarConfs(ConfSpace confSpace, int num) {

		CPUConfEnergyCalculator confEcalc = new CPUConfEnergyCalculator(confSpace);

		var posInterDist = PosInterDist.DesmetEtAl1992;
		var emat = new EmatCalculator.Builder(confEcalc)
			.setPosInterDist(posInterDist)
			.setMinimize(Minimize)
			.setIncludeStaticStatic(true)
			.build()
			.calc();

		return new ConfAStarTree.Builder(emat, confSpace)
			.setTraditional()
			.build()
			.nextConfs(num)
			.stream()
			.map(c -> {
				var conf = c.getAssignments();
				var inters = PosInterDist.all(confSpace, conf);
				var energy = confEcalc.calcEnergy(conf, inters);
				return new ConfSearch.EnergiedConf(c, energy);
			})
			.collect(Collectors.toList());
	}

	private static void assertConfs(List<ConfSearch.EnergiedConf> classic, List<ConfSearch.EnergiedConf> compiled, double delta) {

		assertThat(compiled.size(), is(classic.size()));
		int n = classic.size();

		for (int i=0; i<n; i++) {
			assertEnergy(classic.get(i).getScore(), compiled.get(i).getScore(), delta);
			assertEnergy(classic.get(i).getEnergy(), compiled.get(i).getEnergy(), delta);
		}
	}

	@Test
	public void astar_2RL0() {

		TestConfSpace.AffinityClassic classic = TestConfSpace.Design2RL0Interface7Mut.makeClassic();
		TestConfSpace.AffinityCompiled compiled = TestConfSpace.Design2RL0Interface7Mut.makeCompiled();

		final int NumConfs = 10;

		// the energies won't match exactly, there are way too many tiny differences for that
		// but the energies should be within a few kcal/mol
		assertConfs(
			astarConfs(classic.chainA, NumConfs),
			astarConfs(compiled.chainA, NumConfs),
			3.0
		);
		assertConfs(
			astarConfs(classic.chainB, NumConfs),
			astarConfs(compiled.chainB, NumConfs),
			3.0
		);
		assertConfs(
			astarConfs(classic.complex, NumConfs),
			astarConfs(compiled.complex, NumConfs),
			9.0
		);
	}


	private static DoubleBounds pfunc(edu.duke.cs.osprey.energy.ConfEnergyCalculator confEcalc, EnergyMatrix emat, RCs rcs, double epsilon) {

		var pfunc = new GradientDescentPfunc(
			confEcalc,
			new ConfAStarTree.Builder(emat, rcs)
				.setTraditional()
				.build(),
			new ConfAStarTree.Builder(emat, rcs)
				.setTraditional()
				.build(),
			rcs.getNumConformations()
		);
		pfunc.setPreciseBcalc(true);
		try (var ctx = confEcalc.tasks.contextGroup()) {
			pfunc.setInstanceId(0);
			pfunc.putTaskContexts(ctx);
			pfunc.init(epsilon);
			pfunc.compute();
		}
		var result = pfunc.makeResult();

		assertThat(result.status, is(PartitionFunction.Status.Estimated));
		return result.values.calcFreeEnergyBoundsPrecise();
	}

	private static DoubleBounds pfunc(SimpleConfSpace confSpace, double epsilon) {
		try (var ecalc = new EnergyCalculator.Builder(new ForcefieldParams())
			.setIsMinimizing(Minimize)
			.build()) {

			var confEcalc = new edu.duke.cs.osprey.energy.ConfEnergyCalculator.Builder(confSpace, ecalc)
				.setEnergyPartition(EnergyPartition.Traditional)
				.setAddShellInters(true)
				.build();
			var emat = new SimplerEnergyMatrixCalculator.Builder(confEcalc)
				.setCalcConstantTerm(true)
				.build()
				.calcEnergyMatrix();

			var seq = confSpace.seqSpace.makeWildTypeSequence();
			var rcs = seq.makeRCs(confSpace);

			return pfunc(confEcalc, emat, rcs, epsilon);
		}
	}

	private static DoubleBounds pfunc(ConfSpace confSpace, double epsilon) {

		CPUConfEnergyCalculator confEcalc = new CPUConfEnergyCalculator(confSpace);

		var posInterDist = PosInterDist.DesmetEtAl1992;
		var emat = new EmatCalculator.Builder(confEcalc)
			.setPosInterDist(posInterDist)
			.setMinimize(Minimize)
			.setIncludeStaticStatic(true)
			.build()
			.calc();

		var seq = confSpace.seqSpace.makeWildTypeSequence();
		var rcs = seq.makeRCs(confSpace);

		try (var tasks = new TaskExecutor()) {
			var adapter = new ConfEnergyCalculatorAdapter.Builder(confEcalc, tasks)
				.setPosInterDist(posInterDist)
				.setIncludeStaticStatic(false)
				.setMinimize(Minimize)
				.setIncludeStaticStatic(true)
				.build();

			return pfunc(adapter, emat, rcs, epsilon);
		}
	}

	private static void assertPfuncs(DoubleBounds classic, DoubleBounds compiled, double delta) {
		assertEnergy(classic.lower, compiled.lower, delta);
		assertEnergy(classic.upper, compiled.upper, delta);
	}

	@Test
	public void pfunc_2RL0() {

		TestConfSpace.AffinityClassic classic = TestConfSpace.Design2RL0Interface7Mut.makeClassic();
		TestConfSpace.AffinityCompiled compiled = TestConfSpace.Design2RL0Interface7Mut.makeCompiled();

		final double Epsilon = 0.1;

		// the energies won't match exactly, there are way too many tiny differences for that
		// but the energies should be within a few kcal/mol
		assertPfuncs(
			pfunc(classic.chainA, Epsilon),
			pfunc(compiled.chainA, Epsilon),
			3.0
		);
		assertPfuncs(
			pfunc(classic.chainB, Epsilon),
			pfunc(compiled.chainB, Epsilon),
			3.0
		);
		assertPfuncs(
			pfunc(classic.complex, Epsilon),
			pfunc(compiled.complex, Epsilon),
			9.0
		);
	}
}
