package edu.duke.cs.osprey.ematrix.compiled;

import static edu.duke.cs.osprey.TestBase.isAbsoluteBound;
import static edu.duke.cs.osprey.TestBase.isAbsolutely;
import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.compiled.ConfSpace;
import edu.duke.cs.osprey.confspace.compiled.PosInterDist;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimpleReferenceEnergies;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.compiled.CPUConfEnergyCalculator;
import edu.duke.cs.osprey.energy.compiled.ConfEnergyCalculatorAdapter;
import edu.duke.cs.osprey.gmec.SimpleGMECFinder;
import edu.duke.cs.osprey.kstar.pfunc.GradientDescentPfunc;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.parallelism.TaskExecutor;
import edu.duke.cs.osprey.tools.FileTools;
import edu.duke.cs.osprey.tools.MathTools;
import org.junit.Test;


@SuppressWarnings("deprecation") // yes, we're using the deprecated adapter class: we're testing it
public class TestConfEnergyCalculatorAdapter {

	private static final ConfSpace confSpace = ConfSpace.fromBytes(FileTools.readResourceBytes("/confSpaces/dipeptide.5hydrophobic.ccsx"));

	@Test
	public void energyMatrixRigid() {

		try (TaskExecutor tasks = new TaskExecutor()) {
			CPUConfEnergyCalculator confEcalc = new CPUConfEnergyCalculator(confSpace);

			// compute the true energy matrix using the new calculator
			PosInterDist posInterDist = PosInterDist.DesmetEtAl1992;
			EnergyMatrix emat = new EmatCalculator.Builder(confEcalc)
				.setPosInterDist(posInterDist)
				.setMinimize(false)
				.build()
				.calc();

			// compare to the energy matrix computed using the adapted old calculator
			ConfEnergyCalculatorAdapter adapter = new ConfEnergyCalculatorAdapter.Builder(confEcalc, tasks)
				.setPosInterDist(posInterDist)
				.setMinimize(false)
				.build();
			EnergyMatrix adaptedEmat = new SimplerEnergyMatrixCalculator.Builder(adapter)
				.build()
				.calcEnergyMatrix();

			assertThat(adaptedEmat, is(emat));
		}
	}

	private void energyMatrixMinimized(Parallelism parallelism) {

		try (TaskExecutor tasks = parallelism.makeTaskExecutor()) {
			CPUConfEnergyCalculator confEcalc = new CPUConfEnergyCalculator(confSpace);

			// compute the true energy matrix using the new calculator
			PosInterDist posInterDist = PosInterDist.DesmetEtAl1992;
			EnergyMatrix emat = new EmatCalculator.Builder(confEcalc)
				.setPosInterDist(posInterDist)
				.setMinimize(true)
				.build()
				.calc();

			// compare to the energy matrix computed using the adapted old calculator
			ConfEnergyCalculatorAdapter adapter = new ConfEnergyCalculatorAdapter.Builder(confEcalc, tasks)
				.setPosInterDist(posInterDist)
				.setMinimize(true)
				.build();
			EnergyMatrix adaptedEmat = new SimplerEnergyMatrixCalculator.Builder(adapter)
				.build()
				.calcEnergyMatrix();

			assertThat(adaptedEmat, is(emat));
		}
	}
	@Test public void energyMatrixMinimized_CPU1() { energyMatrixMinimized(Parallelism.makeCpu(1)); }
	@Test public void energyMatrixMinimized_CPU2() { energyMatrixMinimized(Parallelism.makeCpu(2)); }
	@Test public void energyMatrixMinimized_CPU4() { energyMatrixMinimized(Parallelism.makeCpu(4)); }

	@Test
	public void referenceEnergyRigid() {

		try (TaskExecutor tasks = new TaskExecutor()) {
			CPUConfEnergyCalculator confEcalc = new CPUConfEnergyCalculator(confSpace);

			// compute the true reference energies using the new calculator
			SimpleReferenceEnergies eref = new ErefCalculator.Builder(confEcalc)
				.setMinimize(false)
				.build()
				.calc();

			// compare to the reference energies computed using the adapted old calculator
			ConfEnergyCalculatorAdapter adapter = new ConfEnergyCalculatorAdapter.Builder(confEcalc, tasks)
				.setMinimize(false)
				.build();
			SimpleReferenceEnergies adaptedEref = new SimplerEnergyMatrixCalculator.Builder(adapter)
				.build()
				.calcReferenceEnergies();

			assertThat(adaptedEref, is(eref));
		}
	}

	@Test
	public void referenceEnergyMinimized() {

		try (TaskExecutor tasks = new TaskExecutor()) {
			CPUConfEnergyCalculator confEcalc = new CPUConfEnergyCalculator(confSpace);

			// compute the true reference energies using the new calculator
			SimpleReferenceEnergies eref = new ErefCalculator.Builder(confEcalc)
				.setMinimize(true)
				.build()
				.calc();

			// compare to the reference energies computed using the adapted old calculator
			ConfEnergyCalculatorAdapter adapter = new ConfEnergyCalculatorAdapter.Builder(confEcalc, tasks)
				.setMinimize(true)
				.build();
			SimpleReferenceEnergies adaptedEref = new SimplerEnergyMatrixCalculator.Builder(adapter)
				.build()
				.calcReferenceEnergies();

			assertThat(adaptedEref, is(eref));
		}
	}

	@Test
	public void energyMatrixWithEref() {

		try (TaskExecutor tasks = new TaskExecutor()) {
			CPUConfEnergyCalculator confEcalc = new CPUConfEnergyCalculator(confSpace);

			// compute the true energy matrix using the new calculators
			SimpleReferenceEnergies eref = new ErefCalculator.Builder(confEcalc)
				.setMinimize(false)
				.build()
				.calc();
			PosInterDist posInterDist = PosInterDist.DesmetEtAl1992;
			EnergyMatrix emat = new EmatCalculator.Builder(confEcalc)
				.setPosInterDist(posInterDist)
				.setReferenceEnergies(eref)
				.setMinimize(false)
				.build()
				.calc();

			// compare to the energy matrix computed using the adapted old calculator
			ConfEnergyCalculatorAdapter adapter = new ConfEnergyCalculatorAdapter.Builder(confEcalc, tasks)
				.setPosInterDist(posInterDist)
				.setReferenceEnergies(eref)
				.setMinimize(false)
				.build();
			EnergyMatrix adaptedEmat = new SimplerEnergyMatrixCalculator.Builder(adapter)
				.build()
				.calcEnergyMatrix();

			assertThat(adaptedEmat, is(emat));
		}
	}

	private void checkAstar(EnergyMatrix emat, ConfEnergyCalculatorAdapter adapter, double[][] exp) {

		// enumerate conformations with A*, and hopefully don't crash
		// the conf space only has few hundred confs, so we can enumerate them all easily
		ConfAStarTree astar = new ConfAStarTree.Builder(emat, confSpace)
			.setTraditional()
			.build();
		int numConfs = 0;
		while (true) {
			ConfSearch.ScoredConf conf = astar.nextConf();
			if (conf == null) {
				break;
			}
			ConfSearch.EnergiedConf econf = adapter.calcEnergy(conf);

			// make sure the scores and energies match
			if (numConfs < 10) {
				//log("\t{ %10.6f, %10.6f },", econf.getScore(), econf.getEnergy());
				assertThat(econf.getScore(), isAbsolutely(exp[numConfs][0], 1e-6));
				assertThat(econf.getEnergy(), isAbsolutely(exp[numConfs][1], 1e-6));
			}

			// make sure the lower bound is valid
			double error = econf.getScore() - econf.getEnergy();
			assertThat(error, lessThanOrEqualTo(1e-12));

			numConfs++;
		}

		// make sure we got all the confs
		assertThat(numConfs, is(289));
	}

	@Test
	public void astar() {

		try (TaskExecutor tasks = new TaskExecutor()) {
			CPUConfEnergyCalculator confEcalc = new CPUConfEnergyCalculator(confSpace);

			// compute the energy matrix
			PosInterDist posInterDist = PosInterDist.DesmetEtAl1992;
			EnergyMatrix emat = new EmatCalculator.Builder(confEcalc)
				.setPosInterDist(posInterDist)
				.setMinimize(true)
				.build()
				.calc();

			double[][] exp = {
				{ -23.810803, -23.810803 },
				{ -22.992313, -22.990484 },
				{ -22.925150, -22.819235 },
				{ -22.748507, -22.745493 },
				{ -22.661132, -22.552039 },
				{ -22.316248, -22.308067 },
				{ -22.292735, -22.292735 },
				{ -22.258313, -22.136844 },
				{ -22.203470, -22.195088 },
				{ -22.198951, -22.194634 }
			};

			ConfEnergyCalculatorAdapter adapter = new ConfEnergyCalculatorAdapter.Builder(confEcalc, tasks)
				.setPosInterDist(posInterDist)
				.setMinimize(true)
				.build();
			checkAstar(emat, adapter, exp);
		}
	}

	@Test
	public void astarWithEref() {

		try (TaskExecutor tasks = new TaskExecutor()) {
			CPUConfEnergyCalculator confEcalc = new CPUConfEnergyCalculator(confSpace);

			// compute the reference energies
			SimpleReferenceEnergies eref = new ErefCalculator.Builder(confEcalc)
				.setMinimize(true)
				.build()
				.calc();

			// compute the energy matrix with reference energies
			PosInterDist posInterDist = PosInterDist.DesmetEtAl1992;
			EnergyMatrix emat = new EmatCalculator.Builder(confEcalc)
				.setPosInterDist(posInterDist)
				.setReferenceEnergies(eref)
				.setMinimize(true)
				.build()
				.calc();

			double[][] exp = {
				{   0.755802,   0.869753 },
				{   0.870800,   0.992268 },
				{   1.026572,   1.144280 },
				{   1.138343,   1.269227 },
				{   1.205257,   1.266192 },
				{   1.208800,   1.262447 },
				{   1.261006,   1.261006 },
				{   1.330516,   1.390310 },
				{   1.333593,   1.399927 },
				{   1.387427,   1.435387 }
			};

			ConfEnergyCalculatorAdapter adapter = new ConfEnergyCalculatorAdapter.Builder(confEcalc, tasks)
				.setPosInterDist(posInterDist)
				.setMinimize(true)
				.setReferenceEnergies(eref)
				.build();
			checkAstar(emat, adapter, exp);
		}
	}

	public void findGMEC(Parallelism parallelism) {

		try (TaskExecutor tasks = parallelism.makeTaskExecutor()) {
			CPUConfEnergyCalculator confEcalc = new CPUConfEnergyCalculator(confSpace);

			// compute the energy matrix, with reference energies
			ConfEnergyCalculatorAdapter adapter = new ConfEnergyCalculatorAdapter.Builder(confEcalc, tasks)
				.setMinimize(true)
				.build();
			SimpleReferenceEnergies eref = new SimplerEnergyMatrixCalculator.Builder(adapter)
				.build()
				.calcReferenceEnergies();
			adapter = new ConfEnergyCalculatorAdapter.Builder(confEcalc, tasks)
				.setReferenceEnergies(eref)
				.setMinimize(true)
				.build();
			EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(adapter)
				.build()
				.calcEnergyMatrix();

			// define the conf search function
			ConfAStarTree astar = new ConfAStarTree.Builder(emat, confSpace)
				.setTraditional()
				.build();

			// find the GMEC, hopefully without crashing
			SimpleGMECFinder gmecFinder = new SimpleGMECFinder.Builder(astar, adapter).build();
			ConfSearch.EnergiedConf gmec = gmecFinder.find();

			// make sure we got the right conformation and energy
			assertThat(gmec.getScore(), isAbsolutely(0.755802, 1e-6));
			assertThat(gmec.getEnergy(), isAbsolutely(0.869753, 1e-6));
			assertThat(gmec.getAssignments(), is(new int[] { 12, 0 }));
		}
	}
	@Test public void findGMEC_CPU1() { findGMEC(Parallelism.makeCpu(1)); }
	@Test public void findGMEC_CPU2() { findGMEC(Parallelism.makeCpu(2)); }
	@Test public void findGMEC_CPU4() { findGMEC(Parallelism.makeCpu(4)); }

	public void calcPfunc(Parallelism parallelism) {

		try (TaskExecutor tasks = parallelism.makeTaskExecutor()) {
			CPUConfEnergyCalculator confEcalc = new CPUConfEnergyCalculator(confSpace);

			// compute the energy matrix
			ConfEnergyCalculatorAdapter adapter = new ConfEnergyCalculatorAdapter.Builder(confEcalc, tasks)
				.setMinimize(true)
				.build();
			EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(adapter)
				.build()
				.calcEnergyMatrix();

			try (TaskExecutor.ContextGroup contexts = tasks.contextGroup()) {

				// compute the partition function
				RCs rcs = new RCs(confSpace);
				GradientDescentPfunc pfunc = new GradientDescentPfunc(
					adapter,
					new ConfAStarTree.Builder(emat, rcs)
						.setTraditional()
						.build(),
					new ConfAStarTree.Builder(emat, rcs)
						.setTraditional()
						.build(),
					rcs.getNumConformations()
				);
				pfunc.setInstanceId(0);
				pfunc.putTaskContexts(contexts);
				pfunc.init(0.68);
				pfunc.compute();
				PartitionFunction.Result result = pfunc.makeResult();

				// make sure we got the right answer
				assertThat(result.status, is(PartitionFunction.Status.Estimated));
				assertThat(result.values.calcFreeEnergyBounds(), isAbsoluteBound(new MathTools.DoubleBounds(
					-24.381189,
					-24.375781
				), 1e-6));
			}
		}
	}
	@Test public void calcPfunc_CPU1() { calcPfunc(Parallelism.makeCpu(1)); }
	@Test public void calcPfunc_CPU2() { calcPfunc(Parallelism.makeCpu(2)); }
	@Test public void calcPfunc_CPU4() { calcPfunc(Parallelism.makeCpu(4)); }
}
