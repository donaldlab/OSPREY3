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
import edu.duke.cs.osprey.energy.compiled.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.compiled.ConfEnergyCalculatorAdapter;
import edu.duke.cs.osprey.energy.compiled.PosInterGen;
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

	private static ConfSpace confSpace = ConfSpace.fromBytes(FileTools.readResourceBytes("/confSpaces/dipeptide.5hydrophobic.ccs.xz"));

	@Test
	public void energyMatrixRigid() {

		try (CPUConfEnergyCalculator confEcalc = new CPUConfEnergyCalculator(confSpace, 1)) {

			// compute the true energy matrix using the new calculator
			PosInterGen posInterGen = new PosInterGen(PosInterDist.DesmetEtAl1992, null);
			EnergyMatrix emat = new EmatCalculator.Builder(confEcalc, posInterGen)
				.setMinimize(false)
				.build()
				.calc();

			// compare to the energy matrix computed using the adapted old calculator
			ConfEnergyCalculatorAdapter adapter = new ConfEnergyCalculatorAdapter.Builder(confEcalc)
				.setPosInterDist(posInterGen.dist)
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
			try (ConfEnergyCalculator confEcalc = ConfEnergyCalculator.build(confSpace, parallelism, tasks)) {

				// compute the true energy matrix using the new calculator
				PosInterGen posInterGen = new PosInterGen(PosInterDist.DesmetEtAl1992, null);
				EnergyMatrix emat = new EmatCalculator.Builder(confEcalc, posInterGen)
					.setMinimize(true)
					.build()
					.calc();

				// compare to the energy matrix computed using the adapted old calculator
				ConfEnergyCalculatorAdapter adapter = new ConfEnergyCalculatorAdapter.Builder(confEcalc)
					.setPosInterDist(posInterGen.dist)
					.setMinimize(true)
					.build();
				EnergyMatrix adaptedEmat = new SimplerEnergyMatrixCalculator.Builder(adapter)
					.build()
					.calcEnergyMatrix();

				assertThat(adaptedEmat, is(emat));
			}
		}
	}
	@Test public void energyMatrixMinimized_CPU1() { energyMatrixMinimized(Parallelism.makeCpu(1)); }
	@Test public void energyMatrixMinimized_CPU2() { energyMatrixMinimized(Parallelism.makeCpu(2)); }
	@Test public void energyMatrixMinimized_CPU4() { energyMatrixMinimized(Parallelism.makeCpu(4)); }

	@Test
	public void referenceEnergyRigid() {

		try (CPUConfEnergyCalculator confEcalc = new CPUConfEnergyCalculator(confSpace, 1)) {

			// compute the true reference energies using the new calculator
			SimpleReferenceEnergies eref = new ErefCalculator.Builder(confEcalc)
				.setMinimize(false)
				.build()
				.calc();

			// compare to the reference energies computed using the adapted old calculator
			ConfEnergyCalculatorAdapter adapter = new ConfEnergyCalculatorAdapter.Builder(confEcalc)
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

		try (CPUConfEnergyCalculator confEcalc = new CPUConfEnergyCalculator(confSpace, 1)) {

			// compute the true reference energies using the new calculator
			SimpleReferenceEnergies eref = new ErefCalculator.Builder(confEcalc)
				.setMinimize(true)
				.build()
				.calc();

			// compare to the reference energies computed using the adapted old calculator
			ConfEnergyCalculatorAdapter adapter = new ConfEnergyCalculatorAdapter.Builder(confEcalc)
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

		try (CPUConfEnergyCalculator confEcalc = new CPUConfEnergyCalculator(confSpace, 1)) {

			// compute the true energy matrix using the new calculators
			SimpleReferenceEnergies eref = new ErefCalculator.Builder(confEcalc)
				.setMinimize(false)
				.build()
				.calc();
			PosInterGen posInterGen = new PosInterGen(PosInterDist.DesmetEtAl1992, eref);
			EnergyMatrix emat = new EmatCalculator.Builder(confEcalc, posInterGen)
				.setMinimize(false)
				.build()
				.calc();

			// compare to the energy matrix computed using the adapted old calculator
			ConfEnergyCalculatorAdapter adapter = new ConfEnergyCalculatorAdapter.Builder(confEcalc)
				.setPosInterDist(posInterGen.dist)
				.setReferenceEnergies(eref)
				.setMinimize(false)
				.build();
			EnergyMatrix adaptedEmat = new SimplerEnergyMatrixCalculator.Builder(adapter)
				.build()
				.calcEnergyMatrix();

			assertThat(adaptedEmat, is(emat));
		}
	}

	@Test
	public void astar() {

		try (CPUConfEnergyCalculator confEcalc = new CPUConfEnergyCalculator(confSpace, 1)) {

			// compute the energy matrix
			ConfEnergyCalculatorAdapter adapter = new ConfEnergyCalculatorAdapter.Builder(confEcalc)
				.setMinimize(true)
				.build();
			EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(adapter)
				.build()
				.calcEnergyMatrix();

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

				// make sure the lower bound is valid
				double error = econf.getScore() - econf.getEnergy();
				assertThat(error, lessThanOrEqualTo(1e-2));
				// TODO: we're seeing some non-trivial error here (eg, 0.0029, 0.0017)
				//  is this just a minimizer failure... or is this a bug?

				numConfs++;
			}

			// make sure we got all the confs
			assertThat(numConfs, is(289));
		}
	}

	public void findGMEC(Parallelism parallelism) {

		try (TaskExecutor tasks = parallelism.makeTaskExecutor()) {
			try (ConfEnergyCalculator confEcalc = ConfEnergyCalculator.build(confSpace, parallelism, tasks)) {

				// compute the energy matrix, with reference energies
				ConfEnergyCalculatorAdapter adapter = new ConfEnergyCalculatorAdapter.Builder(confEcalc)
					.setMinimize(true)
					.build();
				SimpleReferenceEnergies eref = new SimplerEnergyMatrixCalculator.Builder(adapter)
					.build()
					.calcReferenceEnergies();
				adapter = new ConfEnergyCalculatorAdapter.Builder(confEcalc)
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
				assertThat(gmec.getAssignments(), is(new int[] { 3, 0 }));
				assertThat(gmec.getScore(), isAbsolutely(20.890101, 1e-6));
				assertThat(gmec.getEnergy(), isAbsolutely(20.944495, 1e-6));
			}
		}
	}
	@Test public void findGMEC_CPU1() { findGMEC(Parallelism.makeCpu(1)); }
	@Test public void findGMEC_CPU2() { findGMEC(Parallelism.makeCpu(2)); }
	@Test public void findGMEC_CPU4() { findGMEC(Parallelism.makeCpu(4)); }

	public void calcPfunc(Parallelism parallelism) {

		try (TaskExecutor tasks = parallelism.makeTaskExecutor()) {
			try (ConfEnergyCalculator confEcalc = ConfEnergyCalculator.build(confSpace, parallelism, tasks)) {

				// compute the energy matrix
				ConfEnergyCalculatorAdapter adapter = new ConfEnergyCalculatorAdapter.Builder(confEcalc)
					.setMinimize(true)
					.build();
				EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(adapter)
					.build()
					.calcEnergyMatrix();

				// compute the partition function
				GradientDescentPfunc pfunc = new GradientDescentPfunc(adapter);
				RCs rcs = new RCs(confSpace);
				pfunc.init(
					new ConfAStarTree.Builder(emat, rcs)
						.setTraditional()
						.build(),
					new ConfAStarTree.Builder(emat, rcs)
						.setTraditional()
						.build(),
					rcs.getNumConformations(),
					0.68
				);
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

	// TODO: GPU ecalcs?
}
