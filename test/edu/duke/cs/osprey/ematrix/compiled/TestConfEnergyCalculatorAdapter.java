package edu.duke.cs.osprey.ematrix.compiled;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
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
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.tools.FileTools;
import org.junit.Test;


@SuppressWarnings("deprecation") // yes, we're using the deprecated adapter class: we're testing it
public class TestConfEnergyCalculatorAdapter {

	private static ConfSpace confSpace = new ConfSpace(FileTools.readResourceBytes("/confSpaces/dipeptide.5hydrophobic.ccs.toml.xz"));

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
				.setMinimizing(false)
				.build();
			EnergyMatrix adaptedEmat = new SimplerEnergyMatrixCalculator.Builder(adapter)
				.build()
				.calcEnergyMatrix();

			assertThat(adaptedEmat, is(emat));
		}
	}

	private void energyMatrixMinimized(Parallelism parallelism) {

		try (ConfEnergyCalculator confEcalc = ConfEnergyCalculator.build(confSpace, parallelism)) {

			// compute the true energy matrix using the new calculator
			PosInterGen posInterGen = new PosInterGen(PosInterDist.DesmetEtAl1992, null);
			EnergyMatrix emat = new EmatCalculator.Builder(confEcalc, posInterGen)
				.setMinimize(true)
				.build()
				.calc();

			// compare to the energy matrix computed using the adapted old calculator
			ConfEnergyCalculatorAdapter adapter = new ConfEnergyCalculatorAdapter.Builder(confEcalc)
				.setPosInterDist(posInterGen.dist)
				.setMinimizing(true)
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

		try (CPUConfEnergyCalculator confEcalc = new CPUConfEnergyCalculator(confSpace, 1)) {

			// compute the true reference energies using the new calculator
			SimpleReferenceEnergies eref = new ErefCalculator.Builder(confEcalc)
				.setMinimize(false)
				.build()
				.calc();

			// compare to the reference energies computed using the adapted old calculator
			ConfEnergyCalculatorAdapter adapter = new ConfEnergyCalculatorAdapter.Builder(confEcalc)
				.setMinimizing(false)
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
				.setMinimizing(true)
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
				.setMinimizing(false)
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

			// compute the energy matrix, with reference energies
			SimpleReferenceEnergies eref = new ErefCalculator.Builder(confEcalc)
				.setMinimize(true)
				.build()
				.calc();
			PosInterGen posInterGen = new PosInterGen(PosInterDist.DesmetEtAl1992, eref);
			EnergyMatrix emat = new EmatCalculator.Builder(confEcalc, posInterGen)
				.setMinimize(true)
				.build()
				.calc();

			// enumerate conformations with A*, and hopefully don't crash
			// the conf space only has few hundred confs, so we can enumerate them all easily
			ConfAStarTree astar = new ConfAStarTree.Builder(emat, confSpace)
				.setTraditional()
				.build();
			int numConfs = 0;
			while (astar.nextConf() != null) {
				numConfs++;
			}

			assertThat(numConfs, is(289));
		}
	}

	public void findGMEC(Parallelism parallelism) {

		try (ConfEnergyCalculator confEcalc = ConfEnergyCalculator.build(confSpace, parallelism)) {

			// compute the energy matrix, with reference energies
			SimpleReferenceEnergies eref = new ErefCalculator.Builder(confEcalc)
				.setMinimize(true)
				.build()
				.calc();
			PosInterGen posInterGen = new PosInterGen(PosInterDist.DesmetEtAl1992, eref);
			EnergyMatrix emat = new EmatCalculator.Builder(confEcalc, posInterGen)
				.setMinimize(true)
				.build()
				.calc();

			// define the conf search function
			ConfAStarTree astar = new ConfAStarTree.Builder(emat, confSpace)
				.setTraditional()
				.build();

			// find the GMEC, hopefully without crashing
			ConfEnergyCalculatorAdapter adapter = new ConfEnergyCalculatorAdapter.Builder(confEcalc)
				.setPosInterDist(posInterGen.dist)
				.setMinimizing(true)
				.build();
			SimpleGMECFinder gmecFinder = new SimpleGMECFinder.Builder(astar, adapter).build();
			gmecFinder.find();
		}
	}
	@Test public void findGMEC_CPU1() { findGMEC(Parallelism.makeCpu(1)); }
	@Test public void findGMEC_CPU2() { findGMEC(Parallelism.makeCpu(2)); }
	@Test public void findGMEC_CPU4() { findGMEC(Parallelism.makeCpu(4)); }

	// TODO: GPU ecalcs?
	// TODO: K*?
	// TODO: SOFEA?
}
