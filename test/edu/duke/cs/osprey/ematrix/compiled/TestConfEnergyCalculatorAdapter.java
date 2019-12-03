package edu.duke.cs.osprey.ematrix.compiled;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import edu.duke.cs.osprey.confspace.compiled.ConfSpace;
import edu.duke.cs.osprey.confspace.compiled.PosInterDist;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimpleReferenceEnergies;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.compiled.CPUConfEnergyCalculator;
import edu.duke.cs.osprey.energy.compiled.ConfEnergyCalculatorAdapter;
import edu.duke.cs.osprey.energy.compiled.PosInterGen;
import edu.duke.cs.osprey.tools.FileTools;
import org.junit.Test;


public class TestConfEnergyCalculatorAdapter {

	private static ConfSpace confSpace = new ConfSpace(FileTools.readResourceBytes("/confSpaces/dipeptide.5hydrophobic.ccs.toml.xz"));
	private static CPUConfEnergyCalculator confEcalc = new CPUConfEnergyCalculator(confSpace);

	@Test
	public void energyMatrixRigid() {

		// compute the true energy matrix using the new calculator
		PosInterGen posInterGen = new PosInterGen(PosInterDist.DesmetEtAl1992, null);
		EnergyMatrix emat = new EmatCalculator.Builder(confEcalc, posInterGen)
			.setMinimize(false)
			.build()
			.calc();

		// compare to the energy matrix computed using the adapted old calculator
		ConfEnergyCalculatorAdapter adapter = new ConfEnergyCalculatorAdapter(confEcalc, posInterGen, false);
		EnergyMatrix adaptedEmat = new SimplerEnergyMatrixCalculator.Builder(adapter)
			.build()
			.calcEnergyMatrix();

		assertThat(adaptedEmat, is(emat));
	}

	@Test
	public void energyMatrixMinimized() {

		// compute the true energy matrix using the new calculator
		PosInterGen posInterGen = new PosInterGen(PosInterDist.DesmetEtAl1992, null);
		EnergyMatrix emat = new EmatCalculator.Builder(confEcalc, posInterGen)
			.setMinimize(true)
			.build()
			.calc();

		// compare to the energy matrix computed using the adapted old calculator
		ConfEnergyCalculatorAdapter adapter = new ConfEnergyCalculatorAdapter(confEcalc, posInterGen, true);
		EnergyMatrix adaptedEmat = new SimplerEnergyMatrixCalculator.Builder(adapter)
			.build()
			.calcEnergyMatrix();

		assertThat(adaptedEmat, is(emat));
	}

	@Test
	public void referenceEnergyRigid() {

		// compute the true reference energies using the new calculator
		SimpleReferenceEnergies eref = new ErefCalculator.Builder(confEcalc)
			.setMinimize(false)
			.build()
			.calc();

		// compare to the reference energies computed using the adapted old calculator
		PosInterGen posInterGen = new PosInterGen(PosInterDist.DesmetEtAl1992, null);
		ConfEnergyCalculatorAdapter adapter = new ConfEnergyCalculatorAdapter(confEcalc, posInterGen, false);
		SimpleReferenceEnergies adaptedEref = new SimplerEnergyMatrixCalculator.Builder(adapter)
			.build()
			.calcReferenceEnergies();

		assertThat(adaptedEref, is(eref));
	}

	@Test
	public void referenceEnergyMinimized() {

		// compute the true reference energies using the new calculator
		SimpleReferenceEnergies eref = new ErefCalculator.Builder(confEcalc)
			.setMinimize(true)
			.build()
			.calc();

		// compare to the reference energies computed using the adapted old calculator
		PosInterGen posInterGen = new PosInterGen(PosInterDist.DesmetEtAl1992, null);
		ConfEnergyCalculatorAdapter adapter = new ConfEnergyCalculatorAdapter(confEcalc, posInterGen, true);

		SimpleReferenceEnergies adaptedEref = new SimplerEnergyMatrixCalculator.Builder(adapter)
			.build()
			.calcReferenceEnergies();

		assertThat(adaptedEref, is(eref));
	}

	@Test
	public void energyMatrixWithEref() {

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
		ConfEnergyCalculatorAdapter adapter = new ConfEnergyCalculatorAdapter(confEcalc, posInterGen, false);
		EnergyMatrix adaptedEmat = new SimplerEnergyMatrixCalculator.Builder(adapter)
			.build()
			.calcEnergyMatrix();

		assertThat(adaptedEmat, is(emat));
	}

	// TODO: parallelism?
	// TODO: GMEC finder?
	// TODO: K*?
	// TODO: SOFEA?
}
