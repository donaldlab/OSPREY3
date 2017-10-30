package edu.duke.cs.osprey.energy;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.*;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.energy.forcefield.ResPairCache;
import edu.duke.cs.osprey.energy.forcefield.ResidueForcefieldEnergy;
import edu.duke.cs.osprey.energy.forcefield.TestForcefieldEnergy;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.structure.AtomConnectivity;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.structure.Residues;
import org.junit.BeforeClass;
import org.junit.Test;

import static edu.duke.cs.osprey.TestBase.isAbsolutely;
import static org.junit.Assert.assertThat;

public class TestResidueForcefieldBreakdown {

	@BeforeClass
	public static void before() {
		TestForcefieldEnergy.before();
	}

	public void checkByResidue(Residues residues) {

		// use all pairs residue interactions
		ResidueInteractions inters = TestForcefieldEnergy.IntersType.AllPairs.makeInters(residues);

		// build the energy function
		ForcefieldParams ffparams = new ForcefieldParams();
		AtomConnectivity connectivity = new AtomConnectivity.Builder()
			.addTemplates(residues)
			.setParallelism(Parallelism.makeCpu(4))
			.build();
		ResPairCache resPairCache = new ResPairCache(ffparams, connectivity);
		ResidueForcefieldEnergy efunc = new ResidueForcefieldEnergy(resPairCache, inters, residues);

		for (TestForcefieldEnergy.FFType fftype : TestForcefieldEnergy.FFType.values()) {

			// add up all the breakdown energies
			ResidueForcefieldBreakdown.ByResidue breakdown = new ResidueForcefieldBreakdown.ByResidue(efunc);

			double sum = 0.0;
			for (ResidueForcefieldBreakdown.Type type : ResidueForcefieldBreakdown.Type.atomics()) {
				sum += breakdown.breakdownForcefield(type).sum();
			}

			// should match the all energy
			double allEnergy = breakdown.breakdownForcefield(ResidueForcefieldBreakdown.Type.All).sum();
			assertThat("forcefield type: " + fftype,
				sum, isAbsolutely(allEnergy, 1e-12)
			);

			// should match the total energy
			double totalEnergy = efunc.getEnergy();
			assertThat("forcefield type: " + fftype,
				sum, isAbsolutely(totalEnergy, 1e-12)
			);
		}
	}

	@Test
	public void singleGly() {
		TestForcefieldEnergy.TestResidues r = new TestForcefieldEnergy.TestResidues();
		checkByResidue(new Residues(r.gly15));
	}

	@Test
	public void glyPair() {
		TestForcefieldEnergy.TestResidues r = new TestForcefieldEnergy.TestResidues();
		checkByResidue(new Residues(r.gly06, r.gly15));
	}

	@Test
	public void glySerPair() {
		TestForcefieldEnergy.TestResidues r = new TestForcefieldEnergy.TestResidues();
		checkByResidue(new Residues(r.gly15, r.ser17));
	}

	@Test
	public void trpPair() {
		TestForcefieldEnergy.TestResidues r = new TestForcefieldEnergy.TestResidues();
		checkByResidue(new Residues(r.trp18, r.trp25));
	}

	@Test
	public void the4Residues() {
		TestForcefieldEnergy.TestResidues r = new TestForcefieldEnergy.TestResidues();
		checkByResidue(new Residues(r.gly15, r.ser17, r.trp18, r.trp25));
	}

	@Test
	public void the6Residues() {
		TestForcefieldEnergy.TestResidues r = new TestForcefieldEnergy.TestResidues();
		checkByResidue(new Residues(r.gly15, r.ser17, r.trp18, r.trp25, r.arg22, r.ala24));
	}

	@Test
	public void the10Residues() {
		TestForcefieldEnergy.TestResidues r = new TestForcefieldEnergy.TestResidues();
		checkByResidue(new Residues(r.gly15, r.ser17, r.trp18, r.trp25, r.arg22, r.ala24, r.ile26, r.phe31, r.arg32, r.glu34));
	}

	@Test
	public void the14Residues() {
		TestForcefieldEnergy.TestResidues r = new TestForcefieldEnergy.TestResidues();
		checkByResidue(new Residues(r.gly15, r.ser17, r.trp18, r.trp25, r.arg22, r.ala24, r.ile26, r.phe31, r.arg32, r.glu34, r.val36, r.leu39, r.trp47, r.leu48));
	}

	@Test
	public void the24Residues() {
		TestForcefieldEnergy.TestResidues r = new TestForcefieldEnergy.TestResidues();
		checkByResidue(new Residues(
			r.gly06, r.gly15, r.ser17, r.trp18, r.trp25, r.arg22, r.ala24, r.ile26, r.phe31, r.arg32, r.glu34, r.val36,
			r.leu39, r.trp47, r.leu48, r.ile53, r.arg55, r.val56, r.leu57, r.ile59, r.val62, r.leu64, r.val65, r.met66
		));
	}

	private static interface ConfEcalcFactory {
		ConfEnergyCalculator.Builder update(ConfEnergyCalculator.Builder builder, SimpleConfSpace confSpace, EnergyCalculator ecalc);
	}

	public void checkByPosition1CC8(ConfEcalcFactory confEcalcFactory) {

		// get a conf space
		Strand strand = new Strand.Builder(PDBIO.readFile("examples/python.GMEC/1CC8.ss.pdb")).build();
		strand.flexibility.get("A40").setLibraryRotamers(Strand.WildType).setContinuous();
		strand.flexibility.get("A41").setLibraryRotamers(Strand.WildType).setContinuous();
		strand.flexibility.get("A42").setLibraryRotamers(Strand.WildType).setContinuous();
		SimpleConfSpace confSpace = new SimpleConfSpace.Builder().addStrand(strand).build();

		// get the energy calculators
		ForcefieldParams ffparams = new ForcefieldParams();
		EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpace, ffparams)
			.setType(EnergyCalculator.Type.Cpu)
			.setParallelism(Parallelism.makeCpu(4))
			.build();
		ConfEnergyCalculator confEcalc = confEcalcFactory
			.update(new ConfEnergyCalculator.Builder(confSpace, ecalc), confSpace, ecalc)
			.build();

		// get the energy matrix
		EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(confEcalc)
			.build()
			.calcEnergyMatrix();

		// get a low-energy conf
		ConfSearch tree = new ConfAStarTree.Builder(emat, confSpace).build();
		ConfSearch.ScoredConf conf = tree.nextConf();

		// breakdown the forcefield
		ResidueForcefieldBreakdown.ByPosition breakdown = new ResidueForcefieldBreakdown.ByPosition(confEcalc, conf.getAssignments());

		// add up the pieces
		double sum = 0.0;
		for (ResidueForcefieldBreakdown.Type type : ResidueForcefieldBreakdown.Type.atomics()) {
			sum += breakdown.breakdownForcefield(type).sum();
		}

		// should match the all energy
		double allEnergy = breakdown.breakdownForcefield(ResidueForcefieldBreakdown.Type.All).sum();
		assertThat(sum, isAbsolutely(allEnergy, 1e-12));

		// should match the minimized energy
		assertThat(sum, isAbsolutely(breakdown.epmol.energy, 1e-12));
	}

	@Test
	public void checkByPosition1CC8Traditional() {
		checkByPosition1CC8((builder, confSpace, ecalc) -> {
			return builder
				.setEnergyPartition(EnergyPartition.Traditional);
		});
	}

	@Test
	public void checkByPosition1CC8TraditionalEref() {
		checkByPosition1CC8((builder, confSpace, ecalc) -> {
			return builder
				.setEnergyPartition(EnergyPartition.Traditional)
				.setReferenceEnergies(
					new SimplerEnergyMatrixCalculator.Builder(confSpace, ecalc)
						.build()
						.calcReferenceEnergies()
				);
		});
	}

	@Test
	public void checkByPosition1CC8AllOnPairs() {
		checkByPosition1CC8((builder, confSpace, ecalc) -> {
			return builder
				.setEnergyPartition(EnergyPartition.AllOnPairs);
		});
	}

	@Test
	public void checkByPosition1CC8AllOnPairsEref() {
		checkByPosition1CC8((builder, confSpace, ecalc) -> {
			return builder
				.setEnergyPartition(EnergyPartition.AllOnPairs)
				.setReferenceEnergies(
					new SimplerEnergyMatrixCalculator.Builder(confSpace, ecalc)
						.build()
						.calcReferenceEnergies()
				);
		});
	}
}
