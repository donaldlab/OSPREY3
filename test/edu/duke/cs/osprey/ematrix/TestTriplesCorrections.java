package edu.duke.cs.osprey.ematrix;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyPartition;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import org.junit.Test;

import java.util.*;

import static edu.duke.cs.osprey.TestBase.isAbsolutely;
import static edu.duke.cs.osprey.tools.Log.log;


public class TestTriplesCorrections {

	// conf spaces with < 3 positions won't get corrected, but should still be accurate
	@Test
	public void testPerfectCorrections_1CC8_2N_Traditional() {
		assertAccurateCorrections(
			makeConfSpaceNoFixed(
				"/1CC8.ss.pdb",
				"A45", // val
				"A47" // val
			),
			EnergyPartition.Traditional
		);
	}

	// AllOnPairs conf spaces with 2 positions won't get corrected either, but the bounds should be perfect anyway! =F
	@Test
	public void testPerfectCorrections_1CC8_2N_AllOnPairs() {
		assertPerfectCorrections(
			makeConfSpaceNoFixed(
				"/1CC8.ss.pdb",
				"A45", // val
				"A47" // val
			),
			EnergyPartition.AllOnPairs
		);
	}
	@Test
	public void testPerfectCorrections_1CC8_2F_AllOnPairs() {
		assertPerfectCorrections(
			makeConfSpace(
				"/1CC8.ss.pdb",
				"A45", // val
				"A47" // val
			),
			EnergyPartition.AllOnPairs
		);
	}

	// conf spaces with >= 3 positions should get corrected, and the corrections should be accurate
	@Test
	public void testPerfectCorrections_1CC8_3N_Traditional() {
		assertAccurateCorrections(
			makeConfSpaceNoFixed(
				"/1CC8.ss.pdb",
				"A26", // leu
				"A45", // val
				"A47" // val
			),
			EnergyPartition.Traditional
		);
	}
	@Test
	public void testPerfectCorrections_1CC8_3F_Traditional() {
		assertAccurateCorrections(
			makeConfSpace(
				"/1CC8.ss.pdb",
				"A26", // leu
				"A45", // val
				"A47" // val
			),
			EnergyPartition.Traditional
		);
	}

	// except corrections should be perfect for AllOnPairs with 3 positions! =D
	@Test
	public void testPerfectCorrections_1CC8_3N_AllOnPairs() {
		assertPerfectCorrections(
			makeConfSpaceNoFixed(
				"/1CC8.ss.pdb",
				"A26", // leu
				"A45", // val
				"A47" // val
			),
			EnergyPartition.AllOnPairs
		);
	}
	@Test
	public void testPerfectCorrections_1CC8_3F_AllOnPairs() {
		assertPerfectCorrections(
			makeConfSpace(
				"/1CC8.ss.pdb",
				"A26", // leu
				"A45", // val
				"A47" // val
			),
			EnergyPartition.AllOnPairs
		);
	}

	@Test
	public void testPerfectCorrections_1CC8_4N_Traditional() {
		assertAccurateCorrections(
			makeConfSpaceNoFixed(
				"/1CC8.ss.pdb",
				"A26", // leu
				"A45", // val
				"A47", // val
				"A60" // ile
			),
			EnergyPartition.Traditional
		);
	}
	@Test
	public void testPerfectCorrections_1CC8_4F_Traditional() {
		assertAccurateCorrections(
			makeConfSpace(
				"/1CC8.ss.pdb",
				"A26", // leu
				"A45", // val
				"A47", // val
				"A60" // ile
			),
			EnergyPartition.Traditional
		);
	}
	@Test
	public void testPerfectCorrections_1CC8_4N_AllOnPairs() {
		assertAccurateCorrections(
			makeConfSpaceNoFixed(
				"/1CC8.ss.pdb",
				"A26", // leu
				"A45", // val
				"A47", // val
				"A60" // ile
			),
			EnergyPartition.AllOnPairs
		);
	}
	@Test
	public void testPerfectCorrections_1CC8_4F_AllOnPairs() {
		assertAccurateCorrections(
			makeConfSpace(
				"/1CC8.ss.pdb",
				"A26", // leu
				"A45", // val
				"A47", // val
				"A60" // ile
			),
			EnergyPartition.AllOnPairs
		);
	}

	private static SimpleConfSpace makeConfSpace(String pdbPath, String ... resNums) {

		// share the template library
		ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder().build();

		// build strands, one for each residue
		// so we don't end up with any fixed residues
		Molecule pdb = PDBIO.readResource(pdbPath);
		Strand strand = new Strand.Builder(pdb)
				.setTemplateLibrary(templateLib)
				.build();
		for (String resNum : resNums) {
			strand.flexibility.get(resNum)
				.setLibraryRotamers(Strand.WildType)
				.addWildTypeRotamers()
				.setContinuous();
		}

		return new SimpleConfSpace.Builder()
			.addStrand(strand)
			.build();
	}

	private static SimpleConfSpace makeConfSpaceNoFixed(String pdbPath, String ... resNums) {

		// share the template library
		ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder().build();

		// build strands, one for each residue
		// so we don't end up with any fixed residues
		Molecule pdb = PDBIO.readResource(pdbPath);
		List<Strand> strands = new ArrayList<>();
		for (String resNum : resNums) {
			Strand strand = new Strand.Builder(pdb)
				.setTemplateLibrary(templateLib)
				.setResidues(resNum, resNum)
				.build();
			strands.add(strand);
			strand.flexibility.get(resNum)
				.setLibraryRotamers(Strand.WildType)
				.addWildTypeRotamers()
				.setContinuous();
		}

		return new SimpleConfSpace.Builder()
			.addStrands(strands)
			.build();
	}

	private static void assertPerfectCorrections(SimpleConfSpace confSpace, EnergyPartition epart) {
		forEachTopConf(confSpace, epart, (index, conf, energy, lowerBound, correctedBound) -> {
			assertThat(energy - correctedBound, isAbsolutely(0, 1e-4));
		});
	}

	private static void assertAccurateCorrections(SimpleConfSpace confSpace, EnergyPartition epart) {
		forEachTopConf(confSpace, epart, (index, conf, energy, lowerBound, correctedBound) -> {
			assertThat(energy - correctedBound, greaterThanOrEqualTo(-1e-4));
		});
	}

	interface ConfListener {
		void onConf(int index, int[] conf, double energy, double lowerBound, double correctedBound);
	}

	private static void forEachTopConf(SimpleConfSpace confSpace, EnergyPartition epart, ConfListener listener) {

		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpace, new ForcefieldParams())
			.setParallelism(Parallelism.makeCpu(4))
			.build()) {

			ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(confSpace, ecalc)
				.setEnergyPartition(epart)
				.build();

			// calc an emat
			EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(confEcalc)
				.setTripleCorrectionThreshold(Double.POSITIVE_INFINITY) // TODO: try lower?
				.build()
				.calcEnergyMatrix();

			// get the top 10 confs by energy lower bound
			ConfAStarTree astar = new ConfAStarTree.Builder(emat, confSpace)
				.setTraditional()
				.build();
			for (int i=0; i<10; i++) {

				// get the conf and the lower bound
				ConfSearch.ScoredConf conf = astar.nextConf();
				if (conf == null) {
					break;
				}

				// calculate the energy
				ConfSearch.EnergiedConf econf = confEcalc.calcEnergy(conf);

				// calculate the corrected bound
				double[] corrected = { econf.getScore() };
				emat.forEachHigherOrderTupleIn(econf.getAssignments(), (tuple, tupleEnergy) -> {
					corrected[0] += tupleEnergy;
				});

				// TEMP
				log("conf %2d   e=%9.3f   lb=%9.3f (%6.3f)   clb=%9.3f (%6.3f)",
					i,
					econf.getEnergy(),
					econf.getScore(), econf.getEnergy() - econf.getScore(),
					corrected[0], econf.getEnergy() - corrected[0]
				);

				listener.onConf(
					i,
					econf.getAssignments(),
					econf.getEnergy(),
					econf.getScore(),
					corrected[0]
				);
			}
		}
	}
}
