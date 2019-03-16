package edu.duke.cs.osprey.energy.forcefield;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;
import static edu.duke.cs.osprey.TestBase.*;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimpleReferenceEnergies;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyPartition;
import edu.duke.cs.osprey.energy.approximation.ApproximatedObjectiveFunction.Approximator;
import edu.duke.cs.osprey.energy.approximation.ApproximatorMatrix;
import edu.duke.cs.osprey.energy.approximation.ApproximatorMatrixCalculator;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.structure.PDBIO;
import org.junit.Test;

import java.util.List;
import java.util.function.BiConsumer;


public class TestApproximatedForcefields {

	private static final SimpleConfSpace confSpace;

	static {

		// get a strand with a variety of amino acids
		Strand strand = new Strand.Builder(PDBIO.readResource("/1CC8.ss.pdb")).build();
		strand.flexibility.get("A16").setLibraryRotamers("GLY").setContinuous();
		strand.flexibility.get("A17").setLibraryRotamers("ALA").setContinuous();
		strand.flexibility.get("A18").setLibraryRotamers("VAL").setContinuous();
		strand.flexibility.get("A19").setLibraryRotamers("SER").setContinuous();
		strand.flexibility.get("A20").setLibraryRotamers("LYS").setContinuous();

		// make the conf space
		confSpace = new SimpleConfSpace.Builder()
			.addStrand(strand)
			.build();
	}

	private static void withCPUConfEcalcs(EnergyPartition epart, boolean useEref, BiConsumer<ConfEnergyCalculator,ConfEnergyCalculator> f) {

		// get an energy calculator
		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpace, new ForcefieldParams())
			.setParallelism(Parallelism.makeCpu(4))
			.build()
		) {

			SimpleReferenceEnergies eref = null;
			if (useEref) {
				eref = new SimplerEnergyMatrixCalculator.Builder(confSpace, ecalc)
					.build()
					.calcReferenceEnergies();
			}

			ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(confSpace, ecalc)
				.setEnergyPartition(epart)
				.setReferenceEnergies(eref)
				.build();

			// calc the approximator matrix
			ApproximatorMatrix amat = new ApproximatorMatrixCalculator(confEcalc)
				.setNumSamplesPerDoF(5)
				.calc();

			// check the sorted orders
			checkAmatOrders(amat);

			ConfEnergyCalculator confEcalcApprox = new ConfEnergyCalculator.Builder(confSpace, ecalc)
				.setEnergyPartition(confEcalc.epart)
				.setReferenceEnergies(confEcalc.eref)
				.setApproximatorMatrix(amat)
				.setApproximationErrorBudget(1e-3)
				.build();

			f.accept(confEcalc, confEcalcApprox);
		}
	}

	private static void checkAmatOrders(ApproximatorMatrix amat) {

		for (SimpleConfSpace.Position pos : confSpace.positions) {
			for (SimpleConfSpace.ResidueConf rc : pos.resConfs) {

				List<ApproximatorMatrix.Entry> entries = amat.get(pos.index, rc.index);

				if (entries.size() > 1) {

					// amats should store their approximators in order of weakly increasing error
					for (int i=1; i<entries.size(); i++) {
						double errori = entries.get(i).approximator.error();
						double errorim1 = entries.get(i - 1).approximator.error();
						assertThat(errori, greaterThanOrEqualTo(errorim1));
					}
				}
			}
		}
	}

	private static void check(EnergyPartition epart, boolean useEref) {
		withCPUConfEcalcs(epart, useEref, (confEcalc, confEcalcApprox) -> {

			double epsilon = confEcalcApprox.approximationErrorBudget;

			// HACKHACK: the minimizer is somewhat ill-conditioned in general
			// small changes in the energies can cause it to follow slightly different paths
			// so unfortunately we need to add some slack to our error budget
			// to deal with the real-world side-effects of using an imperfect minimizer
			epsilon *= 2.5;

			// compare energy matrices
			EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(confEcalc).build().calcEnergyMatrix();
			EnergyMatrix ematApprox = new SimplerEnergyMatrixCalculator.Builder(confEcalcApprox).build().calcEnergyMatrix();

			for (SimpleConfSpace.Position pos1 : confSpace.positions) {
				for (SimpleConfSpace.ResidueConf rc1 : pos1.resConfs) {

					double energy = emat.getEnergy(pos1, rc1);
					double energyApprox = ematApprox.getEnergy(pos1, rc1);
					assertThat(
						String.format("%d:%d", pos1.index, rc1.index),
						energyApprox,
						isAbsolutely(energy, epsilon)
					);

					for (SimpleConfSpace.Position pos2 : confSpace.positions.subList(0, pos1.index)) {
						for (SimpleConfSpace.ResidueConf rc2 : pos2.resConfs) {

							energy = emat.getEnergy(pos1, rc1, pos2, rc2);
							energyApprox = ematApprox.getEnergy(pos1, rc1, pos2, rc2);
							assertThat(
								String.format("%d:%d,%d:%d", pos1.index, rc1.index, pos2.index, rc2.index),
								energyApprox,
								isAbsolutely(energy, epsilon)
							);
						}
					}
				}
			}

			// compare full conformations
			ConfAStarTree astar = new ConfAStarTree.Builder(emat, confSpace)
				.setTraditional()
				.build();
			for (int i=0; i<100; i++) {

				ConfSearch.ScoredConf conf = astar.nextConf();

				double energy = confEcalc.calcEnergy(conf).getEnergy();
				double energyApprox = confEcalcApprox.calcEnergy(conf).getEnergy();
				assertThat(energyApprox, isAbsolutely(energy, epsilon));
			}
		});
	}

	@Test public void traditional() { check(EnergyPartition.Traditional, false); }
	@Test public void allOnPairs() { check(EnergyPartition.AllOnPairs, false); }
	@Test public void traditional_eref() { check(EnergyPartition.Traditional, true); }
	@Test public void allOnPairs_eref() { check(EnergyPartition.AllOnPairs, true); }

	@Test
	public void ioRoundtrip() {

		// get an energy calculator
		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpace, new ForcefieldParams())
			.setParallelism(Parallelism.makeCpu(4))
			.build()
		) {
			ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(confSpace, ecalc).build();

			try (TempFile file = new TempFile("amat")) {

				// calc the approximator matrix
				ApproximatorMatrix amat1 = new ApproximatorMatrixCalculator(confEcalc)
					.setNumSamplesPerDoF(3)
					.setCacheFile(file)
					.calc();

				// read it out of the cache
				ApproximatorMatrix amat2 = new ApproximatorMatrixCalculator(confEcalc)
					.setCacheFile(file)
					.calc();

				for (String fixedResNum : confSpace.shellResNumbers) {
					for (SimpleConfSpace.Position pos1 : confSpace.positions) {
						for (SimpleConfSpace.ResidueConf rc1 : pos1.resConfs) {
							Approximator.Addable approximator1 = amat1.get(pos1.index, rc1.index, fixedResNum);
							Approximator.Addable approximator2 = amat2.get(pos1.index, rc1.index, fixedResNum);
							assertThat(approximator1, is(approximator2));
						}
					}
				}
			}
		}
	}
}
