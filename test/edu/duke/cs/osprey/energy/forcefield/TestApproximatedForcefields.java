/*
** This file is part of OSPREY 3.0
** 
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
** 
** OSPREY relies on grants for its development, and since visibility
** in the scientific literature is essential for our success, we
** ask that users of OSPREY cite our papers. See the CITING_OSPREY
** document in this distribution for more information.
** 
** Contact Info:
**    Bruce Donald
**    Duke University
**    Department of Computer Science
**    Levine Science Research Center (LSRC)
**    Durham
**    NC 27708-0129
**    USA
**    e-mail: www.cs.duke.edu/brd/
** 
** <signature of Bruce Donald>, Mar 1, 2018
** Bruce Donald, Professor of Computer Science
*/

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
import edu.duke.cs.osprey.energy.approximation.ApproximatorMatrix;
import edu.duke.cs.osprey.energy.approximation.ApproximatorMatrixCalculator;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.structure.PDBIO;
import org.junit.Test;

import java.util.function.BiConsumer;


public class TestApproximatedForcefields {

	private static final SimpleConfSpace confSpace;

	static {

		// get a strand with a variety of amino acids
		Strand strand = new Strand.Builder(PDBIO.readResource("/1CC8.ss.pdb")).build();

		strand.flexibility.get("A2" ).setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous(); // ala
		strand.flexibility.get("A11").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous(); // val
		strand.flexibility.get("A17").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous(); // gly
		strand.flexibility.get("A24").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous(); // lys
		strand.flexibility.get("A69").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous(); // ser

		// make the conf space
		confSpace = new SimpleConfSpace.Builder()
			.addStrand(strand)
			.build();
	}

	private static void withCPUConfEcalcs(EnergyPartition epart, boolean useEref, BiConsumer<ConfEnergyCalculator,ConfEnergyCalculator> f) {

		// get a (minimizing) energy calculator
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

			// calc the approximator matrix with the minimizing energy calculator
			ApproximatorMatrix amat = new ApproximatorMatrixCalculator(confEcalc).calc();

			// do the tests with non-minimizing energy calculators
			// to avoid having the accuracy tests depend on poorly conditioned numerical minimization
			EnergyCalculator ecalcRigid = new EnergyCalculator.SharedBuilder(ecalc)
				.setIsMinimizing(false)
				.build();
			ConfEnergyCalculator confEcalcRigid = new ConfEnergyCalculator.Builder(confSpace, ecalcRigid)
				.setEnergyPartition(confEcalc.epart)
				.setReferenceEnergies(confEcalc.eref)
				.build();
			ConfEnergyCalculator confEcalcRigidApprox = new ConfEnergyCalculator.Builder(confSpace, ecalcRigid)
				.setEnergyPartition(confEcalc.epart)
				.setReferenceEnergies(confEcalc.eref)
				.setApproximatorMatrix(amat)
				.setApproximationErrorBudget(1e-3)
				.build();

			f.accept(confEcalcRigid, confEcalcRigidApprox);
		}
	}

	private static void check(EnergyPartition epart, boolean useEref) {
		withCPUConfEcalcs(epart, useEref, (confEcalcRigid, confEcalcRigidApprox) -> {

			double epsilon = confEcalcRigidApprox.approximationErrorBudget;

			// compare energy matrices
			EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(confEcalcRigid).build().calcEnergyMatrix();
			EnergyMatrix ematApprox = new SimplerEnergyMatrixCalculator.Builder(confEcalcRigidApprox).build().calcEnergyMatrix();

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
				if (conf == null) {
					break;
				}

				double lowerBound = conf.getScore();
				double lowerBoundApprox = ematApprox.confE(conf.getAssignments());
				assertThat(lowerBoundApprox, isAbsolutely(lowerBound, epsilon));

				double energy = confEcalcRigid.calcEnergy(conf).getEnergy();
				double energyApprox = confEcalcRigidApprox.calcEnergy(conf).getEnergy();
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
					.setCacheFile(file)
					.calc();

				// read it out of the cache
				ApproximatorMatrix amat2 = new ApproximatorMatrixCalculator(confEcalc)
					.setCacheFile(file)
					.calc();

				// check that all the approximators are the same
				for (SimpleConfSpace.Position pos1 : confSpace.positions) {
					for (SimpleConfSpace.ResidueConf rc1 : pos1.resConfs) {

						// singles
						assertThat(amat2.get(pos1, rc1), is(amat1.get(pos1, rc1)));

						// pairs
						for (SimpleConfSpace.Position pos2 : confSpace.positions.subList(0, pos1.index)) {
							for (SimpleConfSpace.ResidueConf rc2 : pos2.resConfs) {
								assertThat(amat2.get(pos1, rc1, pos2, rc2), is(amat1.get(pos1, rc1, pos2, rc2)));
							}
						}

						// inters with fixed residues
						for (String fixedResNum : confSpace.shellResNumbers) {
							assertThat(amat2.get(pos1, rc1, fixedResNum), is(amat1.get(pos1, rc1, fixedResNum)));
						}
					}
				}
			}
		}
	}
}
