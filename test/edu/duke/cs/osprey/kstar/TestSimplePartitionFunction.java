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

package edu.duke.cs.osprey.kstar;

import static edu.duke.cs.osprey.TestBase.TempFile;
import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.*;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimpleReferenceEnergies;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.externalMemory.ExternalMemory;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.kstar.pfunc.SimplePartitionFunction;
import edu.duke.cs.osprey.kstar.pfunc.GradientDescentPfunc;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import org.junit.Test;

import java.math.BigDecimal;
import java.util.concurrent.atomic.AtomicReference;
import java.util.function.Consumer;

public class TestSimplePartitionFunction {

	public static class TestInfo {
		public ForcefieldParams ffparams;
		public Molecule mol;
		public ResidueTemplateLibrary templateLib;
		public Strand protein;
		public Strand ligand;
	}

	private static interface PfuncFactory {
		PartitionFunction make(ConfEnergyCalculator confEcalc);
	}

	public static EnergyMatrix calcEmat(ForcefieldParams ffparams, SimpleConfSpace confSpace, Parallelism parallelism) {
		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpace, ffparams)
			.setParallelism(parallelism)
			.build()
		) {

			// define conf energies
			SimpleReferenceEnergies eref = new SimplerEnergyMatrixCalculator.Builder(confSpace, ecalc)
				.build()
				.calcReferenceEnergies();
			ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(confSpace, ecalc)
				.setReferenceEnergies(eref)
				.build();

			// compute the energy matrix
			return new SimplerEnergyMatrixCalculator.Builder(confEcalc)
				.build()
				.calcEnergyMatrix();
		}
	}

	public static PartitionFunction calcPfunc(ForcefieldParams ffparams, SimpleConfSpace confSpace, Parallelism parallelism, double targetEpsilon, EnergyMatrix emat, PfuncFactory pfuncs, Consumer<PartitionFunction> pfuncComputer) {

		AtomicReference<PartitionFunction> pfuncRef = new AtomicReference<>(null);

		new EnergyCalculator.Builder(confSpace, ffparams)
			.setParallelism(parallelism)
			.use((ecalc) -> {

				// define conf energies
				SimpleReferenceEnergies eref = new SimplerEnergyMatrixCalculator.Builder(confSpace, ecalc)
					.build()
					.calcReferenceEnergies();
				ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(confSpace, ecalc)
					.setReferenceEnergies(eref)
					.build();

				// make the A* search
				ConfAStarTree astar = new ConfAStarTree.Builder(emat, confSpace)
					.setTraditional()
					.build();

				// make the partition function
				PartitionFunction pfunc = pfuncs.make(confEcalc);
				pfunc.setReportProgress(true);

				// compute pfunc for protein
				pfunc.init(astar, astar.getNumConformations(), targetEpsilon);
				pfuncComputer.accept(pfunc);
				pfuncRef.set(pfunc);
			});

		return pfuncRef.get();
	}

	private static PfuncFactory simplePfuncs = (confEcalc) -> new SimplePartitionFunction(confEcalc);
	private static PfuncFactory gdPfuncs = (confEcalc) -> new GradientDescentPfunc(confEcalc);

	public static void testStrand(ForcefieldParams ffparams, SimpleConfSpace confSpace, Parallelism parallelism, double targetEpsilon, String approxQStar, EnergyMatrix emat, PfuncFactory pfuncs) {

		PartitionFunction pfunc;

		// calculate the pfunc all at once, like for K*
		pfunc = calcPfunc(ffparams, confSpace, parallelism, targetEpsilon, emat, pfuncs, (p) -> {
			p.compute();
		});
		assertPfunc(pfunc, PartitionFunction.Status.Estimated, targetEpsilon, approxQStar);

		// calculate the pfunc in waves, like for BBK*
		pfunc = calcPfunc(ffparams, confSpace, parallelism, targetEpsilon, emat, pfuncs, (p) -> {
			while (p.getStatus().canContinue()) {
				p.compute(4);
			}
		});
		assertPfunc(pfunc, PartitionFunction.Status.Estimated, targetEpsilon, approxQStar);
	}

	public static void assertPfunc(PartitionFunction pfunc, PartitionFunction.Status status, double targetEpsilon, String approxQstar) {
		assertThat(pfunc.getStatus(), is(status));
		assertThat(pfunc.getValues().getEffectiveEpsilon(), lessThanOrEqualTo(targetEpsilon));
		double qbound = new BigDecimal(approxQstar).doubleValue()*(1.0 - targetEpsilon);
		assertThat(pfunc.getValues().qstar.doubleValue(), greaterThanOrEqualTo(qbound));
	}


	public static TestInfo make2RL0TestInfo() {

		TestInfo info = new TestInfo();

		// choose a forcefield
		info.ffparams = new ForcefieldParams();

		// choose a molecule
		info.mol = PDBIO.readFile("examples/2RL0.kstar/2RL0.min.reduce.pdb");

		// make sure all strands share the same template library
		info.templateLib = new ResidueTemplateLibrary.Builder(info.ffparams.forcefld).build();

		info.protein = new Strand.Builder(info.mol)
			.setTemplateLibrary(info.templateLib)
			.setResidues("G648", "G654")
			.build();
		info.protein.flexibility.get("G649").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		info.protein.flexibility.get("G650").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		info.protein.flexibility.get("G651").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		info.protein.flexibility.get("G654").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

		info.ligand = new Strand.Builder(info.mol)
			.setTemplateLibrary(info.templateLib)
			.setResidues("A155", "A194")
			.build();
		info.ligand.flexibility.get("A156").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		info.ligand.flexibility.get("A172").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		info.ligand.flexibility.get("A192").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		info.ligand.flexibility.get("A193").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

		return info;
	}

	private static EnergyMatrix calc2RL0ProteinEmat = null;
	public void calc2RL0Protein(PfuncFactory pfuncs, Parallelism parallelism) {
		TestInfo info = make2RL0TestInfo();
		SimpleConfSpace confSpace = new SimpleConfSpace.Builder()
				.addStrand(info.protein)
				.build();
		if (calc2RL0ProteinEmat == null) {
			calc2RL0ProteinEmat = calcEmat(info.ffparams, confSpace, parallelism);
		}
		final double targetEpsilon = 0.05;
		final String approxQStar = "4.370068e+04"; // e=0.001
		testStrand(info.ffparams, confSpace, parallelism, targetEpsilon, approxQStar, calc2RL0ProteinEmat, pfuncs);
	}
	@Test public void test2RL0ProteinSimple1Cpu() { calc2RL0Protein(simplePfuncs, Parallelism.make(1, 0, 0)); }
	@Test public void test2RL0ProteinSimple2Cpus() { calc2RL0Protein(simplePfuncs, Parallelism.make(2, 0, 0)); }
	@Test public void test2RL0ProteinSimple1GpuStream() { calc2RL0Protein(simplePfuncs, Parallelism.make(1, 1, 1)); }
	@Test public void test2RL0ProteinSimple4GpuStreams() { calc2RL0Protein(simplePfuncs, Parallelism.make(2, 1, 4)); }
	@Test public void test2RL0ProteinGD1Cpu() { calc2RL0Protein(gdPfuncs, Parallelism.make(1, 0, 0)); }
	@Test public void test2RL0ProteinGD2Cpus() { calc2RL0Protein(gdPfuncs, Parallelism.make(2, 0, 0)); }
	@Test public void test2RL0ProteinGD1GpuStream() { calc2RL0Protein(gdPfuncs, Parallelism.make(1, 1, 1)); }
	@Test public void test2RL0ProteinGD4GpuStreams() { calc2RL0Protein(gdPfuncs, Parallelism.make(2, 1, 4)); }

	private static EnergyMatrix calc2RL0LigandEmat = null;
	public void calc2RL0LigandPfunc(PfuncFactory pfuncs, Parallelism parallelism) {
		TestInfo info = make2RL0TestInfo();
		SimpleConfSpace confSpace = new SimpleConfSpace.Builder()
			.addStrand(info.ligand)
			.build();
		if (calc2RL0LigandEmat == null) {
			calc2RL0LigandEmat = calcEmat(info.ffparams, confSpace, parallelism);
		}
		final double targetEpsilon = 0.05;
		final String approxQStar = "4.467797e+30"; // e=0.001
		testStrand(info.ffparams, confSpace, parallelism, targetEpsilon, approxQStar, calc2RL0LigandEmat, pfuncs);
	}
	@Test public void test2RL0LigandSimple1Cpu() { calc2RL0LigandPfunc(simplePfuncs, Parallelism.make(1, 0, 0)); }
	@Test public void test2RL0LigandSimple2Cpus() { calc2RL0LigandPfunc(simplePfuncs, Parallelism.make(2, 0, 0)); }
	@Test public void test2RL0LigandSimple1GpuStream() { calc2RL0LigandPfunc(simplePfuncs, Parallelism.make(1, 1, 1)); }
	@Test public void test2RL0LigandSimple4GpuStreams() { calc2RL0LigandPfunc(simplePfuncs, Parallelism.make(2, 1, 4)); }
	@Test public void test2RL0LigandGD1Cpu() { calc2RL0LigandPfunc(gdPfuncs, Parallelism.make(1, 0, 0)); }
	@Test public void test2RL0LigandGD2Cpus() { calc2RL0LigandPfunc(gdPfuncs, Parallelism.make(2, 0, 0)); }
	@Test public void test2RL0LigandGD1GpuStream() { calc2RL0LigandPfunc(gdPfuncs, Parallelism.make(1, 1, 1)); }
	@Test public void test2RL0LigandGD4GpuStreams() { calc2RL0LigandPfunc(gdPfuncs, Parallelism.make(2, 1, 4)); }

	private static EnergyMatrix calc2RL0ComplexEmat = null;
	public void calc2RL0Complex(PfuncFactory pfuncs, Parallelism parallelism) {

		// NOTE: to match the old code precisely, we need to match the old conf space exactly too
		// which means we need to keep the same design position order (CCD is of course sensitive to this)
		// and also add the extra residues in the PDB file that aren't in the strands
		TestInfo info = make2RL0TestInfo();
		SimpleConfSpace confSpace = new SimpleConfSpace.Builder()
			.addStrand(new Strand.Builder(info.mol).setTemplateLibrary(info.templateLib).setResidues("A153", "A154").build())
			.addStrand(info.ligand)
			.addStrand(new Strand.Builder(info.mol).setTemplateLibrary(info.templateLib).setResidues("A195", "A241").build())
			.addStrand(new Strand.Builder(info.mol).setTemplateLibrary(info.templateLib).setResidues("G638", "G647").build())
			.addStrand(info.protein)
			.build();
		if (calc2RL0ComplexEmat == null) {
			calc2RL0ComplexEmat = calcEmat(info.ffparams, confSpace, parallelism);
		}
		final double targetEpsilon = 0.8;
		final String approxQStar = "3.5213742379e+54"; // e=0.05
		testStrand(info.ffparams, confSpace, parallelism, targetEpsilon, approxQStar, calc2RL0ComplexEmat, pfuncs);
	}
	@Test public void test2RL0ComplexSimple1Cpu() { calc2RL0Complex(simplePfuncs, Parallelism.make(1, 0, 0)); }
	@Test public void test2RL0ComplexSimple2Cpus() { calc2RL0Complex(simplePfuncs, Parallelism.make(2, 0, 0)); }
	@Test public void test2RL0ComplexSimple4Cpus() { calc2RL0Complex(simplePfuncs, Parallelism.make(4, 0, 0)); }
	@Test public void test2RL0ComplexSimple1GpuStream() { calc2RL0Complex(simplePfuncs, Parallelism.make(1, 1, 1)); }
	@Test public void test2RL0ComplexSimple4GpuStreams() { calc2RL0Complex(simplePfuncs, Parallelism.make(2, 1, 4)); }
	@Test public void test2RL0ComplexGD1Cpu() { calc2RL0Complex(gdPfuncs, Parallelism.make(1, 0, 0)); }
	@Test public void test2RL0ComplexGD2Cpus() { calc2RL0Complex(gdPfuncs, Parallelism.make(2, 0, 0)); }
	@Test public void test2RL0ComplexGD4Cpus() { calc2RL0Complex(gdPfuncs, Parallelism.make(4, 0, 0)); }
	@Test public void test2RL0ComplexGD1GpuStream() { calc2RL0Complex(gdPfuncs, Parallelism.make(1, 1, 1)); }
	@Test public void test2RL0ComplexGD4GpuStreams() { calc2RL0Complex(gdPfuncs, Parallelism.make(2, 1, 4)); }


	public static TestInfo make1GUA11TestInfo() {

		TestInfo info = new TestInfo();

		// choose a forcefield
		info.ffparams = new ForcefieldParams();

		// choose a molecule
		info.mol = PDBIO.readFile("test-resources/1gua_adj.min.pdb");

		// make sure all strands share the same template library
		info.templateLib = new ResidueTemplateLibrary.Builder(info.ffparams.forcefld).build();

		info.protein = new Strand.Builder(info.mol)
			.setTemplateLibrary(info.templateLib)
			.setResidues("1", "180")
			.build();
		info.protein.flexibility.get("21").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		info.protein.flexibility.get("24").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		info.protein.flexibility.get("25").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		info.protein.flexibility.get("27").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		info.protein.flexibility.get("29").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		info.protein.flexibility.get("40").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

		info.ligand = new Strand.Builder(info.mol)
			.setTemplateLibrary(info.templateLib)
			.setResidues("181", "215")
			.build();
		info.ligand.flexibility.get("209").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		info.ligand.flexibility.get("213").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

		return info;
	}

	private static EnergyMatrix calc1GUA11ProteinEmat = null;
	public void calc1GUA11Protein(PfuncFactory pfuncs, Parallelism parallelism) {
		TestInfo info = make1GUA11TestInfo();
		SimpleConfSpace confSpace = new SimpleConfSpace.Builder()
			.addStrand(info.protein)
			.build();
		if (calc1GUA11ProteinEmat == null) {
			calc1GUA11ProteinEmat = calcEmat(info.ffparams, confSpace, parallelism);
		}
		final double targetEpsilon = 0.9;
		final String approxQStar = "1.1838e+42"; // e = 0.1
		testStrand(info.ffparams, confSpace, parallelism, targetEpsilon, approxQStar, calc1GUA11ProteinEmat, pfuncs);
	}
	@Test public void calc1GUA11ProteinSimple() { calc1GUA11Protein(simplePfuncs, Parallelism.makeCpu(4)); }
	@Test public void calc1GUA11ProteinGD() { calc1GUA11Protein(gdPfuncs, Parallelism.makeCpu(4)); }

	private static EnergyMatrix calc1GUA11LigandEmat = null;
	public void calc1GUA11Ligand(PfuncFactory pfuncs, Parallelism parallelism) {
		TestInfo info = make1GUA11TestInfo();
		SimpleConfSpace confSpace = new SimpleConfSpace.Builder()
			.addStrand(info.ligand)
			.build();
		if (calc1GUA11LigandEmat == null) {
			calc1GUA11LigandEmat = calcEmat(info.ffparams, confSpace, parallelism);
		}
		final double targetEpsilon = 0.9;
		final String approxQStar = "2.7098e+7"; // e = 0.1
		testStrand(info.ffparams, confSpace, Parallelism.makeCpu(1), targetEpsilon, approxQStar, calc1GUA11LigandEmat, pfuncs);
	}
	@Test public void calc1GUA11LigandSimple() { calc1GUA11Ligand(simplePfuncs, Parallelism.makeCpu(4)); }
	@Test public void calc1GUA11LigandGD() { calc1GUA11Ligand(gdPfuncs, Parallelism.makeCpu(4)); }

	private static EnergyMatrix calc1GUA11ComplexEmat = null;
	public void calc1GUA11Complex(PfuncFactory pfuncs, Parallelism parallelism) {
		TestInfo info = make1GUA11TestInfo();
		SimpleConfSpace confSpace = new SimpleConfSpace.Builder()
			.addStrands(info.protein, info.ligand)
			.build();
		if (calc1GUA11ComplexEmat == null) {
			calc1GUA11ComplexEmat = calcEmat(info.ffparams, confSpace, parallelism);
		}
		final double targetEpsilon = 0.9;
		final String approxQStar = "1.1195e+66"; // e = 0.1
		testStrand(info.ffparams, confSpace, Parallelism.makeCpu(1), targetEpsilon, approxQStar, calc1GUA11ComplexEmat, pfuncs);
	}
	@Test public void calc1GUA11ComplexSimple() { calc1GUA11Complex(simplePfuncs, Parallelism.makeCpu(4)); }
	@Test public void calc1GUA11ComplexGD() { calc1GUA11Complex(gdPfuncs, Parallelism.makeCpu(4)); }

	public void calcWithConfDB(PfuncFactory pfuncs) {

		try (TempFile confdbFile = new TempFile("pfunc.conf.db")) {

			TestInfo info = make2RL0TestInfo();
			SimpleConfSpace confSpace = new SimpleConfSpace.Builder()
				.addStrand(info.ligand)
				.build();
			final Parallelism parallelism = Parallelism.makeCpu(4);
			EnergyMatrix emat = calcEmat(info.ffparams, confSpace, parallelism);
			final double targetEpsilon = 0.002;
			final String approxQStar = "4.467797e+30"; // e=0.001

			// calc the pfunc with an empty db
			PartitionFunction pfunc = calcPfunc(new ForcefieldParams(), confSpace, parallelism, targetEpsilon, emat, pfuncs, (p) -> {
				try (ConfDB confdb = new ConfDB(confSpace, confdbFile)) {
					((PartitionFunction.WithConfTable)p).setConfTable(confdb.new ConfTable("test"));
					p.compute();
				}
			});
			assertPfunc(pfunc, PartitionFunction.Status.Estimated, targetEpsilon, approxQStar);

			// the db should have stuff in it
			assertThat(confdbFile.exists(), is(true));
			try (ConfDB confdb = new ConfDB(confSpace, confdbFile)) {
				assertThat(confdb.new ConfTable("test").size(), greaterThan(0L));
			}

			// calc the pfunc with a full db
			pfunc = calcPfunc(new ForcefieldParams(), confSpace, parallelism, targetEpsilon, emat, pfuncs, (p) -> {
				try (ConfDB confdb = new ConfDB(confSpace, confdbFile)) {
					((PartitionFunction.WithConfTable)p).setConfTable(confdb.new ConfTable("test"));
					p.compute();
				}
			});
			assertPfunc(pfunc, PartitionFunction.Status.Estimated, targetEpsilon, approxQStar);
		}
	}
	@Test public void calcWithConfDBGD() { calcWithConfDB(gdPfuncs); }

	@Test
	public void withExternalMemory() {

		ExternalMemory.use(16, () -> {

			TestInfo info = make2RL0TestInfo();
			SimpleConfSpace confSpace = new SimpleConfSpace.Builder()
				.addStrand(info.protein)
				.build();

			final double targetEpsilon = 0.05;
			final String approxQStar = "4.370068e+04"; // e=0.001

			try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpace, new ForcefieldParams())
				.setParallelism(Parallelism.makeCpu(4))
				.build()) {

				// define conf energies
				SimpleReferenceEnergies eref = new SimplerEnergyMatrixCalculator.Builder(confSpace, ecalc)
					.build()
					.calcReferenceEnergies();
				ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(confSpace, ecalc)
					.setReferenceEnergies(eref)
					.build();

				EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(confEcalc)
					.build()
					.calcEnergyMatrix();

				// make the A* search
				RCs rcs = new RCs(confSpace);
				ConfAStarTree astar = new ConfAStarTree.Builder(emat, rcs)
					.setTraditional()
					.build();

				// make the partition function
				PartitionFunction pfunc = new GradientDescentPfunc(confEcalc);
				pfunc.setReportProgress(true);

				// compute pfunc for protein
				pfunc.init(astar, rcs.getNumConformations(), targetEpsilon);
				pfunc.compute();

				assertPfunc(pfunc, PartitionFunction.Status.Estimated, targetEpsilon, approxQStar);
			}
		});
	}
}
