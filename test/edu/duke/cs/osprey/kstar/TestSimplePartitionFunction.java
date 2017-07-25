package edu.duke.cs.osprey.kstar;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimpleReferenceEnergies;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.gmec.ConfSearchFactory;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.kstar.pfunc.SimplePartitionFunction;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import org.junit.Test;

import java.math.BigDecimal;
import java.util.concurrent.atomic.AtomicReference;

public class TestSimplePartitionFunction {

	public static class TestInfo {
		public ForcefieldParams ffparams;
		public Molecule mol;
		public ResidueTemplateLibrary templateLib;
		public Strand protein;
		public Strand ligand;
	}

	public static SimplePartitionFunction calcPfunc(ForcefieldParams ffparams, SimpleConfSpace confSpace, Parallelism parallelism, double targetEpsilon) {

		AtomicReference<SimplePartitionFunction> pfuncRef = new AtomicReference<>(null);

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

				// compute the energy matrix
				EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(confEcalc)
						.build()
						.calcEnergyMatrix();

				// don't really need pruning here (A* is plenty fast enough already), so use NOP pmat
				PruningMatrix pmat = new PruningMatrix(confSpace, 0);

				// make the partition function
				ConfSearchFactory astarFactory = (ematArg, pmatArg) -> new ConfAStarTree.Builder(ematArg, pmatArg)
						.setTraditional()
						.build();
				SimplePartitionFunction pfunc = new SimplePartitionFunction(emat, pmat, astarFactory, confEcalc);
				pfunc.setReportProgress(true);

				// compute pfunc for protein
				pfunc.init(targetEpsilon);
				pfunc.compute();
				pfuncRef.set(pfunc);
			});

		return pfuncRef.get();
	}

	public static void testStrand(ForcefieldParams ffparams, SimpleConfSpace confSpace, Parallelism parallelism, double targetEpsilon, String approxQStar) {
		PartitionFunction pfunc = calcPfunc(ffparams, confSpace, parallelism, targetEpsilon);
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
		info.templateLib = new ResidueTemplateLibrary.Builder(info.ffparams.forcefld)
			.addMoleculeForWildTypeRotamers(info.mol)
			.build();

		info.protein = new Strand.Builder(info.mol)
			.setTemplateLibrary(info.templateLib)
			.setResidues("648", "654")
			.build();
		info.protein.flexibility.get(649).setLibraryRotamers().addWildTypeRotamers().setContinuous();
		info.protein.flexibility.get(650).setLibraryRotamers().addWildTypeRotamers().setContinuous();
		info.protein.flexibility.get(651).setLibraryRotamers().addWildTypeRotamers().setContinuous();
		info.protein.flexibility.get(654).setLibraryRotamers().addWildTypeRotamers().setContinuous();

		info.ligand = new Strand.Builder(info.mol)
			.setTemplateLibrary(info.templateLib)
			.setResidues("155", "194")
			.build();
		info.ligand.flexibility.get(156).setLibraryRotamers().addWildTypeRotamers().setContinuous();
		info.ligand.flexibility.get(172).setLibraryRotamers().addWildTypeRotamers().setContinuous();
		info.ligand.flexibility.get(192).setLibraryRotamers().addWildTypeRotamers().setContinuous();
		info.ligand.flexibility.get(193).setLibraryRotamers().addWildTypeRotamers().setContinuous();

		return info;
	}

	public void calc2RL0Protein(Parallelism parallelism) {
		TestInfo info = make2RL0TestInfo();
		SimpleConfSpace confSpace = new SimpleConfSpace.Builder()
				.addStrand(info.protein)
				.build();

		final double targetEpsilon = 0.05;
		final String approxQStar = "4.370068e+04"; // e=0.001
		testStrand(info.ffparams, confSpace, parallelism, targetEpsilon, approxQStar);
	}
	@Test public void test2RL0Protein1Cpu() { calc2RL0Protein(Parallelism.make(1, 0, 0)); }
	@Test public void test2RL0Protein2Cpus() { calc2RL0Protein(Parallelism.make(2, 0, 0)); }
	@Test public void test2RL0Protein1GpuStream() { calc2RL0Protein(Parallelism.make(1, 1, 1)); }
	@Test public void test2RL0Protein4GpuStreams() { calc2RL0Protein(Parallelism.make(2, 1, 4)); }

	public void calc2RL0LigandPfunc(Parallelism parallelism) {
		TestInfo info = make2RL0TestInfo();
		SimpleConfSpace confSpace = new SimpleConfSpace.Builder()
			.addStrand(info.ligand)
			.build();
		final double targetEpsilon = 0.05;
		final String approxQStar = "4.467797e+30"; // e=0.001
		testStrand(info.ffparams, confSpace, parallelism, targetEpsilon, approxQStar);
	}
	@Test public void test2RL0Ligand1Cpu() { calc2RL0LigandPfunc(Parallelism.make(1, 0, 0)); }
	@Test public void test2RL0Ligand2Cpus() { calc2RL0LigandPfunc(Parallelism.make(2, 0, 0)); }
	@Test public void test2RL0Ligand1GpuStream() { calc2RL0LigandPfunc(Parallelism.make(1, 1, 1)); }
	@Test public void test2RL0Ligand4GpuStreams() { calc2RL0LigandPfunc(Parallelism.make(2, 1, 4)); }

	public void calc2RL0Complex(Parallelism parallelism) {

		// NOTE: to match the old code precisely, we need to match the old conf space exactly too
		// which means we need to keep the same design position order (CCD is of course sensitive to this)
		// and also add the extra residues in the PDB file that aren't in the strands
		TestInfo info = make2RL0TestInfo();
		SimpleConfSpace confSpace = new SimpleConfSpace.Builder()
			.addStrand(new Strand.Builder(info.mol).setTemplateLibrary(info.templateLib).setResidues(153, 154).build())
			.addStrand(info.ligand)
			.addStrand(new Strand.Builder(info.mol).setTemplateLibrary(info.templateLib).setResidues(195, 241).build())
			.addStrand(new Strand.Builder(info.mol).setTemplateLibrary(info.templateLib).setResidues(638, 647).build())
			.addStrand(info.protein)
			.build();

		final double targetEpsilon = 0.8;
		final String approxQStar = "3.5213742379e+54"; // e=0.05
		testStrand(info.ffparams, confSpace, parallelism, targetEpsilon, approxQStar);
	}
	@Test public void test2RL0Complex1Cpu() { calc2RL0Complex(Parallelism.make(1, 0, 0)); }
	@Test public void test2RL0Complex2Cpus() { calc2RL0Complex(Parallelism.make(2, 0, 0)); }
	@Test public void test2RL0Complex1GpuStream() { calc2RL0Complex(Parallelism.make(1, 1, 1)); }
	@Test public void test2RL0Complex4GpuStreams() { calc2RL0Complex(Parallelism.make(2, 1, 4)); }


	public static TestInfo make1GUA11TestInfo() {

		TestInfo info = new TestInfo();

		// choose a forcefield
		info.ffparams = new ForcefieldParams();

		// choose a molecule
		info.mol = PDBIO.readFile("test-resources/1gua_adj.min.pdb");

		// make sure all strands share the same template library
		info.templateLib = new ResidueTemplateLibrary.Builder(info.ffparams.forcefld)
			.addMoleculeForWildTypeRotamers(info.mol)
			.build();

		info.protein = new Strand.Builder(info.mol)
			.setTemplateLibrary(info.templateLib)
			.setResidues("1", "180")
			.build();
		info.protein.flexibility.get(21).setLibraryRotamers().addWildTypeRotamers().setContinuous();
		info.protein.flexibility.get(24).setLibraryRotamers().addWildTypeRotamers().setContinuous();
		info.protein.flexibility.get(25).setLibraryRotamers().addWildTypeRotamers().setContinuous();
		info.protein.flexibility.get(27).setLibraryRotamers().addWildTypeRotamers().setContinuous();
		info.protein.flexibility.get(29).setLibraryRotamers().addWildTypeRotamers().setContinuous();
		info.protein.flexibility.get(40).setLibraryRotamers().addWildTypeRotamers().setContinuous();

		info.ligand = new Strand.Builder(info.mol)
			.setTemplateLibrary(info.templateLib)
			.setResidues("181", "215")
			.build();
		info.ligand.flexibility.get(209).setLibraryRotamers().addWildTypeRotamers().setContinuous();
		info.ligand.flexibility.get(213).setLibraryRotamers().addWildTypeRotamers().setContinuous();

		return info;
	}

	@Test
	public void calc1GUA11Protein() {
		TestInfo info = make1GUA11TestInfo();
		SimpleConfSpace confSpace = new SimpleConfSpace.Builder()
			.addStrand(info.protein)
			.build();

		final double targetEpsilon = 0.1;
		final String approxQStar = "1.1838e+42"; // e = 0.1
		testStrand(info.ffparams, confSpace, Parallelism.makeCpu(1), targetEpsilon, approxQStar);
	}

	@Test
	public void calc1GUA11Ligand() {
		TestInfo info = make1GUA11TestInfo();
		SimpleConfSpace confSpace = new SimpleConfSpace.Builder()
			.addStrand(info.ligand)
			.build();

		final double targetEpsilon = 0.1;
		final String approxQStar = "2.7098e+7"; // e = 0.1
		testStrand(info.ffparams, confSpace, Parallelism.makeCpu(1), targetEpsilon, approxQStar);
	}

	@Test
	public void calc1GUA11Complex() {
		TestInfo info = make1GUA11TestInfo();
		SimpleConfSpace confSpace = new SimpleConfSpace.Builder()
			.addStrands(info.protein, info.ligand)
			.build();

		final double targetEpsilon = 0.1;
		final String approxQStar = "1.1195e+66"; // e = 0.1
		testStrand(info.ffparams, confSpace, Parallelism.makeCpu(1), targetEpsilon, approxQStar);
	}
}
