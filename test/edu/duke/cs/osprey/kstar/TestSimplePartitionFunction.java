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
import org.junit.BeforeClass;
import org.junit.Test;

import java.math.BigDecimal;
import java.util.Arrays;
import java.util.concurrent.atomic.AtomicReference;

public class TestSimplePartitionFunction {

	private static ForcefieldParams ffparams;
	private static Molecule mol;
	private static ResidueTemplateLibrary templateLib;
	private static Strand protein;
	private static Strand ligand;

	@BeforeClass
	public static void beforeClass() {

		// choose a forcefield
		ffparams = new ForcefieldParams();

		// choose a molecule
		mol = PDBIO.readFile("examples/2RL0.kstar/2RL0.min.reduce.pdb");

		// make sure all strands share the same template library
		templateLib = new ResidueTemplateLibrary.Builder(ffparams.forcefld)
				.addMoleculeForWildTypeRotamers(mol)
				.build();

		protein = new Strand.Builder(mol)
				.setTemplateLibrary(templateLib)
				.setResidues("648", "654")
				.build();
		protein.flexibility.get(649).setLibraryRotamers().addWildTypeRotamers().setContinuous();
		protein.flexibility.get(650).setLibraryRotamers().addWildTypeRotamers().setContinuous();
		protein.flexibility.get(651).setLibraryRotamers().addWildTypeRotamers().setContinuous();
		protein.flexibility.get(654).setLibraryRotamers().addWildTypeRotamers().setContinuous();

		ligand = new Strand.Builder(mol)
				.setTemplateLibrary(templateLib)
				.setResidues("155", "194")
				.build();
		ligand.flexibility.get(156).setLibraryRotamers().addWildTypeRotamers().setContinuous();
		ligand.flexibility.get(172).setLibraryRotamers().addWildTypeRotamers().setContinuous();
		ligand.flexibility.get(192).setLibraryRotamers().addWildTypeRotamers().setContinuous();
		ligand.flexibility.get(193).setLibraryRotamers().addWildTypeRotamers().setContinuous();
	}

	private SimplePartitionFunction calcPfunc(SimpleConfSpace confSpace, Parallelism parallelism, double targetEpsilon) {

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

	public void calcProteinPfunc(Parallelism parallelism) {

		SimpleConfSpace confSpace = new SimpleConfSpace.Builder()
				.addStrand(protein)
				.build();

		double targetEpsilon = 0.05;
		PartitionFunction pfunc = calcPfunc(confSpace, parallelism, targetEpsilon);
		assertPfunc(pfunc, PartitionFunction.Status.Estimated, targetEpsilon, "4.370068e+04" /* e=0.001 */);
	}
	@Test public void protein1Cpu() { calcProteinPfunc(Parallelism.make(1, 0, 0)); }
	@Test public void protein2Cpus() { calcProteinPfunc(Parallelism.make(2, 0, 0)); }
	@Test public void protein1GpuStream() { calcProteinPfunc(Parallelism.make(1, 1, 1)); }
	@Test public void protein4GpuStreams() { calcProteinPfunc(Parallelism.make(2, 1, 4)); }

	public void calcLigandPfunc(Parallelism parallelism) {

		SimpleConfSpace confSpace = new SimpleConfSpace.Builder()
			.addStrand(ligand)
			.build();

		double targetEpsilon = 0.05;
		PartitionFunction pfunc = calcPfunc(confSpace, parallelism, targetEpsilon);
		assertPfunc(pfunc, PartitionFunction.Status.Estimated, targetEpsilon, "4.467797e+30" /* e=0.001 */);
	}
	@Test public void ligand1Cpu() { calcLigandPfunc(Parallelism.make(1, 0, 0)); }
	@Test public void ligand2Cpus() { calcLigandPfunc(Parallelism.make(2, 0, 0)); }
	@Test public void ligand1GpuStream() { calcLigandPfunc(Parallelism.make(1, 1, 1)); }
	@Test public void ligand4GpuStreams() { calcLigandPfunc(Parallelism.make(2, 1, 4)); }

	public void calcComplexPfunc(Parallelism parallelism) {

		// NOTE: to match the old code precisely, we need to match the old conf space exactly too
		// which means we need to keep the same design position order (CCD is of course sensitive to this)
		// and also add the extra residues in the PDB file that aren't in the strands
		SimpleConfSpace confSpace = new SimpleConfSpace.Builder()
			.addStrand(new Strand.Builder(mol).setTemplateLibrary(templateLib).setResidues(153, 154).build())
			.addStrand(ligand)
			.addStrand(new Strand.Builder(mol).setTemplateLibrary(templateLib).setResidues(195, 241).build())
			.addStrand(new Strand.Builder(mol).setTemplateLibrary(templateLib).setResidues(638, 647).build())
			.addStrand(protein)
			.build();

		double targetEpsilon = 0.8;
		PartitionFunction pfunc = calcPfunc(confSpace, parallelism, targetEpsilon);
		assertPfunc(pfunc, PartitionFunction.Status.Estimated, targetEpsilon, "3.5213742379e+54" /* e=0.05 */);
	}
	@Test public void complex1Cpu() { calcComplexPfunc(Parallelism.make(1, 0, 0)); }
	@Test public void complex2Cpus() { calcComplexPfunc(Parallelism.make(2, 0, 0)); }
	@Test public void complex1GpuStream() { calcComplexPfunc(Parallelism.make(1, 1, 1)); }
	@Test public void complex4GpuStreams() { calcComplexPfunc(Parallelism.make(2, 1, 4)); }

	public static void assertPfunc(PartitionFunction pfunc, PartitionFunction.Status status, double targetEpsilon, String approxQstar) {
		assertThat(pfunc.getStatus(), is(status));
		assertThat(pfunc.getValues().getEffectiveEpsilon(), lessThanOrEqualTo(targetEpsilon));
		double qbound = new BigDecimal(approxQstar).doubleValue()*(1.0 - targetEpsilon);
		assertThat(pfunc.getValues().qstar.doubleValue(), greaterThanOrEqualTo(qbound));
	}
}
