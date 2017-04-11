package edu.duke.cs.osprey.energy;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.junit.BeforeClass;
import org.junit.Test;

import edu.duke.cs.osprey.TestBase;
import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.ConfSearch.EnergiedConf;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.confspace.ParametricMolecule;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.minimization.CCDMinimizer;
import edu.duke.cs.osprey.minimization.Minimizer;
import edu.duke.cs.osprey.minimization.MoleculeObjectiveFunction;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.structure.PDBIO;

public class TestMinimizingEnergyCalculators extends TestBase {
	
	private static double[] ExpectedEnergies = {
		-89.40966969379559,     -89.10792031499590,     -89.80959784194750,     -88.63999143591965,
		-89.12813398453119,     -89.50404412354328,     -88.39619842057273,     -88.88944810225355,
		-88.91539256575278,     -88.37401748234860,     -88.72521745747015,     -88.95852827537131,
		-88.56492542985796,     -89.13542390896990,     -88.39342805723787,     -88.61512935924303
	};
	
	private static SimpleConfSpace confSpace;
	private static ForcefieldParams ffparams;
	private static List<ScoredConf> confs;
	
	@BeforeClass
	public static void beforeClass() {
		
		// get a conf space
		Strand strand = new Strand.Builder(PDBIO.readFile("examples/1CC8.python/1CC8.ss.pdb")).build();
		strand.flexibility.get(39).setLibraryRotamers("ALA").setContinuous();
		strand.flexibility.get(43).setLibraryRotamers("ALA").setContinuous();
		strand.flexibility.get(40).setLibraryRotamers().setContinuous();
		strand.flexibility.get(41).setLibraryRotamers().setContinuous();
		strand.flexibility.get(42).setLibraryRotamers().setContinuous();
		strand.flexibility.get(44).setLibraryRotamers().setContinuous();
		strand.flexibility.get(45).setLibraryRotamers().setContinuous();
		confSpace = new SimpleConfSpace.Builder().addStrand(strand).build();
		
		ffparams = new ForcefieldParams();
		
		// get an energy matrix
		MinimizingFragmentEnergyCalculator ecalc = new MinimizingFragmentEnergyCalculator.Builder(confSpace, ffparams)
			.setParallelism(Parallelism.makeCpu(2))
			.build();
		EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(confSpace, ecalc)
			.build()
			.calcEnergyMatrix();
		
		// get the low-energy confs
		ConfSearch tree = new ConfAStarTree.Builder(emat, confSpace).build();
		final int numConfs = 16;
		confs = new ArrayList<>();
		for (int i=0; i<numConfs; i++) {
			confs.add(tree.nextConf());
		}
	}
	
	@Test
	public void defaults() {
		assertEnergies(new MinimizingFragmentEnergyCalculator.Builder(confSpace, ffparams));
	}
	
	@Test
	public void cpuOneThread() {
		assertEnergies(new MinimizingFragmentEnergyCalculator.Builder(confSpace, ffparams)
			.setType(MinimizingFragmentEnergyCalculator.Type.Cpu)
			.setParallelism(Parallelism.makeCpu(1)));
	}
	
	@Test
	public void cpuTwoThreads() {
		assertEnergies(new MinimizingFragmentEnergyCalculator.Builder(confSpace, ffparams)
			.setType(MinimizingFragmentEnergyCalculator.Type.Cpu)
			.setParallelism(Parallelism.makeCpu(2)));
	}
	
	@Test
	public void openclOneStream() {
		assertEnergies(new MinimizingFragmentEnergyCalculator.Builder(confSpace, ffparams)
			.setType(MinimizingFragmentEnergyCalculator.Type.OpenCL)
			.setParallelism(Parallelism.makeGpu(1, 1)));
	}
	
	@Test
	public void openclTwoStreams() {
		assertEnergies(new MinimizingFragmentEnergyCalculator.Builder(confSpace, ffparams)
			.setType(MinimizingFragmentEnergyCalculator.Type.OpenCL)
			.setParallelism(Parallelism.makeGpu(1, 2)));
	}
	
	@Test
	public void cudaOneStream() {
		assertEnergies(new MinimizingFragmentEnergyCalculator.Builder(confSpace, ffparams)
			.setType(MinimizingFragmentEnergyCalculator.Type.Cuda)
			.setParallelism(Parallelism.makeGpu(1, 1)));
	}
	
	@Test
	public void cudaTwoStreams() {
		assertEnergies(new MinimizingFragmentEnergyCalculator.Builder(confSpace, ffparams)
			.setType(MinimizingFragmentEnergyCalculator.Type.Cuda)
			.setParallelism(Parallelism.makeGpu(1, 2)));
	}
	
	@Test
	public void cudaCcdOneStream() {
		assertEnergies(new MinimizingFragmentEnergyCalculator.Builder(confSpace, ffparams)
			.setType(MinimizingFragmentEnergyCalculator.Type.CudaCCD)
			.setParallelism(Parallelism.makeGpu(1, 1)));
	}
	
	@Test
	public void cudaCcdTwoStreams() {
		assertEnergies(new MinimizingFragmentEnergyCalculator.Builder(confSpace, ffparams)
			.setType(MinimizingFragmentEnergyCalculator.Type.CudaCCD)
			.setParallelism(Parallelism.makeGpu(1, 2)));
	}
	
	private void assertEnergies(MinimizingFragmentEnergyCalculator.Builder builder) {
		
		for (EnergyPartition epart : Arrays.asList(new EnergyPartition.Traditional(), new EnergyPartition.AllOnPairs())) {
		
			MinimizingConfEnergyCalculator ecalc = new MinimizingConfEnergyCalculator.Builder(builder.build())
				.setEnergyPartition(epart)
				.build();
			assertEnergies(ecalc);
			ecalc.cleanup();
		}
	}
	
	private void assertEnergies(MinimizingConfEnergyCalculator ecalc) {
		
		List<EnergiedConf> econfs = ecalc.calcAllEnergies(confs);
		
		assertThat(econfs.size(), is(ExpectedEnergies.length));
		
		for (int i=0; i<econfs.size(); i++) {
			assertEnergy(ExpectedEnergies[i], econfs.get(i).getEnergy(), "conf " + i);
		}
	}
	
	private void assertEnergy(double exp, double obs, String desc) {
		
		// for minimized energies, lower observed energy is ok,
		// since improvements to minimizers over time could give us lower energies
		if (obs < exp) {
			return;
		}
		
		assertThat(desc, obs, isAbsolutely(exp, 1e-9));
	}
	
	public static void main(String[] args) {
		
		beforeClass();
		
		// compute the expected energies
		for (int i=0; i<confs.size(); i++) {
			
			ParametricMolecule pmol = confSpace.makeMolecule(confs.get(i).getAssignments());
			Minimizer.Result result = new CCDMinimizer(new MoleculeObjectiveFunction(
				pmol,
				confSpace.makeBounds(confs.get(i).getAssignments()),
				new EnergyFunctionGenerator(new ForcefieldParams()).interactionEnergy(FFInterGen.makeFullConf(confSpace, pmol.mol))	
			), false).minimize();
			
			// print the expected energy
			if (i > 0) {
				System.out.print(",");
			}
			System.out.print(i % 4 == 0 ? "\n" : " ");
			System.out.print(String.format("%22.14f", result.energy));
		}
	}
	
	/* NOTE: These test are designed to fail by running the machine out of memory!
		that can be pretty nasty if you're not expecting it,
		so only enable them to find memory allocation/cleanup issues
		and you should monitor process memory usage very carefully while these tests are running
	*/
	
	//@Test
	public void openclMemoryStress() {
		stressMemory(new MinimizingFragmentEnergyCalculator.Builder(confSpace, ffparams)
			.setType(MinimizingFragmentEnergyCalculator.Type.OpenCL)
			.setParallelism(Parallelism.makeGpu(1, 4)));
	}
	
	//@Test
	public void cudaMemoryStress() {
		stressMemory(new MinimizingFragmentEnergyCalculator.Builder(confSpace, ffparams)
			.setType(MinimizingFragmentEnergyCalculator.Type.Cuda)
			.setParallelism(Parallelism.makeGpu(1, 4)));
	}
	
	//@Test
	public void cudaCCDMemoryStress() {
		stressMemory(new MinimizingFragmentEnergyCalculator.Builder(confSpace, ffparams)
			.setType(MinimizingFragmentEnergyCalculator.Type.CudaCCD)
			.setParallelism(Parallelism.makeGpu(1, 4)));
	}

	private void stressMemory(MinimizingFragmentEnergyCalculator.Builder builder) {
		for (int i=0; i<100; i++) {
			MinimizingConfEnergyCalculator ecalc = new MinimizingConfEnergyCalculator.Builder(builder.build()).build();
			ConfEnergyCalculator.Async.Listener listener = (EnergiedConf econf) -> {};
			ecalc.calcEnergyAsync(confs.get(0), listener);
			ecalc.calcEnergyAsync(confs.get(1), listener);
			ecalc.calcEnergyAsync(confs.get(2), listener);
			ecalc.calcEnergyAsync(confs.get(3), listener);
			ecalc.getTasks().waitForFinish();
			ecalc.cleanup();
		}
	}
}
