package edu.duke.cs.osprey.energy;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import java.util.ArrayList;
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

public class TestMinimizingEnergyCalculator extends TestBase {
	
	private static double[] ExpectedEnergies = {
		-89.40966798362470,     -89.10800825146424,     -89.80970644639089,     -88.64018070810017,
		-89.12862025945904,     -89.50382978119062,     -88.39560468584334,     -88.88990492205976,
		-88.91540644584246,     -88.37401274130363,     -88.72430424227448,     -88.95887096966123,
		-88.56866361565410,     -89.13542402245633,     -88.39445813793533,     -88.61514337032392
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
		MinimizingEnergyCalculator ecalc = new MinimizingEnergyCalculator.Builder(confSpace, ffparams)
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
		assertEnergies(new MinimizingEnergyCalculator.Builder(confSpace, ffparams));
	}
	
	@Test
	public void cpuOneThread() {
		assertEnergies(new MinimizingEnergyCalculator.Builder(confSpace, ffparams)
			.setType(MinimizingEnergyCalculator.Type.Cpu)
			.setParallelism(Parallelism.makeCpu(1)));
	}
	
	@Test
	public void cpuTwoThreads() {
		assertEnergies(new MinimizingEnergyCalculator.Builder(confSpace, ffparams)
			.setType(MinimizingEnergyCalculator.Type.Cpu)
			.setParallelism(Parallelism.makeCpu(2)));
	}
	
	@Test
	public void openclOneStream() {
		assertEnergies(new MinimizingEnergyCalculator.Builder(confSpace, ffparams)
			.setType(MinimizingEnergyCalculator.Type.OpenCL)
			.setParallelism(Parallelism.makeGpu(1, 1)));
	}
	
	@Test
	public void openclTwoStreams() {
		assertEnergies(new MinimizingEnergyCalculator.Builder(confSpace, ffparams)
			.setType(MinimizingEnergyCalculator.Type.OpenCL)
			.setParallelism(Parallelism.makeGpu(1, 2)));
	}
	
	@Test
	public void cudaOneStream() {
		assertEnergies(new MinimizingEnergyCalculator.Builder(confSpace, ffparams)
			.setType(MinimizingEnergyCalculator.Type.Cuda)
			.setParallelism(Parallelism.makeGpu(1, 1)));
	}
	
	@Test
	public void cudaTwoStreams() {
		assertEnergies(new MinimizingEnergyCalculator.Builder(confSpace, ffparams)
			.setType(MinimizingEnergyCalculator.Type.Cuda)
			.setParallelism(Parallelism.makeGpu(1, 2)));
	}
	
	@Test
	public void cudaCcdOneStream() {
		assertEnergies(new MinimizingEnergyCalculator.Builder(confSpace, ffparams)
			.setType(MinimizingEnergyCalculator.Type.CudaCCD)
			.setParallelism(Parallelism.makeGpu(1, 1)));
	}
	
	@Test
	public void cudaCcdTwoStreams() {
		assertEnergies(new MinimizingEnergyCalculator.Builder(confSpace, ffparams)
			.setType(MinimizingEnergyCalculator.Type.CudaCCD)
			.setParallelism(Parallelism.makeGpu(1, 2)));
	}
	
	private void assertEnergies(MinimizingEnergyCalculator.Builder builder) {
		
		MinimizingEnergyCalculator minimizer = builder.build();
		List<EnergiedConf> econfs = minimizer.calcAllEnergies(confs);
		minimizer.cleanup();
		
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
		stressMemory(new MinimizingEnergyCalculator.Builder(confSpace, ffparams)
			.setType(MinimizingEnergyCalculator.Type.OpenCL)
			.setParallelism(Parallelism.makeGpu(1, 4)));
	}
	
	//@Test
	public void cudaMemoryStress() {
		stressMemory(new MinimizingEnergyCalculator.Builder(confSpace, ffparams)
			.setType(MinimizingEnergyCalculator.Type.Cuda)
			.setParallelism(Parallelism.makeGpu(1, 4)));
	}
	
	//@Test
	public void cudaCCDMemoryStress() {
		stressMemory(new MinimizingEnergyCalculator.Builder(confSpace, ffparams)
			.setType(MinimizingEnergyCalculator.Type.CudaCCD)
			.setParallelism(Parallelism.makeGpu(1, 4)));
	}

	private void stressMemory(MinimizingEnergyCalculator.Builder builder) {
		for (int i=0; i<100; i++) {
			MinimizingEnergyCalculator minimizer = builder.build();
			ConfEnergyCalculator.Async.Listener listener = (EnergiedConf econf) -> {};
			minimizer.calcEnergyAsync(confs.get(0), listener);
			minimizer.calcEnergyAsync(confs.get(1), listener);
			minimizer.calcEnergyAsync(confs.get(2), listener);
			minimizer.calcEnergyAsync(confs.get(3), listener);
			minimizer.getTasks().waitForFinish();
			minimizer.cleanup();
		}
	}
}
