package edu.duke.cs.osprey.energy;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import java.lang.reflect.Field;
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
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.gpu.cuda.GpuStreamPool;
import edu.duke.cs.osprey.minimization.CCDMinimizer;
import edu.duke.cs.osprey.minimization.Minimizer;
import edu.duke.cs.osprey.minimization.MoleculeObjectiveFunction;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.parallelism.TaskExecutor.TaskException;
import edu.duke.cs.osprey.parallelism.TaskExecutor.TaskListener;
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
		Strand strand = new Strand.Builder(PDBIO.readFile("examples/python.GMEC/1CC8.ss.pdb")).build();
		strand.flexibility.get("A39").setLibraryRotamers("ALA").setContinuous();
		strand.flexibility.get("A43").setLibraryRotamers("ALA").setContinuous();
		strand.flexibility.get("A40").setLibraryRotamers(Strand.WildType).setContinuous();
		strand.flexibility.get("A41").setLibraryRotamers(Strand.WildType).setContinuous();
		strand.flexibility.get("A42").setLibraryRotamers(Strand.WildType).setContinuous();
		strand.flexibility.get("A44").setLibraryRotamers(Strand.WildType).setContinuous();
		strand.flexibility.get("A45").setLibraryRotamers(Strand.WildType).setContinuous();
		confSpace = new SimpleConfSpace.Builder().addStrand(strand).build();
		
		ffparams = new ForcefieldParams();
		
		// get an energy matrix
		EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpace, ffparams)
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
	public void energyDefaults() {
		assertEnergies(new EnergyCalculator.Builder(confSpace, ffparams));
	}
	
	@Test
	public void energyCpuOneThread() {
		assertEnergies(new EnergyCalculator.Builder(confSpace, ffparams)
			.setType(EnergyCalculator.Type.Cpu)
			.setParallelism(Parallelism.makeCpu(1)));
	}
	
	@Test
	public void energyCpuTwoThreads() {
		assertEnergies(new EnergyCalculator.Builder(confSpace, ffparams)
			.setType(EnergyCalculator.Type.Cpu)
			.setParallelism(Parallelism.makeCpu(2)));
	}
	
	@Test
	public void energyOpenclOneStream() {
		assertEnergies(new EnergyCalculator.Builder(confSpace, ffparams)
			.setType(EnergyCalculator.Type.OpenCL)
			.setParallelism(Parallelism.make(4, 1, 1)));
	}
	
	@Test
	public void energyOpenclTwoStreams() {
		assertEnergies(new EnergyCalculator.Builder(confSpace, ffparams)
			.setType(EnergyCalculator.Type.OpenCL)
			.setParallelism(Parallelism.make(4, 1, 2)));
	}
	
	@Test
	public void energyCudaOneStream() {
		assertEnergies(new EnergyCalculator.Builder(confSpace, ffparams)
			.setType(EnergyCalculator.Type.Cuda)
			.setParallelism(Parallelism.make(4, 1, 1)));
	}
	
	@Test
	public void energyCudaTwoStreams() {
		assertEnergies(new EnergyCalculator.Builder(confSpace, ffparams)
			.setType(EnergyCalculator.Type.Cuda)
			.setParallelism(Parallelism.make(5, 1, 2)));
	}
	
	@Test
	public void energyCudaCcdOneStream() {
		assertEnergies(new EnergyCalculator.Builder(confSpace, ffparams)
			.setType(EnergyCalculator.Type.CudaCCD)
			.setParallelism(Parallelism.make(4, 1, 1)));
	}
	
	@Test
	public void energyCudaCcdTwoStreams() {
		assertEnergies(new EnergyCalculator.Builder(confSpace, ffparams)
			.setType(EnergyCalculator.Type.CudaCCD)
			.setParallelism(Parallelism.make(4, 1, 2)));
	}
	
	@Test
	public void energyResidueCudaOneStream() {
		assertEnergies(new EnergyCalculator.Builder(confSpace, ffparams)
			.setType(EnergyCalculator.Type.ResidueCuda)
			.setParallelism(Parallelism.make(4, 1, 1)));
	}
	
	@Test
	public void energyResidueCudaTwoStreams() {
		assertEnergies(new EnergyCalculator.Builder(confSpace, ffparams)
			.setType(EnergyCalculator.Type.ResidueCuda)
			.setParallelism(Parallelism.make(4, 1, 2)));
	}
	
	@Test
	public void energyResidueCudaCcdOneStream() {
		assertEnergies(new EnergyCalculator.Builder(confSpace, ffparams)
			.setType(EnergyCalculator.Type.ResidueCudaCCD)
			.setParallelism(Parallelism.make(4, 1, 1)));
	}
	
	@Test
	public void energyResidueCudaCcdTwoStreams() {
		assertEnergies(new EnergyCalculator.Builder(confSpace, ffparams)
			.setType(EnergyCalculator.Type.ResidueCudaCCD)
			.setParallelism(Parallelism.make(4, 1, 2)));
	}
	
	private void assertEnergies(EnergyCalculator.Builder builder) {
		
		for (EnergyPartition epart : Arrays.asList(EnergyPartition.Traditional, EnergyPartition.AllOnPairs)) {
			builder.use((ecalc) -> {
				ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(confSpace, ecalc)
					.setEnergyPartition(epart)
					.build();
				assertEnergies(confEcalc);
			});
		}
	}
	
	private void assertEnergies(ConfEnergyCalculator confEcalc) {
		
		List<EnergiedConf> econfs = confEcalc.calcAllEnergies(confs);
		
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
		stressMemory(new EnergyCalculator.Builder(confSpace, ffparams)
			.setType(EnergyCalculator.Type.OpenCL)
			.setParallelism(Parallelism.make(4, 1, 4)));
	}
	
	//@Test
	public void cudaMemoryStress() {
		stressMemory(new EnergyCalculator.Builder(confSpace, ffparams)
			.setType(EnergyCalculator.Type.Cuda)
			.setParallelism(Parallelism.make(4, 1, 4)));
	}
	
	//@Test
	public void cudaCCDMemoryStress() {
		stressMemory(new EnergyCalculator.Builder(confSpace, ffparams)
			.setType(EnergyCalculator.Type.CudaCCD)
			.setParallelism(Parallelism.make(4, 1, 4)));
	}

	private void stressMemory(EnergyCalculator.Builder builder) {
		for (int i=0; i<100; i++) {
			builder.use((ecalc) -> {
				ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(confSpace, builder.build()).build();
				TaskListener<EnergiedConf> listener = (econf) -> { /* don't care */ };
				confEcalc.calcEnergyAsync(confs.get(0), listener);
				confEcalc.calcEnergyAsync(confs.get(1), listener);
				confEcalc.calcEnergyAsync(confs.get(2), listener);
				confEcalc.calcEnergyAsync(confs.get(3), listener);
				confEcalc.tasks.waitForFinish();
			});
		}
	}
	
	@Test
	public void cleanup() {
		
		new EnergyCalculator.Builder(confSpace, ffparams)
			.setType(EnergyCalculator.Type.Cuda)
			.setParallelism(Parallelism.make(4, 1, 2))
			.use((ecalc) -> {
				
				ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(confSpace, ecalc).build();
				RCTuple frag = new RCTuple(0, 0);
				ResidueInteractions inters = ResInterGen.of(confSpace)
					.addIntra(0)
					.make();
				
				for (int r=0; r<100; r++) {
				
					// calculate energies for a few fragments
					for (int i=0; i<10; i++) {
						confEcalc.calcEnergyAsync(frag, inters, (econf) -> {
							// all is well
						});
					}
					confEcalc.tasks.waitForFinish();
					
					// all the streams should have been returned
					GpuStreamPool pool = getPool(ecalc);
					assertThat(pool.getNumStreamsAvailable(), is(pool.getNumStreams()));
				}
			});
	}
	
	@Test
	public void cleanupWithExceptions() {
		
		new EnergyCalculator.Builder(confSpace, ffparams)
			.setType(EnergyCalculator.Type.Cuda)
			.setParallelism(Parallelism.make(4, 1, 2))
			.use((ecalc) -> {
				
				ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(confSpace, ecalc).build();
				
				for (int r=0; r<100; r++) {
					
					try {
						
						// calculate energies for a few fragments
						for (int i=0; i<10; i++) {
							confEcalc.calcEnergyAsync((RCTuple)null, null, (econf) -> {
								fail("we shouldn't make it to the listener");
							});
						}
						confEcalc.tasks.waitForFinish();
						
					} catch (TaskException ex) {
						
						// all the streams should have been returned
						GpuStreamPool pool = getPool(ecalc);
						assertThat(pool.getNumStreamsAvailable(), is(pool.getNumStreams()));
						continue;
					}
					
					fail("should have had a TaskException");
				}
			});
	}
	
	@Test
	public void cleanupWithExceptionsInListener() {
		
		new EnergyCalculator.Builder(confSpace, ffparams)
			.setType(EnergyCalculator.Type.Cuda)
			.setParallelism(Parallelism.make(4, 1, 2))
			.use((ecalc) -> {
				
				ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(confSpace, ecalc).build();
				RCTuple frag = new RCTuple(0, 0);
				ResidueInteractions inters = ResInterGen.of(confSpace)
					.addIntra(0)
					.make();
		
				for (int r=0; r<100; r++) {
					
					try {
						
						// calculate energies for a few fragments
						for (int i=0; i<10; i++) {
							confEcalc.calcEnergyAsync(frag, inters, (econf) -> {
								
								// throw an exception in the listener
								throw new Error("oops, a Bad Thing happened");
							});
						}
						confEcalc.tasks.waitForFinish();
						
					} catch (TaskException ex) {
						
						// all the streams should have been returned
						GpuStreamPool pool = getPool(ecalc);
						assertThat(pool.getNumStreamsAvailable(), is(pool.getNumStreams()));
						continue;
					}
					
					fail("should have had a TaskException");
				}
			});
	}
	
	private GpuStreamPool getPool(EnergyCalculator ecalc) {
		try {
			
			// HACKHACK: get the GpuStreamPool from the ecalc via reflection
			Field fieldPool = ecalc.context.getClass().getDeclaredField("pool");
			fieldPool.setAccessible(true);
			return (GpuStreamPool)fieldPool.get(ecalc.context);
			
		} catch (Exception ex) {
			throw new Error(ex);
		}
	}

	public static void checkFinalPose(EnergyCalculator.Builder builder) {
		builder.use((ecalc) -> {

			ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(confSpace, ecalc)
				.build();

			// calc the minimized energy the usual way
			ScoredConf conf = confs.get(0);
			RCTuple frag = new RCTuple(conf.getAssignments());
			EnergyCalculator.EnergiedParametricMolecule epmol = confEcalc.calcEnergy(frag);

			// to actually check the protein pose, calculate the energy of the final pose
			EnergyFunction efunc = ecalc.makeEnergyFunction(epmol);
			double poseEnergy = efunc.getEnergy();
			EnergyFunction.Tools.cleanIfNeeded(efunc);
			assertThat(poseEnergy, isAbsolutely(epmol.energy, 1e-12));
		});
	}

	@Test
	public void poseDefaults() {
		checkFinalPose(new EnergyCalculator.Builder(confSpace, ffparams));
	}

	@Test
	public void poseCpu() {
		checkFinalPose(new EnergyCalculator.Builder(confSpace, ffparams)
			.setType(EnergyCalculator.Type.Cpu)
			.setParallelism(Parallelism.makeCpu(1)));
	}

	@Test
	public void poseOpencl() {
		checkFinalPose(new EnergyCalculator.Builder(confSpace, ffparams)
			.setType(EnergyCalculator.Type.OpenCL)
			.setParallelism(Parallelism.make(4, 1, 1)));
	}

	@Test
	public void poseCuda() {
		checkFinalPose(new EnergyCalculator.Builder(confSpace, ffparams)
			.setType(EnergyCalculator.Type.Cuda)
			.setParallelism(Parallelism.make(4, 1, 1)));
	}

	@Test
	public void poseCudaCcd() {
		checkFinalPose(new EnergyCalculator.Builder(confSpace, ffparams)
			.setType(EnergyCalculator.Type.CudaCCD)
			.setParallelism(Parallelism.make(4, 1, 1)));
	}

	@Test
	public void poseResidueCuda() {
		checkFinalPose(new EnergyCalculator.Builder(confSpace, ffparams)
			.setType(EnergyCalculator.Type.ResidueCuda)
			.setParallelism(Parallelism.make(4, 1, 1)));
	}

	@Test
	public void poseResidueCudaCcd() {
		checkFinalPose(new EnergyCalculator.Builder(confSpace, ffparams)
			.setType(EnergyCalculator.Type.ResidueCudaCCD)
			.setParallelism(Parallelism.make(4, 1, 1)));
	}
}
