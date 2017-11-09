package edu.duke.cs.osprey.minimization;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.junit.BeforeClass;
import org.junit.Test;

import edu.duke.cs.osprey.TestBase;
import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.ConfSearch.EnergiedConf;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.confspace.ParameterizedMoleculeCopy;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.ematrix.SimpleEnergyMatrixCalculator;
import edu.duke.cs.osprey.ematrix.epic.EPICSettings;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.energy.EnergyFunctionGenerator;
import edu.duke.cs.osprey.energy.FFInterGen;
import edu.duke.cs.osprey.energy.MultiTermEnergyFunction;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldInteractions;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams.SolvationForcefield;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tools.Factory;
import edu.duke.cs.osprey.tupexp.LUTESettings;

public class TestMinimization extends TestBase {
	
	private static Map<Boolean,Info> Infos;
	
	private static final double Epsilon = 1e-7;
	
	private static class Info {
		
		public ForcefieldParams ffparams;
		public SearchProblem search;
		public EnergyFunctionGenerator efuncgen;
		public Factory<ForcefieldInteractions,Molecule> intergen;
		public SimpleConfSpace simpleConfSpace;
		public List<ScoredConf> confs;
		public double[] expectedEnergies;
		
		public Info(boolean doSolv, double[] expectedEnergies, int[][] expectedConfs) {
			
			ffparams = makeDefaultFFParams();
			ffparams.solvationForcefield = doSolv ? SolvationForcefield.EEF1 : null;
			
			this.expectedEnergies = expectedEnergies;
			
			efuncgen = new EnergyFunctionGenerator(ffparams);
			EnvironmentVars.curEFcnGenerator = efuncgen;
			
			ResidueFlexibility resFlex = new ResidueFlexibility();
			resFlex.addMutable("39 43", "ALA");
			resFlex.addFlexible("40 41 42 44 45");
			resFlex.sortPositions();
			boolean doMinimize = true;
			boolean addWt = false;
			boolean useEpic = false;
			boolean useTupleExpansion = false;
			boolean useEllipses = false;
			boolean useERef = false;
			boolean addResEntropy = false;
			boolean addWtRots = false;
			ArrayList<String[]> moveableStrands = new ArrayList<String[]>();
			ArrayList<String[]> freeBBZones = new ArrayList<String[]>();
			
			search = new SearchProblem(
				"test", "examples/1CC8/1CC8.ss.pdb", 
				resFlex.flexResList, resFlex.allowedAAs, addWt, doMinimize, useEpic, new EPICSettings(), useTupleExpansion, new LUTESettings(),
				new DEEPerSettings(), moveableStrands, freeBBZones, useEllipses, useERef, addResEntropy, addWtRots, null,
				false, new ArrayList<>()
			);
			
			// make energy function factories
			intergen = (mol) -> FFInterGen.makeFullConf(search.confSpace, search.shellResidues, mol);
			
			// compute the energy matrix and pruning matrix
			search.emat = new SimpleEnergyMatrixCalculator.Cpu(2, ffparams, search.confSpace, search.shellResidues).calcEnergyMatrix();
			search.pruneMat = new PruningMatrix(search.confSpace, 1000);
			
			// prep new-style emat calculation
			Strand strand = new Strand.Builder(PDBIO.readFile("examples/python.GMEC/1CC8.ss.pdb")).build();
			strand.flexibility.get("A39").setLibraryRotamers("ALA").setContinuous();
			strand.flexibility.get("A43").setLibraryRotamers("ALA").setContinuous();
			strand.flexibility.get("A40").setLibraryRotamers(Strand.WildType).setContinuous();
			strand.flexibility.get("A41").setLibraryRotamers(Strand.WildType).setContinuous();
			strand.flexibility.get("A42").setLibraryRotamers(Strand.WildType).setContinuous();
			strand.flexibility.get("A44").setLibraryRotamers(Strand.WildType).setContinuous();
			strand.flexibility.get("A45").setLibraryRotamers(Strand.WildType).setContinuous();
			simpleConfSpace = new SimpleConfSpace.Builder().addStrand(strand).build();
			assertConfSpacesMatch(search.confSpace, simpleConfSpace);

			// build A* tree
			ConfAStarTree tree = new ConfAStarTree.Builder(search.emat, search.pruneMat)
				.setMPLP(new ConfAStarTree.MPLPBuilder()
					.setNumIterations(1)
				).build();
			
			// get the confs
			final int numConfs = 16;
			confs = new ArrayList<>();
			for (int i=0; i<numConfs; i++) {
				confs.add(tree.nextConf());
			}

			// check the confs
			for (int i=0; i<numConfs; i++) {
				assertThat(confs.get(i).getAssignments(), is(expectedConfs[i]));
			}
		}

		@Override
		public String toString() {
			return String.format("solvation: %s", ffparams.solvationForcefield);
		}
	}
	
	@BeforeClass
	public static void before() {
		initDefaultEnvironment();
		
		// for these small problems, more than one thread is actually slower
		MultiTermEnergyFunction.setNumThreads(1);
		
		Infos = new HashMap<>();
		Infos.put(true, new Info(true, new double[] {
			-89.40966969380000,     -89.10792031501593,     -89.80959784194780,     -88.63999143525771,
			-89.12813398454430,     -89.50404412354364,     -88.39619842044286,     -88.88944810225273,
			-88.91539256575295,     -88.37401748235216,     -88.72521745744406,     -88.95852827535202,
			-88.56492542985737,     -89.13542390896987,     -88.39342805735203,     -88.61512935924370
		}, new int[][] {
			new int[] { 0, 3, 0, 22, 0, 4, 1 },
			new int[] { 0, 4, 0, 22, 0, 4, 1 },
			new int[] { 0, 3, 0, 18, 0, 4, 1 },
			new int[] { 0, 3, 5, 22, 0, 4, 1 },
			new int[] { 0, 3, 5, 21, 0, 4, 1 },
			new int[] { 0, 4, 0, 18, 0, 4, 1 },
			new int[] { 0, 4, 5, 22, 0, 4, 1 },
			new int[] { 0, 4, 5, 21, 0, 4, 1 },
			new int[] { 0, 3, 0, 21, 0, 4, 1 },
			new int[] { 0, 2, 0, 22, 0, 4, 1 },
			new int[] { 0, 3, 5, 18, 0, 4, 1 },
			new int[] { 0, 3, 5, 19, 0, 4, 1 },
			new int[] { 0, 3, 5, 17, 0, 4, 1 },
			new int[] { 0, 3, 0, 15, 0, 4, 1 },
			new int[] { 0, 3, 0, 25, 0, 4, 1 },
			new int[] { 0, 4, 0, 21, 0, 4, 1 }
		}));
		Infos.put(false, new Info(false, new double[] {
			-70.47295175891054,     -70.18816916335444,     -70.29217193373283,     -70.48982613389194,
			-69.42023183578743,     -69.91475960627933,     -70.35460308647797,     -69.19541718369999,
			-69.11087623590001,     -69.43654754385510,     -69.48055945486242,     -70.12045267097507,
			-69.42165499907402,     -68.97467845375986,     -69.34083168223796,     -69.74355354866111
		}, new int[][] {
			new int[] { 0, 3, 5, 22, 0, 4, 1 },
			new int[] { 0, 4, 5, 22, 0, 4, 1 },
			new int[] { 0, 3, 0, 22, 0, 4, 1 },
			new int[] { 0, 3, 5, 21, 0, 4, 1 },
			new int[] { 0, 2, 5, 22, 0, 4, 1 },
			new int[] { 0, 4, 0, 22, 0, 4, 1 },
			new int[] { 0, 4, 5, 21, 0, 4, 1 },
			new int[] { 0, 2, 0, 22, 0, 4, 1 },
			new int[] { 0, 3, 5, 25, 0, 4, 1 },
			new int[] { 0, 2, 5, 21, 0, 4, 1 },
			new int[] { 0, 3, 5, 18, 0, 4, 1 },
			new int[] { 0, 3, 0, 18, 0, 4, 1 },
			new int[] { 0, 3, 5, 17, 0, 4, 1 },
			new int[] { 0, 4, 5, 25, 0, 4, 1 },
			new int[] { 0, 4, 5, 18, 0, 4, 1 },
			new int[] { 0, 4, 0, 18, 0, 4, 1 }
		}));
	}
	
	public static void main(String[] args) {
		
		before();
		
		for (boolean doSolv : Arrays.asList(true, false)) {
		
			Info info = Infos.get(doSolv);
			
			System.out.println("Solvation? " + doSolv);
			
			// compute the expected energies
			for (int i=0; i<info.confs.size(); i++) {
				
				// get the objective function
				ParameterizedMoleculeCopy pmol = new ParameterizedMoleculeCopy(info.search.confSpace);
				EnergyFunction efunc = info.efuncgen.interactionEnergy(info.intergen.make(pmol.getCopiedMolecule()));
				RCTuple tuple = new RCTuple(info.confs.get(i).getAssignments());
				MoleculeModifierAndScorer mof = new MoleculeModifierAndScorer(efunc, info.search.confSpace, tuple, pmol);
				
				// use the original CCD minimizer
				Minimizer.Result result = new CCDMinimizer(mof, false).minimize();
				
				// print the expected energy
				if (i > 0) {
					System.out.print(",");
				}
				System.out.print(i % 4 == 0 ? "\n" : " ");
				System.out.print(String.format("%22.14f", result.energy));
			}
			
			System.out.println();
		}
	}
	
	@Test
	public void testSearchProblemMinimizer() {
		
		for (boolean doSolv : Arrays.asList(true, false)) {
			Info info = Infos.get(doSolv);
			
			for (int i=0; i<info.expectedEnergies.length; i++) {
				double energy = info.search.minimizedEnergy(info.confs.get(i).getAssignments());
				assertThat(energy, isAbsolutely(info.expectedEnergies[i], Epsilon));
			}
		}
	}
	
	@Test
	public void testCpuConfMinimizer1Thread() {
		check((ffparams, intergen, confSpace) -> new CpuConfMinimizer.Builder(ffparams, intergen, confSpace).build());
	}
	
	@Test
	public void testCpuConfMinimizer2Threads() {
		check((ffparams, intergen, confSpace) -> new CpuConfMinimizer.Builder(ffparams, intergen, confSpace).setNumThreads(2).build());
	}
	
	@Test
	public void testCudaConfMinmizer1Stream() {
		check((ffparams, intergen, confSpace) -> new GpuConfMinimizer.Builder(ffparams, intergen, confSpace).setGpuInfo(GpuConfMinimizer.Type.Cuda, 1, 1).build());
	}
	
	@Test
	public void testCudaConfMinmizer2Streams() {
		check((ffparams, intergen, confSpace) -> new GpuConfMinimizer.Builder(ffparams, intergen, confSpace).setGpuInfo(GpuConfMinimizer.Type.Cuda, 1, 2).build());
	}
	
	@Test
	public void testCudaCCDConfMinmizer1Stream() {
		check((ffparams, intergen, confSpace) -> new GpuConfMinimizer.Builder(ffparams, intergen, confSpace).setGpuInfo(GpuConfMinimizer.Type.CudaCCD, 1, 1).build());
	}
	
	@Test
	public void testCudaCCDConfMinmizer2Streams() {
		check((ffparams, intergen, confSpace) -> new GpuConfMinimizer.Builder(ffparams, intergen, confSpace).setGpuInfo(GpuConfMinimizer.Type.CudaCCD, 1, 2).build());
	}
	
	@Test
	public void testOpenCLConfMinmizer1Stream() {
		check((ffparams, intergen, confSpace) -> new GpuConfMinimizer.Builder(ffparams, intergen, confSpace).setGpuInfo(GpuConfMinimizer.Type.OpenCL, 1, 1).build());
	}
	
	@Test
	public void testOpenCLConfMinmizer2Streams() {
		check((ffparams, intergen, confSpace) -> new GpuConfMinimizer.Builder(ffparams, intergen, confSpace).setGpuInfo(GpuConfMinimizer.Type.OpenCL, 1, 2).build());
	}
	
	@Test
	public void testCpuOriginalCCD1Thread() {
		check(EnergyCalculator.Type.CpuOriginalCCD, Parallelism.makeCpu(1));
	}
	@Test
	public void testCpuOriginalCCD2Threads() {
		check(EnergyCalculator.Type.CpuOriginalCCD, Parallelism.makeCpu(2));
	}
	
	@Test
	public void testCpu1Thread() {
		check(EnergyCalculator.Type.Cpu, Parallelism.makeCpu(1));
	}
	@Test
	public void testCpu2Threads() {
		check(EnergyCalculator.Type.Cpu, Parallelism.makeCpu(2));
	}
	
	@Test
	public void testOpenCL1Stream() {
		check(EnergyCalculator.Type.OpenCL, Parallelism.make(4, 1, 1));
	}
	@Test
	public void testOpenCL2Streams() {
		check(EnergyCalculator.Type.OpenCL, Parallelism.make(4, 1, 2));
	}
	
	@Test
	public void testCuda1Stream() {
		check(EnergyCalculator.Type.Cuda, Parallelism.make(4, 1, 1));
	}
	@Test
	public void testCuda2Streams() {
		check(EnergyCalculator.Type.Cuda, Parallelism.make(4, 1, 2));
	}
	
	@Test
	public void testCudaCCD1Stream() {
		check(EnergyCalculator.Type.CudaCCD, Parallelism.make(4, 1, 1));
	}
	@Test
	public void testCudaCCD2Streams() {
		check(EnergyCalculator.Type.CudaCCD, Parallelism.make(4, 1, 2));
	}
	
	@Test
	public void testResidueCuda1Stream() {
		check(EnergyCalculator.Type.ResidueCuda, Parallelism.make(4, 1, 1));
	}
	@Test
	public void testResidueCuda2Streams() {
		check(EnergyCalculator.Type.ResidueCuda, Parallelism.make(4, 1, 2));
	}
	
	@Test
	public void testResidueCudaCCD1Stream() {
		check(EnergyCalculator.Type.ResidueCudaCCD, Parallelism.make(4, 1, 1));
	}
	@Test
	public void testResidueCudaCCD2Streams() {
		check(EnergyCalculator.Type.ResidueCudaCCD, Parallelism.make(4, 1, 2));
	}
	
	private static interface MinimizerFactory {
		ConfMinimizer make(ForcefieldParams ffparams, Factory<ForcefieldInteractions,Molecule> intergen, ConfSpace confSpace);
	}
	
	private void check(MinimizerFactory factory) {
		
		for (boolean doSolv : Arrays.asList(true, false)) {
			
			Info info = Infos.get(doSolv);
			ConfMinimizer minimizer = factory.make(info.ffparams, info.intergen, info.search.confSpace);
			List<EnergiedConf> econfs = minimizer.minimize(info.confs);
			minimizer.cleanup();
			
			checkConfs(info, econfs);
		}
	}
	
	private void check(EnergyCalculator.Type type, Parallelism parallelism) {
		
		for (boolean doSolv : Arrays.asList(true, false)) {
			
			Info info = Infos.get(doSolv);
			new EnergyCalculator.Builder(info.simpleConfSpace, info.ffparams)
				.setType(type)
				.setParallelism(parallelism)
				.use((ecalc) -> {
					
					ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(info.simpleConfSpace, ecalc).build();
					checkConfs(info, confEcalc.calcAllEnergies(info.confs));
				});
		}
	}
	
	private void checkConfs(Info info, List<EnergiedConf> econfs) {
		
		assertThat(info.toString(), econfs.size(), is(info.confs.size()));
		
		for (int i=0; i<info.confs.size(); i++) {
			
			ScoredConf conf = info.confs.get(i);
			EnergiedConf econf = econfs.get(i);
			
			assertThat(info.toString(), econf.getAssignments(), is(conf.getAssignments()));
			assertThat(info.toString(), econf.getScore(), is(conf.getScore()));
			
			// penalize large errors, but not lower energies
			double absErr = econf.getEnergy() - info.expectedEnergies[i];
			assertThat(info.toString(), absErr, lessThanOrEqualTo(Epsilon));
		}
	}
}