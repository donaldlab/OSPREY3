package edu.duke.cs.osprey.kstar;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import java.math.BigDecimal;
import java.util.ArrayList;
import org.junit.BeforeClass;
import org.junit.Test;

import edu.duke.cs.osprey.TestBase;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.epic.EPICSettings;
import edu.duke.cs.osprey.energy.MultiTermEnergyFunction;
import edu.duke.cs.osprey.kstar.pfunc.PFAbstract;
import edu.duke.cs.osprey.kstar.pfunc.PFAbstract.EApproxReached;
import edu.duke.cs.osprey.multistatekstar.ResidueTermini;
import edu.duke.cs.osprey.kstar.pfunc.PFFactory;
import edu.duke.cs.osprey.parallelism.ThreadParallelism;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tupexp.LUTESettings;


/*
 * This is a version of TestPartitionFunction
 * that tests the ability to handle residue numbers with different chain IDs
 * including if the residue numbers overlap
 */

public class TestPartitionFunctionChainID extends TestBase {
	
	private static final int NumThreads = 2;
	private static final int NumGpus = 0;
	
	@BeforeClass
	public static void before() {
		initDefaultEnvironment();
		
		// configure parallelism
		ThreadParallelism.setNumThreadsIfPossible(NumThreads);
		MultiTermEnergyFunction.setNumThreads(1);
	}
	
	public static KSSearchProblem makeSearch(String pdbPath, int strand, String firstResNumber, String lastResNumber, String flexibility) {
		
		// configure residue flexibility
		ResidueFlexibility resFlex = new ResidueFlexibility();
		for (String pos : flexibility.split(" ")) {
			String[] parts = pos.split("-");
			if (parts.length == 1) {
				
				// just a flexible position
				resFlex.addFlexible(parts[0]);
				
			} else {
				
				// mutable position
				resFlex.addMutable(parts[1], parts[0]);
			}
		}
		
		// create the K* search problem
		boolean doMinimize = true;
		boolean addWt = true;
		boolean useEpic = false;
		boolean useTupleExpansion = false;
		boolean useEllipses = false;
		boolean useERef = true;
		boolean addResEntropy = false;
		boolean addWtRots = true;
		ArrayList<String[]> moveableStrands = new ArrayList<String[]>();
		ArrayList<String[]> freeBBZones = new ArrayList<String[]>();
		ResidueTermini termini = null;
		if (firstResNumber != null && lastResNumber != null) {
			termini = new ResidueTermini(strand, firstResNumber, lastResNumber);
		}
		KSSearchProblem search = new KSSearchProblem(
			null, "test", pdbPath,
			resFlex.flexResList, resFlex.allowedAAs, addWt, doMinimize, useEpic, new EPICSettings(), useTupleExpansion, new LUTESettings(),
			new DEEPerSettings(), moveableStrands, freeBBZones, useEllipses, useERef, addResEntropy, addWtRots, termini, false
		);
		
		// calc energy matrix
		search.numEmatThreads = 2;
		search.emat = (EnergyMatrix)search.calcMatrix(SearchProblem.MatrixType.EMAT);
		
		// calc pruning matrix
		// NOTE: don't actually need any pruning, A* is fast enough for this small problem
		final double pruningInterval = 5;
		search.pruneMat = new PruningMatrix(search.confSpace, pruningInterval);
		
		return search;
	}
	
	public static PFAbstract makePfunc(KSSearchProblem search, String pfImpl, int strand, String flexibility) {
		
		// as far as I can tell, partition functions only use the CFP to get the steric threshold
		// and to configure HOT things
		// so create a CFP here with just defaults only for the pfunc
		KSConfigFileParser cfp = new KSConfigFileParser();
		
		// except set the number of threads/gpus for the new minimizers
		cfp.params.setValue("MinimizationThreads", Integer.toString(NumThreads));
		cfp.params.setValue("MinimizationGpus", Integer.toString(NumGpus));
		
		return makePfunc(search, pfImpl, strand, flexibility, cfp);
	}
	
	public static PFAbstract makePfunc(KSSearchProblem search, String pfImpl, int strand, String flexibility, KSConfigFileParser cfp) {
		
		// build the sequence info
		ArrayList<String> sequenceList = new ArrayList<>();
		ArrayList<Integer> positions = new ArrayList<>();
		for (String pos : flexibility.split(" ")) {
			if (pos.indexOf('-') > 0) {
				
				// mutable position
				sequenceList.add(pos);
				
			} else {
				
				// flexible position, get wt amino acid
				int i = sequenceList.size();
				String aa = search.confSpace.posFlex.get(i).RCs.get(0).AAType;
				sequenceList.add(aa + "-" + pos);
			}
		}
		
		assert (flexibility.split(" ").length == search.confSpace.posFlex.size());
		for (int i=0; i<search.confSpace.posFlex.size(); i++) {
			positions.add(i);
		}
		
		// get the pfunc for the wild-type protein
		String checkpointName = null;
		String sequenceSearchName = "jeff"; // it's a good name
		PFAbstract pfunc = PFFactory.getPartitionFunction(
			pfImpl,
			strand,
			sequenceList,
			positions,
			checkpointName,
			sequenceSearchName,
			cfp,
			search
		);
		
		// initialize the partition function internal state
		pfunc.setNumUnPruned();
		pfunc.setNumPruned();
		
		return pfunc;
	}
	
	@SuppressWarnings("unused")
	private void printPfunc(PFAbstract pfunc) {
		System.out.println(String.format("sequence: %s, q*: %.10e, epsilon: %.6f",
			KSAbstract.list1D2String(pfunc.getSequence(), " "),
			pfunc.getQStar().doubleValue(),
			pfunc.getEffectiveEpsilon()
		));
	}
	
	private PFAbstract makeAndComputePfunc(String pfImpl, int strand, String firstRes, String lastRes, String flexibility, double targetEpsilon) {
		KSSearchProblem search = makeSearch("examples/2RL0.kstar/2RL0.min.reduce.numoverlap.pdb", strand, firstRes, lastRes, flexibility);
		PFAbstract pfunc = makePfunc(search, pfImpl, strand, flexibility);
		PFAbstract.targetEpsilon = targetEpsilon;
		pfunc.start();
		pfunc.runToCompletion();
		pfunc.cleanup();
		return pfunc;
	}
	
	private PFAbstract makeAndComputeProteinPfunc(String pfImpl, String flexibility, double targetEpsilon) {
		return makeAndComputePfunc(pfImpl, 0, "G148", "G154", flexibility, targetEpsilon);
	}
	
	private PFAbstract makeAndComputeLigandPfunc(String pfImpl, String flexibility, double targetEpsilon) {
		return makeAndComputePfunc(pfImpl, 1, "A155", "A194", flexibility, targetEpsilon);
	}
	
	private PFAbstract makeAndComputeComplexPfunc(String pfImpl, String flexibility, double targetEpsilon) {
		return makeAndComputePfunc(pfImpl, 2, null, null, flexibility, targetEpsilon);
	}
	
	private void testProteinWildType(String pfImpl, double targetEpsilon) {
		
		PFAbstract pfunc = makeAndComputeProteinPfunc(pfImpl, "G149 G150 G151 G154", targetEpsilon);
	
		// is this the right partition function?
		assertThat(pfunc.getNumUnPruned().intValueExact(), is(80));
		assertThat(pfunc.getNumPruned().intValueExact(), is(144));

		// check the answer
		assertPfunc(pfunc, EApproxReached.TRUE, targetEpsilon, "4.3704590631e+04"); // e=0.05
	}
	
	private void testLigandWildType(String pfImpl, double targetEpsilon) {
		
		PFAbstract pfunc = makeAndComputeLigandPfunc(pfImpl, "A156 A172 A192 A193", targetEpsilon);
		
		// is this the right partition function?
		assertThat(pfunc.getNumUnPruned().intValueExact(), is(896));
		assertThat(pfunc.getNumPruned().intValueExact(), is(924));
		
		// check the answer
		assertPfunc(pfunc, EApproxReached.TRUE, targetEpsilon, "4.4699772362e+30"); // e=0.05
	}
	
	private void testComplexWildType(String pfImpl, double targetEpsilon) {
		
		PFAbstract pfunc = makeAndComputeComplexPfunc(pfImpl, "G149 G150 G151 G154 A156 A172 A192 A193", targetEpsilon);
		
		// is this the right partition function?
		assertThat(pfunc.getNumUnPruned().intValueExact(), is(8064));
		assertThat(pfunc.getNumPruned().intValueExact(), is(2286144));
		
		// check the answer
		assertPfunc(pfunc, EApproxReached.TRUE, targetEpsilon, "3.5213742379e+54"); // e=0.05
	}
	
	public static void assertPfunc(PFAbstract pfunc, EApproxReached epsilonStatus, double targetEpsilon, String approxQstar) {
		assertThat(pfunc.getEpsilonStatus(), is(epsilonStatus));
		assertThat(pfunc.getEffectiveEpsilon(), lessThanOrEqualTo(targetEpsilon));
		double qbound = new BigDecimal(approxQstar).doubleValue()*(1.0 - targetEpsilon);
		assertThat(pfunc.getQStar().doubleValue(), greaterThanOrEqualTo(qbound));
	}
	
	public static void assertPfunc(PFAbstract pfunc, EApproxReached epsilonStatus) {
		assertThat(pfunc.getEpsilonStatus(), is(epsilonStatus));
	}
	
	// wild type tests
	
	@Test
	public void testProteinWildTypeTraditional() {
		testProteinWildType("traditional", 0.05);
	}
	
	@Test
	public void testProteinWildTypeParallel0() {
		testProteinWildType("parallel0", 0.05);
	}
	
	@Test
	public void testLigandWildTypeTraditional() {
		testLigandWildType("traditional", 0.9); // traditional is slow for this one, use a higher epsilon
	}
	
	@Test
	public void testLigandWildTypeParallel0() {
		testLigandWildType("parallel0", 0.05);
	}
	
	// traditional can't get anywhere with the complex pfunc in a reasonable amount of time
	//@Test
	//public void testComplexWildTypeTraditional() {
	//	testComplexWildType("traditional", 0.99);
	//}
	
	@Test
	public void testComplexWildTypeParallel0() {
		testComplexWildType("parallel0", 0.9);
	}
	
	@Test
	public void testProteinWildTypeParallelConf() {
		testProteinWildType("parallelConf", 0.05);
	}
	
	@Test
	public void testLigandWildTypeParallelConf() {
		testLigandWildType("parallelConf", 0.05);
	}
	
	@Test
	public void testComplexWildTypeParallelConf() {
		testComplexWildType("parallelConf", 0.9);
	}
	
	
	// mutant tests
	
	@Test
	public void testLigandNotPossibleTraditional() {
		PFAbstract pfunc = makeAndComputeLigandPfunc("traditional", "PHE-A156 LYS-A172 LEU-A192 THR-A193", 0.95);
		assertPfunc(pfunc, EApproxReached.NOT_POSSIBLE);
	}
	
	@Test
	public void testLigandNotPossibleParallel() {
		PFAbstract pfunc = makeAndComputeLigandPfunc("parallel0", "PHE-A156 LYS-A172 LEU-A192 THR-A193", 0.95);
		assertPfunc(pfunc, EApproxReached.NOT_POSSIBLE);
	}
	
	@Test
	public void testLigandNotPossibleParallelConf() {
		PFAbstract pfunc = makeAndComputeLigandPfunc("parallelConf", "PHE-A156 LYS-A172 LEU-A192 THR-A193", 0.95);
		assertPfunc(pfunc, EApproxReached.NOT_POSSIBLE);
	}
	
	@Test
	public void testComplexMutant1Parallel0() {
		PFAbstract pfunc = makeAndComputeComplexPfunc("parallel0", "PHE-G149 ASP-G150 GLU-G151 THR-G154 PHE-A156 LYS-A172 ILE-A192 THR-A193", 0.95);
		assertPfunc(pfunc, EApproxReached.TRUE, 0.95, "3.5213742379e+54"); // e=0.05
	}
	
	@Test
	public void testComplexMutant1ParallelConf() {
		PFAbstract pfunc = makeAndComputeComplexPfunc("parallelConf", "PHE-G149 ASP-G150 GLU-G151 THR-G154 PHE-A156 LYS-A172 ILE-A192 THR-A193", 0.95);
		assertPfunc(pfunc, EApproxReached.TRUE, 0.95, "3.5213742379e+54"); // e=0.05
	}
	
	@Test
	public void testComplexMutant2Parallel0() {
		PFAbstract pfunc = makeAndComputeComplexPfunc("parallel0", "PHE-G149 ASP-G150 GLU-G151 THR-G154 PHE-A156 LYS-A172 LEU-A192 THR-A193", 0.95);
		assertPfunc(pfunc, EApproxReached.TRUE, 0.95, "3.2878232508e+12"); // e=0.05
	}
	
	@Test
	public void testComplexMutant2ParallelConf() {
		PFAbstract pfunc = makeAndComputeComplexPfunc("parallelConf", "PHE-G149 ASP-G150 GLU-G151 THR-G154 PHE-A156 LYS-A172 LEU-A192 THR-A193", 0.95);
		assertPfunc(pfunc, EApproxReached.TRUE, 0.95, "3.2878232508e+12"); // e=0.05
	}
}
