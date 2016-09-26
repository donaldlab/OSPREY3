package edu.duke.cs.osprey.kstar;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Arrays;

import org.junit.BeforeClass;
import org.junit.Test;

import edu.duke.cs.osprey.TestBase;
import edu.duke.cs.osprey.confspace.PositionConfSpace;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.epic.EPICSettings;
import edu.duke.cs.osprey.energy.MultiTermEnergyFunction;
import edu.duke.cs.osprey.kstar.pfunc.PFAbstract;
import edu.duke.cs.osprey.kstar.pfunc.PFAbstract.EApproxReached;
import edu.duke.cs.osprey.kstar.pfunc.PFFactory;
import edu.duke.cs.osprey.parallelism.ThreadParallelism;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tupexp.LUTESettings;

public class TestPartitionFunction extends TestBase {
	
	@BeforeClass
	public static void before() {
		initDefaultEnvironment();
		
		// configure parallelism
		ThreadParallelism.setNumThreads(1);
		MultiTermEnergyFunction.setNumThreads(1);
	}
	
	public static KSSearchProblem makeSearch(int strand, String firstResNumber, String lastResNumber, String flexibleResNumbers) {
		
		ResidueFlexibility resFlex = new ResidueFlexibility();
		resFlex.addFlexible(flexibleResNumbers);
		
		return makeSearch(strand, firstResNumber, lastResNumber, resFlex);
	}
	
	public static KSSearchProblem makeSearchSequence(int strand, String firstResNumber, String lastResNumber, String sequence) {
		
		ResidueFlexibility resFlex = new ResidueFlexibility();
		for (String pos : sequence.split(" ")) {
			String[] parts = pos.split("-");
			resFlex.addMutable(parts[1], parts[0]);
		}
		
		return makeSearch(strand, firstResNumber, lastResNumber, resFlex);
	}
	
	public static KSSearchProblem makeSearch(int strand, String firstResNumber, String lastResNumber, ResidueFlexibility resFlex) {
		
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
		KSTermini termini = null;
		if (firstResNumber != null && lastResNumber != null) {
			termini = new KSTermini(strand, resFlex.size(), new ArrayList<>(Arrays.asList(firstResNumber, lastResNumber)));
		}
		KSSearchProblem search = new KSSearchProblem(
			"test", "test/2RL0.kstar/2RL0.min.reduce.pdb", 
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
	
	public static PFAbstract makePfunc(KSSearchProblem search, String pfImpl, int strand) {
		return makePfunc(search, pfImpl, strand, null);
	}
	
	public static PFAbstract makePfunc(KSSearchProblem search, String pfImpl, int strand, String sequence) {
		
		// as far as I can tell, partition functions only use the CFP to get the steric threshold
		// and to configure HOT things
		// so create a CFP here with just defaults only for the pfunc
		KSConfigFileParser cfp = new KSConfigFileParser();
		
		return makePfunc(search, pfImpl, strand, sequence, cfp);
	}
	
	public static PFAbstract makePfunc(KSSearchProblem search, String pfImpl, int strand, String sequence, KSConfigFileParser cfp) {
		
		// build the sequence info
		ArrayList<String> sequenceList = new ArrayList<>();
		ArrayList<Integer> positions = new ArrayList<>();
		
		if (sequence != null) {
			sequenceList.addAll(Arrays.asList(sequence.split(" ")));
		} else {
			for (int i=0; i<search.confSpace.posFlex.size(); i++) {
				PositionConfSpace posConfSpace = search.confSpace.posFlex.get(i);
				sequenceList.add(posConfSpace.RCs.get(0).AAType + "-" + posConfSpace.res.getPDBResNumber());
			}
		}
		
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
		System.out.println(String.format("Q*: %.10e, epsilon: %.8f",
			pfunc.getQStar().doubleValue(),
			pfunc.getEffectiveEpsilon()
		));
	}
	
	private PFAbstract makeAndComputePfunc(int strand, String firstResidueNumber, String lastResidueNumber, String sequenceOrFlexibleResidues, String pfImpl, double targetEpsilon) {
		
		// do we have a sequence? or a list of flexible residues?
		PFAbstract pfunc;
		if (sequenceOrFlexibleResidues.indexOf('-') > 0) {
			
			// yup, it's a sequence
			KSSearchProblem search = makeSearchSequence(strand, firstResidueNumber, lastResidueNumber, sequenceOrFlexibleResidues);
			pfunc = makePfunc(search, pfImpl, strand, sequenceOrFlexibleResidues);
			
		} else {
			
			// it's just flexible residues
			KSSearchProblem search = makeSearch(strand, firstResidueNumber, lastResidueNumber, sequenceOrFlexibleResidues);
			pfunc = makePfunc(search, pfImpl, strand);
		}
		
		// compute it
		PFAbstract.targetEpsilon = targetEpsilon;
		pfunc.start();
		pfunc.runToCompletion();
		
		return pfunc;
	}
	
	private void testProteinWildType(String pfImpl) {
		
		PFAbstract pfunc = makeAndComputePfunc(KSTermini.PROTEIN, "648", "654", "649 650 651 654", pfImpl, 0.05);
	
		// is this the right partition function?
		assertThat(pfunc.getNumUnPruned().intValueExact(), is(5130));
		assertThat(pfunc.getNumPruned().intValueExact(), is(0));

		// check the answer
		assertPfunc(pfunc, EApproxReached.TRUE, 0.05, "4.3704590631e+04");
	}
	
	private void testLigandWildType(String pfImpl) {
		
		PFAbstract pfunc = makeAndComputePfunc(KSTermini.LIGAND, "155", "194", "156 172 192 193", pfImpl, 0.05);

		// is this the right partition function?
		assertThat(pfunc.getNumUnPruned().intValueExact(), is(21280));
		assertThat(pfunc.getNumPruned().intValueExact(), is(0));
		
		// check the answer
		assertPfunc(pfunc, EApproxReached.TRUE, 0.05, "4.4699772362e+30");
	}
	
	private void testComplexWildType(String pfImpl) {
		
		PFAbstract pfunc = makeAndComputePfunc(KSTermini.COMPLEX, null, null, "649 650 651 654 156 172 192 193", pfImpl, 0.9);

		// is this the right partition function?
		assertThat(pfunc.getNumUnPruned().intValueExact(), is(109166400));
		assertThat(pfunc.getNumPruned().intValueExact(), is(0));
		
		// check the answer
		assertPfunc(pfunc, EApproxReached.TRUE, 0.9, "3.5178662402e+54");
	}
	
	private void assertPfunc(PFAbstract pfunc, EApproxReached epsilonStatus, double targetEpsilon, String qstar) {
		assertThat(pfunc.getEpsilonStatus(), is(epsilonStatus));
		assertThat(pfunc.getEffectiveEpsilon(), lessThanOrEqualTo(targetEpsilon));
		assertThat(pfunc.getQStar(), isRelatively(new BigDecimal(qstar), PFAbstract.targetEpsilon));
	}
	
	private void assertPfunc(PFAbstract pfunc, EApproxReached epsilonStatus) {
		assertThat(pfunc.getEpsilonStatus(), is(epsilonStatus));
	}
	
	
	// wild type tests
	
	@Test
	public void testProteinWildTypeTraditional() {
		testProteinWildType("traditional");
	}
	
	@Test
	public void testProteinWildTypeParallel0() {
		testProteinWildType("parallel0");
	}
	
	@Test
	public void testLigandWildTypeTraditional() {
		
		// use a higher epsilon, since traditional is slow
		PFAbstract pfunc = makeAndComputePfunc(KSTermini.LIGAND, "155", "194", "156 172 192 193", "traditional", 0.9);
		assertPfunc(pfunc, EApproxReached.TRUE, 0.9, "4.4699772362e+30");
	}
	
	@Test
	public void testLigandWildTypeParallel0() {
		testLigandWildType("parallel0");
	}
	
	// traditional can't get anywhere with this pfunc in a reasonable amount of time
	//@Test
	//public void testComplexWildTypeTraditional() {
	//	
	//	// use a higher epsilon, since traditional is slow
	//	PFAbstract pfunc = makeAndComputePfunc(KSTermini.COMPLEX, null, null, "649 650 651 654 156 172 192 193", "traditional", 0.99);
	//	assertPfunc(pfunc, EApproxReached.TRUE, 0.99, "???");
	//}
	
	@Test
	public void testComplexWildTypeParallel0() {
		testComplexWildType("parallel0");
	}
	
	@Test
	public void testProteinWildTypeSimple() {
		testProteinWildType("simple");
	}
	
	@Test
	public void testLigandWildTypeSimple() {
		testLigandWildType("simple");
	}
	
	@Test
	public void testComplexWildTypeSimple() {
		testComplexWildType("simple");
	}
	
	
	// other sequence tests
	
	@Test
	public void testLigandNotPossibleTraditional() {
		PFAbstract pfunc = makeAndComputePfunc(KSTermini.LIGAND, "155", "194", "PHE-156 LYS-172 LEU-192 THR-193", "traditional", 0.95);
		assertPfunc(pfunc, EApproxReached.NOT_POSSIBLE);
	}
	
	@Test
	public void testLigandNotPossibleParallel() {
		PFAbstract pfunc = makeAndComputePfunc(KSTermini.LIGAND, "155", "194", "PHE-156 LYS-172 LEU-192 THR-193", "parallel0", 0.95);
		assertPfunc(pfunc, EApproxReached.NOT_POSSIBLE);
	}
	
	@Test
	public void testLigandNotPossibleSimple() {
		PFAbstract pfunc = makeAndComputePfunc(KSTermini.LIGAND, "155", "194", "PHE-156 LYS-172 LEU-192 THR-193", "simple", 0.95);
		assertPfunc(pfunc, EApproxReached.NOT_POSSIBLE);
	}
	
	@Test
	public void testComplexMutantParallel0() {
		PFAbstract pfunc = makeAndComputePfunc(KSTermini.COMPLEX, null, null, "PHE-649 ASP-650 GLU-651 THR-654 PHE-156 LYS-172 ILE-192 THR-193", "parallel0", 0.9);
		assertPfunc(pfunc, EApproxReached.TRUE, 0.9, "3.519112e+54");
	}
	
	@Test
	public void testComplexMutantSimple() {
		PFAbstract pfunc = makeAndComputePfunc(KSTermini.COMPLEX, null, null, "PHE-649 ASP-650 GLU-651 THR-654 PHE-156 LYS-172 ILE-192 THR-193", "simple", 0.9);
		assertPfunc(pfunc, EApproxReached.TRUE, 0.9, "3.519112e+54");
	}
}
