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
import edu.duke.cs.osprey.parallelism.ThreadParallelism;
import edu.duke.cs.osprey.kstar.pfunc.PFFactory;
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
		
		// create the K* search problem
		ResidueFlexibility resFlex = new ResidueFlexibility();
		resFlex.addFlexible(flexibleResNumbers);
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
		
		// build the sequence info
		ArrayList<String> sequence = new ArrayList<>();
		ArrayList<Integer> positions = new ArrayList<>();
		for (int i=0; i<search.confSpace.posFlex.size(); i++) {
			PositionConfSpace posConfSpace = search.confSpace.posFlex.get(i);
			sequence.add(posConfSpace.RCs.get(0).AAType + "-" + posConfSpace.res.getPDBResNumber());
			positions.add(i);
		}
		
		// as far as I can tell, partition functions only use the CFP to get the steric threshold
		// and to configure HOT things
		// so create a CFP here with just defaults only for the pfunc
		KSConfigFileParser cfp = new KSConfigFileParser();
		
		// get the pfunc for the wild-type protein
		String checkpointName = null;
		String sequenceSearchName = "jeff"; // it's a good name
		PFAbstract pfunc = PFFactory.getPartitionFunction(
			pfImpl,
			strand,
			sequence,
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
	
	private void testProtein(String pfImpl) {
		
		KSSearchProblem search = makeSearch(KSTermini.PROTEIN, "648", "654", "649 650 651 654");
		PFAbstract pfunc = makePfunc(search, pfImpl, KSTermini.PROTEIN);

		// is this the right partition function?
		assertThat(pfunc.getNumUnPruned().intValueExact(), is(5130));
		assertThat(pfunc.getNumPruned().intValueExact(), is(0));
		
		// compute it
		PFAbstract.targetEpsilon = 0.05;
		pfunc.start();
		pfunc.runToCompletion();
	
		// check the answer
		assertThat(pfunc.getEpsilonStatus(), is(EApproxReached.TRUE));
		assertThat(pfunc.getEffectiveEpsilon(), lessThanOrEqualTo(PFAbstract.targetEpsilon));
		assertThat(pfunc.getQStar(), isRelatively(new BigDecimal("4.3704590631e+04"), PFAbstract.targetEpsilon));
	}
	
	private void testLigand(String pfImpl) {
		
		KSSearchProblem search = makeSearch(KSTermini.LIGAND, "155", "194", "156 172 192 193");
		PFAbstract pfunc = makePfunc(search, pfImpl, KSTermini.LIGAND);

		// is this the right partition function?
		assertThat(pfunc.getNumUnPruned().intValueExact(), is(21280));
		assertThat(pfunc.getNumPruned().intValueExact(), is(0));
		
		// compute it
		PFAbstract.targetEpsilon = 0.8;
		pfunc.start();
		pfunc.runToCompletion();
	
		// check the answer
		assertThat(pfunc.getEpsilonStatus(), is(EApproxReached.TRUE));
		assertThat(pfunc.getEffectiveEpsilon(), lessThanOrEqualTo(PFAbstract.targetEpsilon));
		assertThat(pfunc.getQStar(), isRelatively(new BigDecimal("4.4699772362e+30"), PFAbstract.targetEpsilon));
	}
	
	private void testComplex(String pfImpl) {
		
		KSSearchProblem search = makeSearch(KSTermini.COMPLEX, null, null, "649 650 651 654 156 172 192 193");
		PFAbstract pfunc = makePfunc(search, pfImpl, KSTermini.COMPLEX);

		// is this the right partition function?
		assertThat(pfunc.getNumUnPruned().intValueExact(), is(109166400));
		assertThat(pfunc.getNumPruned().intValueExact(), is(0));
		
		// compute it
		PFAbstract.targetEpsilon = 0.90;
		pfunc.start();
		pfunc.runToCompletion();
	
		// check the answer
		assertThat(pfunc.getEpsilonStatus(), is(EApproxReached.TRUE));
		assertThat(pfunc.getEffectiveEpsilon(), lessThanOrEqualTo(PFAbstract.targetEpsilon));
		assertThat(pfunc.getQStar(), isRelatively(new BigDecimal("3.5178662402e+54"), PFAbstract.targetEpsilon));
	}
	
	@Test
	public void testProteinTraditional() {
		testProtein("traditional");
	}
	
	@Test
	public void testProteinParallel0() {
		testProtein("parallel0");
	}
	
	@Test
	public void testProteinParallel1() {
		testProtein("parallel1");
	}
	
	@Test
	public void testProteinParallel2() {
		testProtein("parallel2");
	}
	
	@Test
	public void testLigandTraditional() {
		testLigand("traditional");
	}
	
	// TODO: the parallel implementations don't seem to be working for the ligand 
	@Test
	public void testLigandParallel0() {
		testLigand("parallel0");
	}
	
	@Test
	public void testLigandParallel1() {
		testLigand("parallel1");
	}
	
	@Test
	public void testLigandParallel2() {
		testLigand("parallel2");
	}
	
	@Test
	public void testComplexParallel2() {
		testComplex("parallel2");
	}
}
