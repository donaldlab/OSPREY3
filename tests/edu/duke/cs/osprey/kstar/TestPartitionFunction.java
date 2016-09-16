package edu.duke.cs.osprey.kstar;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Arrays;

import org.junit.BeforeClass;
import org.junit.Test;

import edu.duke.cs.osprey.TestBase;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.epic.EPICSettings;
import edu.duke.cs.osprey.kstar.pfunc.PFAbstract;
import edu.duke.cs.osprey.kstar.pfunc.PFAbstract.EApproxReached;
import edu.duke.cs.osprey.kstar.pfunc.impl.PFTraditional;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tupexp.LUTESettings;

public class TestPartitionFunction extends TestBase {
	
	@BeforeClass
	public static void before() {
		initDefaultEnvironment();
	}
	
	@Test
	public void testTraditional() {
		
		// create the K* search problem
		ResidueFlexibility resFlex = new ResidueFlexibility();
		resFlex.addFlexible("649 650 651 654");
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
		KSTermini termini = new KSTermini(KSTermini.PROTEIN, 4, new ArrayList<>(Arrays.asList("648", "654")));
		KSSearchProblem search = new KSSearchProblem(
			"test", "test/2RL0.kstar/2RL0.min.reduce.pdb", 
			resFlex.flexResList, resFlex.allowedAAs, addWt, doMinimize, useEpic, new EPICSettings(), useTupleExpansion, new LUTESettings(),
			new DEEPerSettings(), moveableStrands, freeBBZones, useEllipses, useERef, addResEntropy, addWtRots, termini, false
		);
		
		// calc energy matrix
		search.emat = (EnergyMatrix)search.calcMatrix(SearchProblem.MatrixType.EMAT);
		
		// calc pruning matrix
		// NOTE: don't actually need any pruning, A* is fast enough for this small problem
		final double pruningInterval = 5;
		search.pruneMat = new PruningMatrix(search.confSpace, pruningInterval);
		
		// build the sequence info
		ArrayList<String> sequence = new ArrayList<>();
		ArrayList<Integer> positions = new ArrayList<>();
		for (int i=0; i<resFlex.size(); i++) {
			sequence.add(resFlex.allowedAAs.get(i).get(0) + "-" + resFlex.flexResList.get(i));
			positions.add(i);
		}
		
		// as far as I can tell, partition functions only use the CFP to get the steric threshold
		// and to configure HOT things
		// so create a CFP here with just defaults only for the pfunc
		KSConfigFileParser cfp = new KSConfigFileParser();
		
		// get the pfunc for the wild-type protein
		PFAbstract.targetEpsilon = 0.05;
		String checkpointName = null;
		String sequenceSearchName = "jeff"; // it's a good name
		PFAbstract pfunc = new PFTraditional(
		//PFAbstract pfunc = new PFParallel2(
			KSTermini.PROTEIN,
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
		
		// do some checks first to make sure we have a sane partition function
		assertThat(pfunc.getNumUnPruned().intValueExact(), is(5130));
		assertThat(pfunc.getNumPruned().intValueExact(), is(0));
		
		// compute the partition function!
		pfunc.start();
		pfunc.runToCompletion();
	
		printPfunc(pfunc);
		
		// check the answer
		assertThat(pfunc.getEpsilonStatus(), is(EApproxReached.TRUE));
		assertThat(pfunc.getEffectiveEpsilon(), lessThanOrEqualTo(PFAbstract.targetEpsilon));
		assertThat(pfunc.getQStar(), isRelatively(new BigDecimal("4.3700523336e+04"), PFAbstract.targetEpsilon));
	}
	
	private void printPfunc(PFAbstract pfunc) {
		System.out.println(String.format("Q*: %.10e, epsilon: %.8f",
			pfunc.getQStar().doubleValue(),
			pfunc.getEffectiveEpsilon()
		));
	}
}
