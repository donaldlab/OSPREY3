package edu.duke.cs.osprey.kstar;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Arrays;

import org.junit.BeforeClass;
import org.junit.Test;

import edu.duke.cs.osprey.TestBase;
import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.astar.conf.order.AStarOrder;
import edu.duke.cs.osprey.astar.conf.order.StaticScoreHMeanAStarOrder;
import edu.duke.cs.osprey.astar.conf.scoring.AStarScorer;
import edu.duke.cs.osprey.astar.conf.scoring.MPLPPairwiseHScorer;
import edu.duke.cs.osprey.astar.conf.scoring.PairwiseGScorer;
import edu.duke.cs.osprey.astar.conf.scoring.mplp.NodeUpdater;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.ConfEnergyCalculator;
import edu.duke.cs.osprey.control.ConfSearchFactory;
import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.control.MinimizingEnergyCalculator;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.epic.EPICSettings;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.kstar.pfunc.SimplePartitionFunction;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.tools.Factory;
import edu.duke.cs.osprey.tupexp.LUTESettings;

public class TestSimplePartitionFunction extends TestBase {
	
	@BeforeClass
	public static void before() {
		initDefaultEnvironment();
	}
	
	public static SimplePartitionFunction makePfunc(int strand, String firstResNumber, String lastResNumber, String flexibleResNumbers) {
		
		// create the search problem
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
		SearchProblem search = new SearchProblem(
			"test", "test/2RL0.kstar/2RL0.min.reduce.pdb", 
			resFlex.flexResList, resFlex.allowedAAs, addWt, doMinimize, useEpic, new EPICSettings(), useTupleExpansion, new LUTESettings(),
			new DEEPerSettings(), moveableStrands, freeBBZones, useEllipses, useERef, addResEntropy, addWtRots, termini,
			false, new ArrayList<>()
		);
		
		// calc energy matrix
		search.emat = (EnergyMatrix)search.calcMatrix(SearchProblem.MatrixType.EMAT);
		
		// calc pruning matrix
		// NOTE: don't actually need any pruning, A* is fast enough for this small problem
		final double pruningInterval = 5;
		search.pruneMat = new PruningMatrix(search.confSpace, pruningInterval);
		
		// make the A* tree factory
		ConfSearchFactory confSearchFactory = new ConfSearchFactory() {
			@Override
			public ConfSearch make(EnergyMatrix emat, PruningMatrix pmat) {
				
                AStarScorer gscorer = new PairwiseGScorer(emat);
				AStarScorer hscorer = new MPLPPairwiseHScorer(new NodeUpdater(), emat, 1, 0.0001);
				AStarOrder order = new StaticScoreHMeanAStarOrder();
                RCs rcs = new RCs(pmat);
                
                return new ConfAStarTree(order, gscorer, hscorer, rcs);
			}
		};
		
		// make the energy calculator
		Factory<EnergyFunction,Molecule> efuncs = new Factory<EnergyFunction,Molecule>() {
			@Override
			public EnergyFunction make(Molecule mol) {
				return EnvironmentVars.curEFcnGenerator.fullConfEnergy(search.confSpace, search.shellResidues, mol);
			}
		};
		ConfEnergyCalculator.Async ecalc = new MinimizingEnergyCalculator(search, efuncs);
	
		// make the pfunc
		return new SimplePartitionFunction(search.emat, search.pruneMat, confSearchFactory, ecalc);
	}
	
	@Test
	public void testProtein() {
		
		SimplePartitionFunction pfunc = makePfunc(KSTermini.PROTEIN, "648", "654", "649 650 651 654");

		// compute it
		final double targetEpsilon = 0.05;
		pfunc.init(targetEpsilon);
		PartitionFunction.Tools.stepUntilDone(pfunc);
	
		// check the answer
		assertThat(pfunc.getStatus(), is(PartitionFunction.Status.Estimated));
		assertThat(pfunc.getValues().getEffectiveEpsilon(), lessThanOrEqualTo(targetEpsilon));
		assertThat(pfunc.getValues().qstar, isRelatively(new BigDecimal("4.3704590631e+04"), targetEpsilon));
	}
	
	@Test
	public void testLigand() {
		
		SimplePartitionFunction pfunc = makePfunc(KSTermini.LIGAND, "155", "194", "156 172 192 193");

		// compute it
		final double targetEpsilon = 0.05;
		pfunc.init(targetEpsilon);
		PartitionFunction.Tools.stepUntilDone(pfunc);
	
		// check the answer
		assertThat(pfunc.getStatus(), is(PartitionFunction.Status.Estimated));
		assertThat(pfunc.getValues().getEffectiveEpsilon(), lessThanOrEqualTo(targetEpsilon));
		assertThat(pfunc.getValues().qstar, isRelatively(new BigDecimal("4.4699772362e+30"), targetEpsilon));
	}
	
	@Test
	public void testComplex() {
		
		SimplePartitionFunction pfunc = makePfunc(KSTermini.COMPLEX, null, null, "649 650 651 654 156 172 192 193");

		// compute it
		final double targetEpsilon = 0.8;
		pfunc.init(targetEpsilon);
		PartitionFunction.Tools.stepUntilDone(pfunc);
	
		// check the answer
		assertThat(pfunc.getStatus(), is(PartitionFunction.Status.Estimated));
		assertThat(pfunc.getValues().getEffectiveEpsilon(), lessThanOrEqualTo(targetEpsilon));
		assertThat(pfunc.getValues().qstar, isRelatively(new BigDecimal("3.5178662402e+54"), targetEpsilon));
	}
}
