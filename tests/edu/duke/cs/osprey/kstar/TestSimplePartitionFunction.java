package edu.duke.cs.osprey.kstar;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import java.math.BigDecimal;

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
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.kstar.pfunc.SimplePartitionFunction;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.tools.Factory;

public class TestSimplePartitionFunction extends TestBase {
	
	@BeforeClass
	public static void before() {
		initDefaultEnvironment();
	}
	
	public static SimplePartitionFunction makePfunc(SearchProblem search) {
		
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
		
		KSSearchProblem search = TestPartitionFunction.makeSearch(KSTermini.PROTEIN, "648", "654", "649 650 651 654"); 
		SimplePartitionFunction pfunc = makePfunc(search);

		// compute it
		final double targetEpsilon = 0.05;
		pfunc.init(targetEpsilon);
		pfunc.compute();
	
		// check the answer
		assertThat(pfunc.getStatus(), is(PartitionFunction.Status.Estimated));
		assertThat(pfunc.getValues().getEffectiveEpsilon(), lessThanOrEqualTo(targetEpsilon));
		assertThat(pfunc.getValues().qstar, isRelatively(new BigDecimal("4.3704590631e+04"), targetEpsilon));
	}
	
	@Test
	public void testLigand() {
		
		KSSearchProblem search = TestPartitionFunction.makeSearch(KSTermini.LIGAND, "155", "194", "156 172 192 193");
		SimplePartitionFunction pfunc = makePfunc(search);

		// compute it
		final double targetEpsilon = 0.05;
		pfunc.init(targetEpsilon);
		pfunc.compute();
	
		// check the answer
		assertThat(pfunc.getStatus(), is(PartitionFunction.Status.Estimated));
		assertThat(pfunc.getValues().getEffectiveEpsilon(), lessThanOrEqualTo(targetEpsilon));
		assertThat(pfunc.getValues().qstar, isRelatively(new BigDecimal("4.4699772362e+30"), targetEpsilon));
	}
	
	@Test
	public void testComplex() {
		
		KSSearchProblem search = TestPartitionFunction.makeSearch(KSTermini.COMPLEX, null, null, "649 650 651 654 156 172 192 193");
		SimplePartitionFunction pfunc = makePfunc(search);

		// compute it
		final double targetEpsilon = 0.8;
		pfunc.init(targetEpsilon);
		pfunc.compute();
	
		// check the answer
		assertThat(pfunc.getStatus(), is(PartitionFunction.Status.Estimated));
		assertThat(pfunc.getValues().getEffectiveEpsilon(), lessThanOrEqualTo(targetEpsilon));
		assertThat(pfunc.getValues().qstar, isRelatively(new BigDecimal("3.5178662402e+54"), targetEpsilon));
	}
}
