package edu.duke.cs.osprey.multistatekstar;

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
import edu.duke.cs.osprey.control.MinimizingEnergyCalculator;
import edu.duke.cs.osprey.control.ParamSet;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.kstar.pfunc.ParallelConfPartitionFunction;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.pruning.PruningMatrix;

public class MultiStateKStarSettings {
	
	public double targetEpsilon;
	public double pruningWindow;
	public double stericThreshold;
	public MultiStateConfigFileParser cfp;
	public MultiStateSearchProblem[] sps;
	public LMV[] constraints;
	
	public MultiStateKStarSettings() {}
	
	public static ConfEnergyCalculator.Async makeEnergyCalculator(MultiStateConfigFileParser stateCfp,
			SearchProblem search) {
		// make the conf energy calculator
		ConfEnergyCalculator.Async ecalc = MinimizingEnergyCalculator.make(makeDefaultFFParams(stateCfp.getParams()),
				search, Parallelism.makeFromConfig(stateCfp), true);
		return ecalc;
	}

	public static PartitionFunction makePartitionFunction(MultiStateConfigFileParser stateCfp,
			SearchProblem search, ConfEnergyCalculator.Async ecalc) {

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

		return new ParallelConfPartitionFunction(search.emat, search.pruneMat, confSearchFactory, ecalc);
	}

	public static ForcefieldParams makeDefaultFFParams(ParamSet sParams) {
		// values from default config file
		String forceField = sParams.getValue("forcefield");
		boolean distDepDielect = sParams.getBool("distDepDielect");
		double dielectConst = sParams.getDouble("dielectConst");
		double vdwMult = sParams.getDouble("vdwMult");
		boolean doSolv = sParams.getBool("DoSolvationE");
		double solvScale = sParams.getDouble("SolvScale");
		boolean useHForElectrostatics = sParams.getBool("HElect");
		boolean useHForVdw = sParams.getBool("HVDW");
		return new ForcefieldParams(
				forceField, distDepDielect, dielectConst, vdwMult,
				doSolv, solvScale, useHForElectrostatics, useHForVdw
				);
	}
}
