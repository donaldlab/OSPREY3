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
import edu.duke.cs.osprey.kstar.pfunc.ParallelPartitionFunction;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.pruning.PruningMatrix;

public class MultiStateKStarSettings {
	
	public boolean isReportingProgress;
	public double targetEpsilon;
	public double pruningWindow;
	public double stericThreshold;
	public MultiStateConfigFileParser cfp;
	public MultiStateSearchProblem[] search;
	public LMV[] constraints;
	public ConfEnergyCalculator.Async[] ecalcs;
	public String[] sequence;
	
	public MultiStateKStarSettings() {}
	
	public String getFormattedSequence() {
		if(sequence==null) return null;
		StringBuilder sb = new StringBuilder();
		for(String aa : sequence) sb.append(aa+" ");
		return sb.toString().trim();
	}
	
	public static ConfEnergyCalculator.Async makeEnergyCalculator(MultiStateConfigFileParser cfp,
			SearchProblem search) {
		// make the conf energy calculator
		ConfEnergyCalculator.Async ecalc = MinimizingEnergyCalculator.make(makeDefaultFFParams(cfp.getParams()),
				search, Parallelism.makeFromConfig(cfp), true);
		return ecalc;
	}

	public static PartitionFunction makePartitionFunction(MultiStateConfigFileParser cfp,
			EnergyMatrix emat, PruningMatrix pruneMat, ConfEnergyCalculator.Async ecalc) {

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

		return new ParallelPartitionFunction(emat, pruneMat, confSearchFactory, ecalc);
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
