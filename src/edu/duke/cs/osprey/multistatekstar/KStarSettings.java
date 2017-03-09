package edu.duke.cs.osprey.multistatekstar;

import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.ConfEnergyCalculator;
import edu.duke.cs.osprey.control.ConfSearchFactory;
import edu.duke.cs.osprey.control.MinimizingEnergyCalculator;
import edu.duke.cs.osprey.control.ParamSet;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.pruning.PruningMatrix;

public class KStarSettings {
	
	public enum ScoreType {
	    Continuous,
	    Discrete,
	    DiscreteLowerBound,
	    DiscreteUpperBound;
	}
	
	public boolean isReportingProgress;
	public double targetEpsilon;
	public int state;
	public int numTopConfsToSave;
	public MSConfigFileParser cfp;
	public MSSearchProblem[] search;
	public ScoreType scoreType;
	public LMV[] constraints;
	public ConfEnergyCalculator.Async[] ecalcs;

	public KStarSettings() {}
	
	public static ConfEnergyCalculator.Async makeEnergyCalculator(
			MSConfigFileParser cfp,
			SearchProblem multiSeqSearch
			) {
		// make the conf energy calculator
		ConfEnergyCalculator.Async ecalc = MinimizingEnergyCalculator.make(
				makeDefaultFFParams(cfp.getParams()),
				multiSeqSearch, 
				Parallelism.makeFromConfig(cfp), 
				true
				);
		return ecalc;
	}

	public static ConfSearchFactory makeConfSearchFactory(
			MSSearchProblem singleSeqSearch, 
			MSConfigFileParser cfp
			) {
		ConfSearchFactory confSearchFactory = ConfSearchFactory.Tools.makeFromConfig(singleSeqSearch, cfp);
		return confSearchFactory;
	}

	public static PartitionFunction makePartitionFunction(
			EnergyMatrix emat,
			PruningMatrix pruneMat, 
			ConfSearchFactory confSearchFactory,
			ConfEnergyCalculator.Async ecalc,
			boolean cont
			) {
		if(!cont) return new SerialPartitionFunction(emat, pruneMat, confSearchFactory, ecalc);
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
