package edu.duke.cs.osprey.multistatekstar;

import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.ConfEnergyCalculator;
import edu.duke.cs.osprey.control.ConfSearchFactory;
import edu.duke.cs.osprey.control.MinimizingEnergyCalculator;
import edu.duke.cs.osprey.control.ParamSet;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.multistatekstar.KStarScore.KStarScoreType;
import edu.duke.cs.osprey.multistatekstar.KStarScore.PartitionFunctionType;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.pruning.PruningMatrix;

/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 * 
 */
public class KStarFactory {

	public static KStarScore makeKStarScore(
			ParamSet msParams,
			int state,
			MSConfigFileParser cfp,
			MSSearchProblem[] searchCont,
			MSSearchProblem[] searchDisc,
			ConfEnergyCalculator.Async[] ecalcsCont,
			ConfEnergyCalculator.Async[] ecalcsDisc,
			KStarScoreType scoreType
			) {

		ParamSet sParams = cfp.getParams();
		KStarSettings settings = new KStarSettings();
		settings.state = state;
		settings.cfp = cfp;
		settings.targetEpsilon = sParams.getDouble("EPSILON");
		settings.numTopConfsToSave = sParams.getInt("NumTopConfsToSave");
		settings.isReportingProgress = msParams.getBool("ISREPORTINGPROGRESS");
		settings.scoreType = scoreType;

		//make LMVs
		int numUbConstr = sParams.getInt("NUMUBCONSTR");
		int numPartFuncs = sParams.getInt("NUMUBSTATES")+1;
		settings.constraints = new LMV[numUbConstr];
		for(int constr=0;constr<numUbConstr;constr++)
			settings.constraints[constr] = new LMV(sParams.getValue("UBCONSTR"+constr), numPartFuncs);

		settings.pfTypes = new PartitionFunctionType[numPartFuncs];
		settings.ecalcs = new ConfEnergyCalculator.Async[numPartFuncs];
		settings.search = new MSSearchProblem[numPartFuncs];

		switch(settings.scoreType) {

		case Continuous:
			for(int subState=0;subState<numPartFuncs;++subState){
				settings.pfTypes[subState] = PartitionFunctionType.Continuous;
				settings.search[subState] = searchCont[subState];
				settings.ecalcs[subState] = ecalcsCont[subState];
			}
			return new ContinuousKStarScore(settings);
			
		case Discrete:
			for(int subState=0;subState<numPartFuncs;++subState){
				settings.pfTypes[subState] = PartitionFunctionType.Discrete;
				settings.search[subState] = searchDisc[subState];
				settings.ecalcs[subState] = ecalcsDisc[subState];
			}
			return new DiscreteKStarScore(settings);
			
		case DiscretePairWise:
			for(int subState=0;subState<numPartFuncs;++subState){
				settings.pfTypes[subState] = PartitionFunctionType.Discrete;
				settings.search[subState] = searchCont[subState];
				settings.ecalcs[subState] = ecalcsCont[subState];
			}
			settings.numTopConfsToSave = 0;
			return new DiscreteKStarScore(settings);
			
		case DiscreteUpperBound:
			for(int subState=0;subState<numPartFuncs-1;++subState){
				settings.pfTypes[subState] = PartitionFunctionType.Discrete;
				settings.search[subState] = searchDisc[subState];
				settings.ecalcs[subState] = ecalcsDisc[subState];
			}
			settings.pfTypes[numPartFuncs-1] = PartitionFunctionType.DiscreteUpperBound;
			settings.search[numPartFuncs-1] = searchCont[numPartFuncs-1];
			settings.ecalcs[numPartFuncs-1] = ecalcsCont[numPartFuncs-1];
			settings.numTopConfsToSave = 0;
			return new DiscreteKStarScore(settings);
			
		case DiscreteLowerBound:
			for(int subState=0;subState<numPartFuncs-1;++subState){
				settings.pfTypes[subState] = PartitionFunctionType.DiscreteUpperBound;
				settings.search[subState] = searchCont[subState];
				settings.ecalcs[subState] = ecalcsCont[subState];
			}
			settings.pfTypes[numPartFuncs-1] = PartitionFunctionType.Discrete;
			settings.search[numPartFuncs-1] = searchDisc[numPartFuncs-1];
			settings.ecalcs[numPartFuncs-1] = ecalcsDisc[numPartFuncs-1];
			settings.numTopConfsToSave = 0;
			return new DiscreteKStarScore(settings);
			
		default:
			throw new UnsupportedOperationException("ERROR: unsupported K* score type"+settings.scoreType);
		}
	}

	public static ConfEnergyCalculator.Async makeEnergyCalculator(
			MSConfigFileParser cfp,
			SearchProblem multiSeqSearch,
			Parallelism parallelism
			) {
		// make the conf energy calculator
		ConfEnergyCalculator.Async ecalc = MinimizingEnergyCalculator.make(
				makeDefaultFFParams(cfp.getParams()),
				multiSeqSearch, 
				parallelism, 
				multiSeqSearch.contSCFlex
				);
		return ecalc;
	}

	public static PartitionFunction makePartitionFunction(
			PartitionFunctionType type, 
			EnergyMatrix emat,
			PruningMatrix pruneMat,
			ConfSearchFactory confSearchFactory,
			ConfEnergyCalculator.Async ecalc
			) {
		switch(type) {
		case Continuous:
			return new ContinuousPartitionFunction(emat, pruneMat, confSearchFactory, ecalc);
		case Discrete:
			return new DiscretePartitionFunction(emat, pruneMat, confSearchFactory, ecalc);
		case DiscreteUpperBound:
			return new DiscreteUpperBoundPartitionFunction(emat, pruneMat, confSearchFactory, ecalc);
		default:
			throw new UnsupportedOperationException("ERROR: unsupported partition function type "+type);
		}
	}

	public static ConfSearchFactory makeConfSearchFactory(
			MSSearchProblem singleSeqSearch, 
			MSConfigFileParser cfp
			) {
		ConfSearchFactory confSearchFactory = ConfSearchFactory.Tools.makeFromConfig(singleSeqSearch, cfp);
		return confSearchFactory;
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
