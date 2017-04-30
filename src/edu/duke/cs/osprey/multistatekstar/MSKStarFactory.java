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
public class MSKStarFactory {
	
	public static KStarScoreType getKStarScoreType(ParamSet sParams) {
		boolean doMinimize = sParams.getBool("DOMINIMIZE");
		String scoreType = sParams.getValue("SCORETYPE", "FINAL").toLowerCase();
		switch(scoreType) {
		case "final": return doMinimize ? KStarScoreType.Minimized : KStarScoreType.Discrete;
		case "pairwisebound": return KStarScoreType.PairWiseMinimized;
		case "discretelowerbound": return KStarScoreType.DiscreteLowerBound;
		case "discreteupperbound": return KStarScoreType.DiscreteUpperBound;
		case "minimizedlowerbound": return KStarScoreType.MinimizedLowerBound;
		case "minimizedupperbound": return KStarScoreType.MinimizedUpperBound;
		default:
			throw new UnsupportedOperationException("ERROR: unsupported K* score type "+scoreType);
		}
	}
	
	public static KStarScore makeKStarScore(MSKStarSettings settings, PartitionFunction[] pfs) {
		switch(settings.scoreType) {
		case Minimized: return new KStarScoreMinimized(settings, pfs);
		case PairWiseMinimized: return new KStarScoreDiscrete(settings, pfs);
		case MinimizedUpperBound: return new KStarScoreUpperBound(settings, pfs);
		case MinimizedLowerBound: return new KStarScoreLowerBound(settings, pfs);
		case Discrete: return new KStarScoreDiscrete(settings, pfs);
		case DiscreteUpperBound: return new KStarScoreUpperBound(settings, pfs);
		case DiscreteLowerBound: return new KStarScoreLowerBound(settings, pfs);
		default: throw new UnsupportedOperationException("ERROR: unsupported K* score type"+settings.scoreType);
		}
	}
	
	public static KStarScore makeKStarScore(
			ParamSet msParams,
			int state,
			MSConfigFileParser cfp,
			LMB[] sConstr,
			MSSearchProblem[] searchCont,
			MSSearchProblem[] searchDisc,
			ConfEnergyCalculator.Async[] ecalcsCont,
			ConfEnergyCalculator.Async[] ecalcsDisc,
			KStarScoreType scoreType
			) {

		ParamSet sParams = cfp.getParams();
		MSKStarSettings settings = new MSKStarSettings();
		settings.state = state;
		settings.cfp = cfp;
		settings.targetEpsilon = sParams.getDouble("EPSILON");
		settings.numTopConfsToSave = sParams.getInt("NumTopConfsToSave");
		settings.isReportingProgress = msParams.getBool("ISREPORTINGPROGRESS");
		settings.scoreType = scoreType;
		settings.constraints = sConstr;
		int numPartFuncs = sParams.getInt("NUMOFSTRANDS")+1;
		settings.pfTypes = new PartitionFunctionType[numPartFuncs];
		settings.ecalcs = new ConfEnergyCalculator.Async[numPartFuncs];
		settings.search = new MSSearchProblem[numPartFuncs];

		switch(settings.scoreType) {

		case Minimized:
			for(int subState=0;subState<numPartFuncs;++subState){
				settings.pfTypes[subState] = PartitionFunctionType.Minimized;
				settings.search[subState] = searchCont[subState];
				settings.ecalcs[subState] = ecalcsCont[subState];
			}
			settings.isFinal = true;
			return new KStarScoreMinimized(settings);
			
		case PairWiseMinimized:
			for(int subState=0;subState<numPartFuncs;++subState){
				settings.pfTypes[subState] = PartitionFunctionType.Discrete;
				settings.search[subState] = searchCont[subState];
				settings.ecalcs[subState] = ecalcsCont[subState];
			}
			settings.isFinal = true;
			settings.numTopConfsToSave = 0;
			return new KStarScoreDiscrete(settings);
			
		case MinimizedUpperBound:
			for(int subState=0;subState<numPartFuncs-1;++subState){
				settings.pfTypes[subState] = PartitionFunctionType.Discrete;
				settings.search[subState] = searchDisc[subState];
				settings.search[subState].settings.energyLBs = false;
				settings.ecalcs[subState] = ecalcsDisc[subState];
			}
			settings.pfTypes[numPartFuncs-1] = PartitionFunctionType.UpperBound;
			settings.search[numPartFuncs-1] = searchCont[numPartFuncs-1];
			settings.search[numPartFuncs-1].settings.energyLBs = true;
			settings.ecalcs[numPartFuncs-1] = ecalcsCont[numPartFuncs-1];
			settings.isFinal = false;
			settings.isReportingProgress = false;
			settings.numTopConfsToSave = 0;
			return new KStarScoreUpperBound(settings);
			
		case MinimizedLowerBound:
			for(int subState=0;subState<numPartFuncs-1;++subState){
				settings.pfTypes[subState] = PartitionFunctionType.UpperBound;
				settings.search[subState] = searchCont[subState];
				settings.search[subState].settings.energyLBs = true;
				settings.ecalcs[subState] = ecalcsCont[subState];
			}
			settings.pfTypes[numPartFuncs-1] = PartitionFunctionType.Discrete;
			settings.search[numPartFuncs-1] = searchDisc[numPartFuncs-1];
			settings.search[numPartFuncs-1].settings.energyLBs = false;
			settings.ecalcs[numPartFuncs-1] = ecalcsDisc[numPartFuncs-1];
			settings.isFinal = false;
			settings.isReportingProgress = false;
			settings.numTopConfsToSave = 0;
			return new KStarScoreLowerBound(settings);
			
		case Discrete:
			for(int subState=0;subState<numPartFuncs;++subState){
				settings.pfTypes[subState] = PartitionFunctionType.Discrete;
				settings.search[subState] = searchDisc[subState];
				settings.ecalcs[subState] = ecalcsDisc[subState];
			}
			settings.isFinal = true;
			return new KStarScoreDiscrete(settings);
			
		case DiscreteUpperBound:
			for(int subState=0;subState<numPartFuncs-1;++subState){
				settings.pfTypes[subState] = PartitionFunctionType.Discrete;
				settings.search[subState] = searchDisc[subState];
				settings.search[subState].settings.energyLBs = false;
				settings.ecalcs[subState] = ecalcsDisc[subState];
			}
			settings.pfTypes[numPartFuncs-1] = PartitionFunctionType.UpperBound;
			settings.search[numPartFuncs-1] = searchDisc[numPartFuncs-1];
			settings.search[numPartFuncs-1].settings.energyLBs = true;
			settings.ecalcs[numPartFuncs-1] = ecalcsDisc[numPartFuncs-1];
			settings.isFinal = false;
			settings.isReportingProgress = false;
			settings.numTopConfsToSave = 0;
			return new KStarScoreUpperBound(settings);
			
		case DiscreteLowerBound:
			for(int subState=0;subState<numPartFuncs-1;++subState){
				settings.pfTypes[subState] = PartitionFunctionType.UpperBound;
				settings.search[subState] = searchDisc[subState];
				settings.search[subState].settings.energyLBs = true;
				settings.ecalcs[subState] = ecalcsDisc[subState];
			}
			settings.pfTypes[numPartFuncs-1] = PartitionFunctionType.Discrete;
			settings.search[numPartFuncs-1] = searchDisc[numPartFuncs-1];
			settings.search[numPartFuncs-1].settings.energyLBs = false;
			settings.ecalcs[numPartFuncs-1] = ecalcsDisc[numPartFuncs-1];
			settings.isFinal = false;
			settings.isReportingProgress = false;
			settings.numTopConfsToSave = 0;
			return new KStarScoreLowerBound(settings);
			
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
				parallelism 
				);
		return ecalc;
	}

	public static PartitionFunction makePartitionFunction(
			PartitionFunctionType type, 
			EnergyMatrix emat,
			PruningMatrix pmat,
			PruningMatrix invmat,
			ConfSearchFactory confSearchFactory,
			ConfEnergyCalculator.Async ecalc
			) {
		switch(type) {
		case Minimized:
			return new PartitionFunctionMinimized(emat, pmat, invmat, confSearchFactory, ecalc);
		case Discrete:
			return new PartitionFunctionDiscrete(emat, pmat, invmat, confSearchFactory, ecalc);
		case UpperBound:
			return new PartitionFunctionDiscreteUppperBound(emat, pmat, invmat, confSearchFactory, ecalc);
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
