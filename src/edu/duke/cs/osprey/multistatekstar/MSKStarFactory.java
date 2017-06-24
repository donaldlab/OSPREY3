package edu.duke.cs.osprey.multistatekstar;

import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.ParamSet;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams.SolvationForcefield;
import edu.duke.cs.osprey.gmec.GMECConfEnergyCalculator;
import edu.duke.cs.osprey.gmec.ConfSearchFactory;
import edu.duke.cs.osprey.gmec.MinimizingConfEnergyCalculator;
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

	public static KStarScore makeKStarScore(
			ParamSet msParams,
			int state,
			MSConfigFileParser cfp,
			LMB[] sConstr,
			MSSearchProblem[] searchCont,
			MSSearchProblem[] searchDisc,
			GMECConfEnergyCalculator.Async[] ecalcsCont,
			GMECConfEnergyCalculator.Async[] ecalcsDisc,
			KStarScoreType scoreType
			) {

		ParamSet sParams = cfp.params;
		MSKStarSettings settings = new MSKStarSettings();
		settings.state = state;
		settings.cfp = cfp;
		settings.targetEpsilon = sParams.getDouble("EPSILON");
		settings.numTopConfsToSave = sParams.getInt("NumTopConfsToSave");
		settings.isReportingProgress = msParams.getBool("ISREPORTINGPROGRESS");
		settings.scoreType = scoreType;
		settings.constraints = sConstr;
		int numPartFuncs = sParams.getInt("NUMUBSTATES")+1;
		settings.pfTypes = new PartitionFunctionType[numPartFuncs];
		settings.ecalcs = new GMECConfEnergyCalculator.Async[numPartFuncs];
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
				settings.search[subState].settings.energyLBs = true;
				settings.ecalcs[subState] = ecalcsDisc[subState];
			}
			settings.pfTypes[numPartFuncs-1] = PartitionFunctionType.UpperBound;
			settings.search[numPartFuncs-1] = searchCont[numPartFuncs-1];
			settings.search[numPartFuncs-1].settings.energyLBs = false;
			settings.ecalcs[numPartFuncs-1] = ecalcsCont[numPartFuncs-1];
			settings.isFinal = false;
			settings.isReportingProgress = false;
			settings.numTopConfsToSave = 0;
			return new KStarScoreUpperBound(settings);
			
		case MinimizedLowerBound:
			for(int subState=0;subState<numPartFuncs-1;++subState){
				settings.pfTypes[subState] = PartitionFunctionType.UpperBound;
				settings.search[subState] = searchCont[subState];
				settings.search[subState].settings.energyLBs = false;
				settings.ecalcs[subState] = ecalcsCont[subState];
			}
			settings.pfTypes[numPartFuncs-1] = PartitionFunctionType.Discrete;
			settings.search[numPartFuncs-1] = searchDisc[numPartFuncs-1];
			settings.search[numPartFuncs-1].settings.energyLBs = true;
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
				settings.search[subState].settings.energyLBs = true;
				settings.ecalcs[subState] = ecalcsDisc[subState];
			}
			settings.pfTypes[numPartFuncs-1] = PartitionFunctionType.UpperBound;
			settings.search[numPartFuncs-1] = searchDisc[numPartFuncs-1];
			settings.search[numPartFuncs-1].settings.energyLBs = false;
			settings.ecalcs[numPartFuncs-1] = ecalcsDisc[numPartFuncs-1];
			settings.isFinal = false;
			settings.isReportingProgress = false;
			settings.numTopConfsToSave = 0;
			return new KStarScoreUpperBound(settings);
			
		case DiscreteLowerBound:
			for(int subState=0;subState<numPartFuncs-1;++subState){
				settings.pfTypes[subState] = PartitionFunctionType.UpperBound;
				settings.search[subState] = searchDisc[subState];
				settings.search[subState].settings.energyLBs = false;
				settings.ecalcs[subState] = ecalcsDisc[subState];
			}
			settings.pfTypes[numPartFuncs-1] = PartitionFunctionType.Discrete;
			settings.search[numPartFuncs-1] = searchDisc[numPartFuncs-1];
			settings.search[numPartFuncs-1].settings.energyLBs = true;
			settings.ecalcs[numPartFuncs-1] = ecalcsDisc[numPartFuncs-1];
			settings.isFinal = false;
			settings.isReportingProgress = false;
			settings.numTopConfsToSave = 0;
			return new KStarScoreLowerBound(settings);
			
		default:
			throw new UnsupportedOperationException("ERROR: unsupported K* score type"+settings.scoreType);
		}
	}

	public static MinimizingConfEnergyCalculator makeEnergyCalculator(
			MSConfigFileParser cfp,
			SearchProblem multiSeqSearch,
			Parallelism parallelism
			) {
		// make the conf energy calculator
		return MinimizingConfEnergyCalculator.make(
				makeDefaultFFParams(cfp.params),
				multiSeqSearch, 
				parallelism 
				);
	}

	public static PartitionFunction makePartitionFunction(
			PartitionFunctionType type, 
			EnergyMatrix emat,
			PruningMatrix pmat,
			PruningMatrix invmat,
			ConfSearchFactory confSearchFactory,
			GMECConfEnergyCalculator.Async ecalc
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
		ForcefieldParams ffparams = new ForcefieldParams(sParams.getValue("forcefield"));
		ffparams.distDepDielect = sParams.getBool("distDepDielect");
		ffparams.dielectric = sParams.getDouble("dielectConst");
		ffparams.vdwMultiplier = sParams.getDouble("vdwMult");
		ffparams.solvationForcefield = sParams.getBool("DoSolvationE") ? SolvationForcefield.EEF1 : null;
		ffparams.solvScale = sParams.getDouble("SolvScale");
		ffparams.hElect = sParams.getBool("HElect");
		ffparams.hVDW = sParams.getBool("HVDW");
		return ffparams;
	}

}
