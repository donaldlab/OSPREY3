package edu.duke.cs.osprey.multistatekstar;

import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.ArrayList;

import edu.duke.cs.osprey.control.ConfEnergyCalculator;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.pruning.PruningMatrix;

public class KStarRatio {

	public MultiStateKStarSettings settings;
	public PartitionFunction[] partitionFunctions;
	public ConfEnergyCalculator.Async[] ecalcs;
	public int numStates;

	public KStarRatio(MultiStateKStarSettings settings) {
		this.settings = settings;
		numStates = settings.sps.length;
		ecalcs = new ConfEnergyCalculator.Async[numStates];
		partitionFunctions = new PartitionFunction[numStates];
		
		for(int state=0;state<numStates;++state) {
			//prune matrices
			settings.sps[state].prunePmat(settings.sps[state], settings.pruningWindow, settings.stericThreshold);
			//create ecalcs
			ecalcs[state] = MultiStateKStarSettings.makeEnergyCalculator(settings.cfp, 
					settings.sps[state]);
			//create partition functions
			partitionFunctions[state] = MultiStateKStarSettings.makePartitionFunction(settings.cfp, 
					settings.sps[state].emat, settings.sps[state].pruneMat, ecalcs[state]);
			partitionFunctions[state].init(settings.targetEpsilon);
		}
	}

	public BigDecimal getKStarRatio() {
		BigDecimal ans = BigDecimal.ONE; int state;
		for(state=0;state<partitionFunctions.length-1;++state)
			ans = ans.multiply(partitionFunctions[state].getValues().qstar);
		ans = partitionFunctions[state].getValues().qstar.divide(ans, RoundingMode.HALF_UP);
		return ans;
	}

	/**
	 * compute until maxNumConfs conformations have been processed
	 * @param maxNumConfs
	 */
	public void compute(int maxNumConfs) {
		for(int state=0;state<partitionFunctions.length;++state){
			compute(state, maxNumConfs);
		}
	}
	
	/**
	 * compute until a conf score boltzmann weight of minbound has been processed.
	 * this is used in the second phase to process confs from p*
	 */
	private BigDecimal phase2(int state, BigDecimal pStar) {
		// we have p* / q* = epsilon1 > target epsilon
		// we want p1* / q* <= target epsilon
		// therefore, p1* <= q* x target epsilon
		// we take p1* as our new value of p* and shift 
		// the pairwise lower bound probability mass 
		// of p* - p1* to q*.
		// this is accomplished by enumerating confs in p*
		// until BoltzmannE(sum_scoreWeights) >= p*-p1*
		
		PartitionFunction pf = partitionFunctions[state];
		BigDecimal qstar = pf.getValues().qstar;
		BigDecimal pstar = pf.getValues().pstar;
		BigDecimal targetScoreWeights = pstar.subtract(BigDecimal.valueOf(settings.targetEpsilon).multiply(qstar));
		
		ConfEnergyCalculator.Async ecalc = MultiStateKStarSettings.makeEnergyCalculator(settings.cfp, 
				settings.sps[state]);
		
		PruningMatrix inVMat = ((QPruningMatrix)settings.sps[state].pruneMat).invert();
		PartitionFunction phase2PF = MultiStateKStarSettings.makePartitionFunction(settings.cfp, 
				settings.sps[state].emat, inVMat, ecalc);
		phase2PF.init(0.0);
		((ParallelConfPartitionFunction2)phase2PF).compute(targetScoreWeights);
		
		ecalc.cleanup();
		
		return phase2PF.getValues().qstar;
	}

	private void compute(int state, int maxNumConfs) {
		PartitionFunction pf = partitionFunctions[state];
		pf.compute(maxNumConfs);
		
		ecalcs[state].cleanup();
		ecalcs[state] = null;
		
		//no more q conformations, and we have not reached epsilon
		if(pf.getValues().qprime.compareTo(BigDecimal.ZERO)==0 && 
				pf.getValues().getEffectiveEpsilon() > settings.targetEpsilon) {
			BigDecimal phase2Qstar = phase2(state, pf.getValues().pstar);
			pf.getValues().qstar = pf.getValues().qstar.add(phase2Qstar);
		}
	}

	private ArrayList<LMV> getLMVsForStateOnly(int state) {
		if(settings.constraints==null) return null;
		ArrayList<LMV> ans = new ArrayList<>();

		for(int l=0;l<settings.constraints.length;++l){
			BigDecimal[] coeffs = settings.constraints[l].coeffs;
			if(coeffs[state].compareTo(BigDecimal.ZERO)==0) continue;
			for(int c=0;c<coeffs.length;++c){
				if(coeffs[state].compareTo(BigDecimal.ZERO)!=0 && c!= state) break;
			}
			ans.add(settings.constraints[l]);
		}

		return ans;
	}
}
