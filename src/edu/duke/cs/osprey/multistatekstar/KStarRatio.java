package edu.duke.cs.osprey.multistatekstar;

import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.ArrayList;

import edu.duke.cs.osprey.control.ConfEnergyCalculator;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;

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
					settings.sps[state], ecalcs[state]);
		}
	}

	public BigDecimal getKStarRatio() {

		BigDecimal ans = BigDecimal.ONE; int state;
		for(state=0;state<partitionFunctions.length-1;++state)
			ans = ans.multiply(partitionFunctions[state].getValues().qstar);
		ans = partitionFunctions[state].getValues().qstar.divide(ans, RoundingMode.HALF_UP);
		return ans;

	}

	public void compute(int maxNumConfs) {
		for(int state=0;state<partitionFunctions.length;++state){
			compute(state, maxNumConfs);
		}
	}

	private void compute(int state, int maxNumConfs) {
		PartitionFunction pf = partitionFunctions[state];
		for(int phase=1;phase<3;++phase) {
			//pf.compute(maxNumConfs);
		}
		ecalcs[state].cleanup();
		ecalcs[state] = null;
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
