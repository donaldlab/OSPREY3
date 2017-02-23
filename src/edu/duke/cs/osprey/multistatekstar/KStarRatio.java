package edu.duke.cs.osprey.multistatekstar;

import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.ArrayList;

import edu.duke.cs.osprey.control.ConfEnergyCalculator;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;

public class KStarRatio {

	public KStarSettings ksSet;
	public PartitionFunction[] partitionFunctions;
	ConfEnergyCalculator.Async[] ecalcs;
	public LMV[] lmvs;

	public KStarRatio(KStarSettings ksSet,
			PartitionFunction[] partitionFunctions,
			ConfEnergyCalculator.Async[] ecalcs,
			LMV[] lmvs) {

			this.ksSet = ksSet;
			this.partitionFunctions = partitionFunctions;
			this.ecalcs = ecalcs;
			this.lmvs = lmvs;
	}

	public BigDecimal getKStarRatio() {
		
		BigDecimal ans = BigDecimal.ONE; int state;
		for(state=0;state<partitionFunctions.length-1;++state)
			ans = ans.multiply(partitionFunctions[state].getValues().qstar);
		ans = partitionFunctions[state].getValues().qstar.divide(ans, RoundingMode.HALF_UP);
		return ans;
		
	}
	
	public void compute() {
		for(int state=0;state<partitionFunctions.length;++state){
			ecalcs[state].cleanup();
		}
	}
	
	private void compute(PartitionFunction pFunc, int maxNumConfs) {
		for(int phase=1;phase<3;++phase) {
			pFunc.compute(maxNumConfs);
		}
	}
	
	private ArrayList<LMV> getLMVsForStateOnly(int state) {
		if(lmvs==null) return null;
		ArrayList<LMV> ans = new ArrayList<>();
		
		for(int l=0;l<lmvs.length;++l){
			BigDecimal[] coeffs = lmvs[l].coeffs;
			if(coeffs[state].compareTo(BigDecimal.ZERO)==0) continue;
			for(int c=0;c<coeffs.length;++c){
				if(coeffs[state].compareTo(BigDecimal.ZERO)!=0 && c!= state) break;
			}
			ans.add(lmvs[l]);
		}
		
		return ans;
	}
}
