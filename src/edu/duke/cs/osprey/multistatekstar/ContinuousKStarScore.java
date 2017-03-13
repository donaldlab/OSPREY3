package edu.duke.cs.osprey.multistatekstar;

import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.PriorityQueue;

import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.control.ConfSearchFactory;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction.Status;
import edu.duke.cs.osprey.pruning.PruningMatrix;

/**
 * @author Adegoke Ojewole (ao68@duke.edu)
 * 
 */

public class ContinuousKStarScore implements KStarScore {

	public KStarSettings settings;
	public ContinuousPartitionFunction[] partitionFunctions;
	protected boolean[] initialized;
	public int numStates;
	private boolean constrSatisfied;

	public ContinuousKStarScore(KStarSettings settings) {
		this.settings = settings;
		numStates = settings.search.length;
		partitionFunctions = new ContinuousPartitionFunction[numStates];
		initialized = new boolean[numStates];
		Arrays.fill(partitionFunctions, null);
		Arrays.fill(initialized, false);
		constrSatisfied = true;
	}

	private BigDecimal getDenom() {
		PartitionFunction pf;
		BigDecimal ans = BigDecimal.ONE.setScale(128, RoundingMode.HALF_UP); int state;
		for(state=0;state<partitionFunctions.length-1;++state) {
			pf = partitionFunctions[state];
			if(pf.getValues().qstar.compareTo(BigDecimal.ZERO)==0)
				return BigDecimal.ZERO;
			ans = ans.multiply(pf.getValues().qstar);
		}
		return ans;
	}
	
	public BigDecimal getScore() {
		BigDecimal den = getDenom();
		if(den.compareTo(BigDecimal.ZERO) == 0) return BigDecimal.ZERO;
		PartitionFunction pf = partitionFunctions[numStates-1];
		return pf.getValues().qstar.divide(den, RoundingMode.HALF_UP);
	}
	
	@Override
	public BigDecimal getLowerBoundScore() {
		return getScore();
	}

	@Override
	public BigDecimal getUpperBoundScore() {
		BigDecimal den = getDenom();
		if(den.compareTo(BigDecimal.ZERO) == 0) return BigDecimal.ZERO;
		PartitionFunction pf = partitionFunctions[numStates-1];
		BigDecimal num = pf.getValues().qstar.add(pf.getValues().qprime).add(pf.getValues().pstar);
		return num.divide(den, RoundingMode.HALF_UP);
	}

	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append("Seq: "+settings.search[numStates-1].settings.getFormattedSequence()+", ");
		sb.append(String.format("score: %12e, ", getScore()));
		for(int state=0;state<numStates;++state) {
			BigDecimal qstar = partitionFunctions[state]==null ? BigDecimal.ZERO : 
				partitionFunctions[state].getValues().qstar;
			sb.append(String.format("pf: %2d, q*: %12e, ", state, qstar));
		}
		String ans = sb.toString().trim();
		return ans.substring(0,ans.length()-1);
	}
	
	protected boolean init(int state) {
		//first prune the pruning matrix
		settings.search[state].prunePmat();
		
		//make conf search factory (i.e. A* tree)
		ConfSearchFactory confSearchFactory = KStarFactory.makeConfSearchFactory(settings.search[state], settings.cfp);
		
		//create partition function
		partitionFunctions[state] = (ContinuousPartitionFunction) KStarFactory.makePartitionFunction( 
				settings.pfTypes[state],
				settings.search[state].emat, 
				settings.search[state].pruneMat,
				confSearchFactory,
				settings.ecalcs[state]
				);

		partitionFunctions[state].setReportProgress(settings.isReportingProgress);

		//init partition function
		partitionFunctions[state].init(settings.targetEpsilon);

		//create priority queue for top confs if requested
		if(settings.search[state].isFullyDefined() && settings.numTopConfsToSave > 0) {

			partitionFunctions[state].topConfs = new PriorityQueue<ScoredConf>(
					settings.numTopConfsToSave, 
					new ConfComparator()
					);

			partitionFunctions[state].maxNumTopConfs = settings.numTopConfsToSave;

			final int pfState = state;
			partitionFunctions[state].setConfListener((ScoredConf conf) -> {
				partitionFunctions[pfState].saveConf(conf);
			});
			
		}
		
		return true;
	}

	/**
	 * compute until maxNumConfs conformations have been processed
	 * @param maxNumConfs
	 */
	public void compute(int maxNumConfs) {

		for(int state=0;state<numStates;++state){
			
			if(!constrSatisfied)
				return;

			if(!initialized[state])
				initialized[state] = init(state);
			
			compute(state, maxNumConfs);
		}
		
		//check all constraints now. technically, we should only check constraints
		//that don't pertain to only one state, which we have already checked in compute(state)
		if(constrSatisfied)
			constrSatisfied = checkConstraints();
	}

	/**
	 * compute until a conf score boltzmann weight of minbound has been processed.
	 * this is used in the second phase to process confs from p*
	 */
	private PartitionFunction phase2(int state) {
		// we have p* / q* = epsilon1 > target epsilon
		// we want p1* / q* <= target epsilon
		// therefore, p1* <= q* x target epsilon
		// we take p1* as our new value of p* and shift 
		// the pairwise lower bound probability mass 
		// of p* - p1* to q*.
		// this is accomplished by enumerating confs in p*
		// until BoltzmannE(sum_scoreWeights) >= p* - p1*

		PartitionFunction pf = partitionFunctions[state];
		BigDecimal targetScoreWeights;
		double epsilon = pf.getValues().getEffectiveEpsilon();
		double targetEpsilon = settings.targetEpsilon;
		BigDecimal qstar = pf.getValues().qstar;
		BigDecimal qprime = pf.getValues().qprime;
		BigDecimal pstar = pf.getValues().pstar;

		if(epsilon==1.0) {
			targetScoreWeights = pstar;
		}

		else {
			targetScoreWeights = BigDecimal.valueOf(targetEpsilon/(1.0-targetEpsilon));
			targetScoreWeights = targetScoreWeights.multiply(qstar);
			targetScoreWeights = (pstar.add(qprime)).subtract(targetScoreWeights);
		}

		PruningMatrix invPmat = ((QPruningMatrix)settings.search[state].pruneMat).invert();
		settings.search[state].pruneMat = invPmat;

		ConfSearchFactory confSearchFactory = KStarFactory.makeConfSearchFactory(settings.search[state], settings.cfp);

		ContinuousPartitionFunction p2pf = (ContinuousPartitionFunction) KStarFactory.makePartitionFunction( 
				settings.pfTypes[state],
				settings.search[state].emat, 
				settings.search[state].pruneMat, 
				confSearchFactory,
				settings.ecalcs[state]
				);

		p2pf.init(0.03);
		p2pf.compute(targetScoreWeights);
		return p2pf;
	}

	private void compute(int state, int maxNumConfs) {
		if(settings.isReportingProgress) 
			System.out.println("state"+state+": "+settings.search[state].settings.getFormattedSequence());
		ContinuousPartitionFunction pf = partitionFunctions[state];
		pf.compute(maxNumConfs);	

		//no more q conformations, and we have not reached epsilon
		if(pf.getValues().getEffectiveEpsilon() > settings.targetEpsilon) {
			ContinuousPartitionFunction p2pf = (ContinuousPartitionFunction) phase2(state);
			pf.getValues().qstar = pf.getValues().qstar.add(p2pf.getValues().qstar);
			pf.setStatus(p2pf.getStatus());
			if(settings.search[state].isFullyDefined() && settings.numTopConfsToSave > 0)
				pf.saveEConfs(p2pf.topConfs);
		}
		
		constrSatisfied = checkConstraints(state);
		
		//print top confs
		if(settings.search[state].isFullyDefined() && settings.numTopConfsToSave > 0)
			pf.writeTopConfs(settings.state, settings.search[state]);
	}

	private ArrayList<LMV> getLMVsForState(int state) {
		ArrayList<LMV> ans = new ArrayList<>();
		if(settings.constraints==null) return ans;

		for(int l=0;l<settings.constraints.length;++l){
			BigDecimal[] coeffs = settings.constraints[l].coeffs;
			if(coeffs[state].compareTo(BigDecimal.ZERO)==0) continue;
			boolean addConstr = true;
			for(int c=0;c<coeffs.length;++c){
				if(c!=state && coeffs[c].compareTo(BigDecimal.ZERO)!=0) {
					addConstr = false;
					break;
				}
			}
			if(addConstr)
				ans.add(settings.constraints[l]);
		}

		return ans;
	}
	
	private boolean checkConstraints(int state) {
		
		//see if partition function satisfies constraints
		for(LMV constr : getLMVsForState(state)){
			
			BigDecimal[] stateVals = new BigDecimal[numStates];
			for(int s=0;s<numStates;++s){
				ContinuousPartitionFunction pf = partitionFunctions[s];
				stateVals[s] = pf == null ? BigDecimal.ZERO : pf.getValues().qstar;
			}
			
			//can short circuit computation of k* score if any of the unbound
			//states does not satisfy constraints
			if(constr.eval(stateVals).compareTo(BigDecimal.ZERO) > 0)
				return false;
		}
		
		return true;
	}
	
	private boolean checkConstraints() {
		if(!constrSatisfied) return constrSatisfied;
		if(settings.constraints==null) return true;
		
		BigDecimal[] stateVals = new BigDecimal[numStates];
		
		for(int c=0;c<settings.constraints.length;++c){
			LMV constr = settings.constraints[c];	
			for(int s=0;s<numStates;++s){
				ContinuousPartitionFunction pf = partitionFunctions[s];
				stateVals[s] = pf == null ? BigDecimal.ZERO : pf.getValues().qstar;
			}

			if(constr.eval(stateVals).compareTo(BigDecimal.ZERO) > 0)
				return false;
		}
		return true;
	}

	@Override
	public boolean constrSatisfied() {
		return constrSatisfied;
	}

	@Override
	public boolean computed() {
		for(int state=0;state<numStates;++state) {
			PartitionFunction pf = partitionFunctions[state];
			if(pf.getStatus()==Status.Estimating) return false;
		}
		return true;
	}
}
