package edu.duke.cs.osprey.multistatekstar;

import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.PriorityQueue;

import edu.duke.cs.osprey.confspace.ConfSearch.EnergiedConf;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.control.ConfSearchFactory;
import edu.duke.cs.osprey.control.GMECFinder;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction.Status;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction.Values;
import edu.duke.cs.osprey.pruning.PruningMatrix;

/**
 * @author Adegoke Ojewole (ao68@duke.edu)
 * 
 */

public class KStarScoreMinimized implements KStarScore {

	public MSKStarSettings settings;
	public PartitionFunctionMinimized[] partitionFunctions;
	public boolean[] initialized;
	public int numStates;
	protected boolean constrSatisfied;

	public KStarScoreMinimized(MSKStarSettings settings) {
		this.settings = settings;
		numStates = settings.search.length;
		partitionFunctions = new PartitionFunctionMinimized[numStates];
		initialized = new boolean[numStates];
		Arrays.fill(partitionFunctions, null);
		Arrays.fill(initialized, false);
		constrSatisfied = true;
	}

	public KStarScoreMinimized(MSKStarSettings settings, PartitionFunction[] other) {
		this(settings);
		for(int state=0;state<numStates;++state) {
			partitionFunctions[state] = (PartitionFunctionMinimized) other[state];
			if(other[state] != null) initialized[state] = true;
		}
	}

	@Override
	public MSKStarSettings getSettings() {
		return settings;
	}

	@Override
	public PartitionFunction getPartitionFunction(int state) {
		return partitionFunctions[state];
	}

	protected BigDecimal getDenom() {
		PartitionFunction pf;
		BigDecimal ans = BigDecimal.ONE.setScale(64, RoundingMode.HALF_UP);
		for(int state=0;state<numStates-1;++state) {
			pf = partitionFunctions[state];
			if(pf==null || pf.getValues().qstar.compareTo(BigDecimal.ZERO)==0)
				return BigDecimal.ZERO.setScale(64, RoundingMode.HALF_UP);
			ans = ans.multiply(pf.getValues().qstar);
		}
		return ans;
	}

	public BigDecimal getScore() {
		BigDecimal den = getDenom();
		if(den.compareTo(BigDecimal.ZERO) == 0) return den;
		PartitionFunction pf = partitionFunctions[numStates-1];
		return pf==null ? BigDecimal.ZERO : pf.getValues().qstar.setScale(64, RoundingMode.HALF_UP).divide(den, RoundingMode.HALF_UP);
	}

	@Override
	public BigDecimal getLowerBoundScore() {
		return getScore();
	}

	@Override
	public BigDecimal getUpperBoundScore() {
		if(isComputed()) return getScore();

		BigDecimal den = getDenom();
		if(den.compareTo(BigDecimal.ZERO) == 0) return den;
		PartitionFunction pf = partitionFunctions[numStates-1];
		if(pf==null) return BigDecimal.ZERO.setScale(64, RoundingMode.HALF_UP);
		BigDecimal num = pf.getValues().qstar.setScale(64, RoundingMode.HALF_UP);
		num = (num.add(pf.getValues().qprime)).add(pf.getValues().pstar);
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

	protected boolean computeMinGMEC() {
		return settings.isFinal && 
			settings.cfp.getParams().getBool("DOMINIMIZE") && 
			(computeMinGMECRatio() || settings.computeGMEC);
	}
	
	protected boolean computeMinGMECRatio() {
		return settings.isFinal &&
			settings.cfp.getParams().getBool("DOMINIMIZE") && 
			settings.cfp.getParams().getBool("COMPUTEGMECRATIO");
	}
	
	protected EnergiedConf findMinGMEC(int state) {
		GMECFinder gmecFinder = new GMECFinder();
		gmecFinder.setLogConfsToConsole(settings.isReportingProgress);
		gmecFinder.isReportingProgress(settings.isReportingProgress);
        gmecFinder.init(settings.cfp, settings.search[state]);
        List<EnergiedConf> energiedConfs = gmecFinder.calcGMEC(settings.search[state]);
        
        if(MSSearchProblem.DEBUG) 
        	settings.search[state].checkPruningMatrix();
        
        return energiedConfs.size() > 0 ? energiedConfs.get(0) : null;
	}
	
	protected boolean init(int state) {		
		if(settings.isReportingProgress) {
			System.out.println();
			System.out.println("state"+state+": "+settings.search[state].settings.getFormattedSequence()+" "+settings.pfTypes[state]);
		}
		
		//find the gmec if we are asked to do so
		//it's important find gmec here, before the prunepmat step, since gmec finding
		//sets the ival to a level required to find the gmec
		EnergiedConf minGMEC = computeMinGMEC() ? findMinGMEC(state) : null;
		
		//first prune the pruning matrix
		boolean doPruning = isFinal() || settings.cfp.getParams().getBool("PRUNEPARTIALSEQCONFS");
		settings.search[state].prunePmat(doPruning, settings.cfp.getParams().getInt("ALGOPTION")>=3);

		//make conf search factory (i.e. A* tree)
		ConfSearchFactory confSearchFactory = MSKStarFactory.makeConfSearchFactory(settings.search[state], settings.cfp);

		//create partition function
		partitionFunctions[state] = (PartitionFunctionMinimized) MSKStarFactory.makePartitionFunction( 
				settings.pfTypes[state],
				settings.search[state].emat, 
				settings.search[state].pruneMat,
				new PruningMatrixInverted(settings.search[state], settings.search[state].pruneMat),
				confSearchFactory,
				settings.ecalcs[state]
				);

		PartitionFunctionMinimized pf = partitionFunctions[state];
		
		pf.setReportProgress(settings.isReportingProgress);

		//skip initialization if we are only interested in the gmec ratio approximation of the k* score
		if(computeMinGMECRatio()) {
			pf.setValues(new Values());
			if(minGMEC != null) pf.getValues().qstar = pf.getBoltzmannCalculator().calc(minGMEC.getEnergy());
			pf.setStatus(Status.Estimated);
		}
		
		//init partition function
		else {
			if(settings.isFinal) {
				pf.setComputeGMECRatio(settings.computeGMECRatio);
				if(minGMEC != null) pf.setMinGMEC(minGMEC);
			}
			
			pf.init(settings.targetEpsilon);
		}
		
		//create priority queue for top confs if requested
		if(settings.search[state].isFullyAssigned() && settings.numTopConfsToSave > 0) {

			pf.topConfs = new PriorityQueue<ScoredConf>(
					settings.numTopConfsToSave, 
					new ConfComparator()
					);

			pf.maxNumTopConfs = settings.numTopConfsToSave;

			pf.setConfListener((ScoredConf conf) -> {
				pf.saveConf(conf);
			});

		}
		
		//save mingmec conf immediately
		if(minGMEC != null && pf.getConfListener() != null) {
			pf.getConfListener().onConf(minGMEC);
			if(computeMinGMECRatio())
				pf.writeTopConfs(settings.state, settings.search[state]);
		}

		return true;
	}

	/**
	 * compute until maxNumConfs conformations have been processed
	 * @param maxNumConfs
	 */
	public void compute(int maxNumConfs) {

		for(int state=0;state<numStates;++state){

			if(!constrSatisfied)//state-specific constraints
				return;

			if(!initialized[state]) {
				initialized[state] = init(state);
			}

			if(partitionFunctions[state].getStatus() != Status.Estimated) {
				compute(state, maxNumConfs);
			}
		}

		//check all constraints now. technically, we should only check constraints
		//that don't pertain to only one state, which we have already checked in compute(state)
		if(settings.isFinal && constrSatisfied) 
			constrSatisfied = checkConstraints();

		if(isComputed()) cleanup();
	}

	/**
	 * compute only unbound states
	 * @param maxNumConfs
	 */
	@Override
	public void computeUnboundStates(int maxNumConfs) {
		for(int state=0;state<numStates-1;++state){

			if(!constrSatisfied)//state-specific constraints
				return;

			if(!initialized[state])
				initialized[state] = init(state);

			if(partitionFunctions[state].getStatus() != Status.Estimated)
				compute(state, maxNumConfs);

			//don't check all constraints, because we are not computing 
			//the bound state partition function
			if(settings.isFinal && constrSatisfied) 
				constrSatisfied = checkConstraints(state);

			partitionFunctions[state].cleanup();
		}
	}

	public void computeBoundState(int maxNumConfs) {
		if(!constrSatisfied)
			return;

		int state = numStates-1;
		if(!initialized[state])
			initialized[state] = init(state);

		if(partitionFunctions[state].getStatus() != Status.Estimated)
			compute(state, maxNumConfs);

		if(partitionFunctions[state].getStatus()==Status.Estimated) {//assumption: unbound states are complete
			if(settings.isFinal && constrSatisfied) 
				constrSatisfied = checkConstraints();
		}

		if(isComputed()) cleanup();
	}

	private void cleanup() {
		for(PartitionFunctionMinimized pf : partitionFunctions) {
			if(pf==null) continue;
			pf.cleanup();
		}
	}

	/**
	 * compute until a conf score boltzmann weight of minbound has been processed.
	 * this is used in the second phase to process confs from p*
	 */
	private PartitionFunction phase2(int state, int maxNumConfs) {
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

		ConfSearchFactory confSearchFactory = MSKStarFactory.makeConfSearchFactory(settings.search[state], settings.cfp);

		PruningMatrix invmat = ((PartitionFunctionMinimized)pf).invmat;

		PartitionFunctionMinimized p2pf = (PartitionFunctionMinimized) MSKStarFactory.makePartitionFunction( 
				settings.pfTypes[state],
				settings.search[state].emat, 
				invmat,
				new PruningMatrixNull(invmat), 
				confSearchFactory,
				settings.ecalcs[state]
				);

		p2pf.setReportProgress(settings.isReportingProgress);

		p2pf.init(targetEpsilon);//enumerating over pstar, energies can be high
		p2pf.getValues().qstar = qstar;//initialize to old qstar
		p2pf.compute(targetScoreWeights, maxNumConfs);

		return p2pf;
	}

	//in the bound state, can override maxNumConfs with value from config
	private boolean overrideMaxNumConfs(int state, int maxNumConfs) {
		return isFinal() && state==numStates-1 && maxNumConfs==Integer.MAX_VALUE;
	}

	protected void compute(int state, int maxNumConfs) {			
		PartitionFunctionMinimized pf = partitionFunctions[state];

		//in the bound state, can override maxNumConfs with value from config
		boolean overrideMaxNumConfs = overrideMaxNumConfs(state, maxNumConfs);
		if(overrideMaxNumConfs) maxNumConfs = settings.cfp.getParams().getInt("MAXNUMCONFS");

		pf.compute(maxNumConfs);

		if(overrideMaxNumConfs && pf.getNumConfsEvaluated()>=maxNumConfs) pf.setStatus(Status.Estimated);
		
		//we are not trying to compute the partition function to completion
		if(pf.getStatus() == Status.Estimating)
			return;

		//no more q conformations, and we have not reached epsilon
		else if(pf.getStatus() == Status.NotEnoughConformations) {
			maxNumConfs = overrideMaxNumConfs ? maxNumConfs-pf.getNumConfsEvaluated() : Integer.MAX_VALUE;
				
			PartitionFunctionMinimized p2pf = (PartitionFunctionMinimized) phase2(state, maxNumConfs);

			pf.getValues().qstar = p2pf.getValues().qstar;
			pf.setNumConfsEvaluated(pf.getNumConfsEvaluated() + p2pf.getNumConfsEvaluated());

			if(settings.search[state].isFullyAssigned() && settings.numTopConfsToSave > 0)
				pf.saveEConfs(p2pf.topConfs);
		}

		pf.setStatus(Status.Estimated);

		if(isFinal()) {//final is a superset of fully defined
			if(constrSatisfied) constrSatisfied = checkConstraints(state);
			if(settings.numTopConfsToSave > 0) pf.writeTopConfs(settings.state, settings.search[state]);
		}
	}

	protected ArrayList<LMB> getLMBsForState(int state, boolean negCoeff) {
		ArrayList<LMB> ans = new ArrayList<>();
		for(LMB constr : getLMBsForState(state)) {
			if(negCoeff && constr.getCoeffs()[state].compareTo(BigDecimal.ZERO)<0) ans.add(constr);
			else if(!negCoeff && constr.getCoeffs()[state].compareTo(BigDecimal.ZERO)>0) ans.add(constr);
		}
		ans.trimToSize();
		return ans;
	}

	/**
	 * returns constraints that ONLY involve the specified state
	 * @param state
	 * @return
	 */
	protected ArrayList<LMB> getLMBsForState(int state) {
		ArrayList<LMB> ans = new ArrayList<>();
		if(settings.constraints==null) return ans;

		for(int l=0;l<settings.constraints.length;++l) {
			BigDecimal[] coeffs = settings.constraints[l].coeffs;
			if(coeffs[state].compareTo(BigDecimal.ZERO)==0) continue;
			boolean addConstr = true;
			for(int c=0;c<coeffs.length;++c) {
				if(c!=state && coeffs[c].compareTo(BigDecimal.ZERO)!=0) {
					addConstr = false;
					break;
				}
			}
			if(addConstr)
				ans.add(settings.constraints[l]);
		}
		ans.trimToSize();
		return ans;
	}

	private boolean checkConstraints(ArrayList<LMB> constraints) {
		for(LMB constr : constraints) {
			BigDecimal[] stateVals = new BigDecimal[numStates];
			for(int s=0;s<numStates;++s){
				PartitionFunctionMinimized pf = partitionFunctions[s];
				stateVals[s] = pf == null ? BigDecimal.ZERO : pf.getValues().qstar;
			}

			//can short circuit computation of k* score if any of the unbound
			//states does not satisfy constraints
			if(constr.eval(stateVals).compareTo(BigDecimal.ZERO) >= 0)
				return false;
		}
		return true;
	}

	/**
	 * see if partition function satisfies either lower or upper bound 
	 * constraints involving this state only
	 * @param state
	 * @param lbConstr: true=lb, false=ub
	 * @return
	 */
	protected boolean checkConstraints(int state, boolean negCoeff) {
		return checkConstraints(getLMBsForState(state, negCoeff));
	}

	/**
	 * see if partition function satisfies constraints involving this state only
	 * @param state
	 * @return
	 */
	protected boolean checkConstraints(int state) {
		return checkConstraints(getLMBsForState(state));
	}

	private boolean checkConstraints() {
		if(!constrSatisfied) return constrSatisfied;
		if(settings.constraints==null) return true;

		BigDecimal[] stateVals = new BigDecimal[numStates];

		for(int c=0;c<settings.constraints.length;++c){
			LMB constr = settings.constraints[c];	
			for(int s=0;s<numStates;++s){
				PartitionFunctionMinimized pf = partitionFunctions[s];
				stateVals[s] = pf == null ? BigDecimal.ZERO : pf.getValues().qstar;
			}

			if(constr.eval(stateVals).compareTo(BigDecimal.ZERO) >= 0)
				return false;
		}
		return true;
	}

	@Override
	public boolean constrSatisfied() {
		return constrSatisfied;
	}

	@Override
	public boolean isFullyProcessed() {
		if(!settings.isFinal) return false;
		return isComputed();
	}

	@Override
	public boolean isComputed() {
		int nulls = 0;
		for(PartitionFunctionMinimized pf : partitionFunctions) {
			if(pf==null) nulls++;
			else if(pf.getStatus() != Status.Estimated) return false;
		}
		//all non-null pfs are estimated; the reason why we skipped a pf must
		//be that a constraint is not satified
		if(nulls>0) {
			if(!constrSatisfied) return true;
			//otherwise, we erroneously skipped a partition function
			else throw new RuntimeException("ERROR: illegally skipped a partition function computation");
		}
		return true;
	}

	@Override
	public boolean isFinal() {
		return settings.isFinal;
	}

	@Override
	public boolean isFullyAssigned() {
		return settings.search[numStates-1].isFullyAssigned();
	}
}
