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
import edu.duke.cs.osprey.tools.BigDecimalUtil;

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

	@Override
	public BigDecimal getDenominator() {
		PartitionFunction pf;
		BigDecimal ans = BigDecimal.ONE.setScale(64, RoundingMode.HALF_UP);
		for(int state=0;state<numStates-1;++state) {
			pf = partitionFunctions[state];
			if(pf==null || pf.getValues().qstar.compareTo(BigDecimal.ZERO)==0) {
				return BigDecimal.ZERO.setScale(64, RoundingMode.HALF_UP);
			}
			ans = ans.multiply(pf.getValues().qstar);
		}
		return ans;
	}

	@Override
	public BigDecimal getNumerator() {
		BigDecimal ans = partitionFunctions[numStates-1].getValues().qstar;
		if(ans == null) ans = BigDecimal.ZERO.setScale(64, RoundingMode.HALF_UP);
		return ans;
	}

	public BigDecimal toLog10(BigDecimal val) {
		if(val.compareTo(BigDecimal.ZERO)==0)
			return KStarScore.MIN_VALUE;

		return BigDecimalUtil.log10(val);
	}

	public BigDecimal getScore() {
		BigDecimal den = getDenominator();
		if(den.compareTo(BigDecimal.ZERO) == 0) return toLog10(den);
		PartitionFunction pf = partitionFunctions[numStates-1];
		BigDecimal ans = pf==null ? BigDecimal.ZERO : pf.getValues().qstar.setScale(64, RoundingMode.HALF_UP).divide(den, RoundingMode.HALF_UP);
		return toLog10(ans);
	}

	@Override
	public BigDecimal getLowerBoundScore() {
		return getScore();
	}

	@Override
	public BigDecimal getUpperBoundScore() {
		if(isComputed()) return getScore();

		BigDecimal den = getDenominator();
		if(den.compareTo(BigDecimal.ZERO) == 0) return toLog10(den);
		PartitionFunction pf = partitionFunctions[numStates-1];
		if(pf==null) return toLog10(BigDecimal.ZERO.setScale(64, RoundingMode.HALF_UP));
		BigDecimal num = pf.getValues().qstar.setScale(64, RoundingMode.HALF_UP);
		num = (num.add(pf.getValues().qprime)).add(pf.getValues().pstar);
		BigDecimal ans = num.divide(den, RoundingMode.HALF_UP);
		return toLog10(ans);
	}

	public String toString() {
		StringBuilder sb = new StringBuilder();
		//sb.append("Seq: "+settings.search[numStates-1].settings.getFormattedSequence()+", ");
		String seq = "";
		for(int state = 0; state < numStates-1; ++state) {
			seq += (settings.search[state].settings.getFormattedSequence()+" | ");
		}
		int ind = seq.lastIndexOf("|");
		if(ind>=0) seq = new StringBuilder(seq).replace(ind, ind+1,"").toString();
		seq = seq.trim();
		sb.append("Seq: "+seq+", ");
		
		sb.append(String.format("log10(score): %12e, ", getScore()));
		for(int state=0;state<numStates;++state) {
			BigDecimal qstar = partitionFunctions[state]==null ? BigDecimal.ZERO : 
				partitionFunctions[state].getValues().qstar;
			sb.append(String.format("pf: %2d, q*: %12e, ", state, qstar));
		}
		String ans = sb.toString().trim();
		return ans.substring(0,ans.length()-1);
	}

	protected boolean computeMinGMEC() {
		return settings.cfp.getParams().getBool("DOMINIMIZE") && 
				(computeMinGMECRatio() || computeGMEC());
	}

	protected boolean computeMinGMECRatio() {
		return computeGMECRatio() && 
				settings.cfp.getParams().getBool("DOMINIMIZE");
	}

	protected boolean computeGMEC() {
		return settings.isFinal && (settings.computeGMEC || settings.cfp.getParams().getBool("ComputeGMEC"));
	}

	protected boolean computeGMECRatio() {
		return settings.isFinal &&
				settings.cfp.getParams().getBool("ComputeGMECRatio");
	}

	protected List<EnergiedConf> findMinGMEC(int state) {
		List<EnergiedConf> energiedConfs = null;
		try {
			if(settings.isReportingProgress) {
				System.out.println();
				System.out.println("computing GMEC");
			}

			GMECFinder gmecFinder = new GMECFinder();
			gmecFinder.setLogConfsToConsole(settings.isReportingProgress);
			gmecFinder.isReportingProgress(settings.isReportingProgress);
			gmecFinder.init(settings.cfp, 
					settings.search[state],
					MSKStarFactory.makeConfSearchFactory(settings.search[state], settings.cfp),
					settings.ecalcs[state]);
			energiedConfs = gmecFinder.calcGMEC(settings.search[state]);

			if(settings.isReportingProgress) {
				System.out.println();
				System.out.println("done computing GMEC");
			}
		} catch (Exception e) {
			e.printStackTrace(System.out);
			throw new RuntimeException(e);
		}

		return energiedConfs != null && energiedConfs.size() > 0 ? energiedConfs : null;
	}

	protected boolean init(int state) {		
		if(settings.isReportingProgress) {
			System.out.println();
			System.out.println("state"+state+": "+settings.search[state].settings.getFormattedSequence()+" "+settings.pfTypes[state]);
		}

		//find the gmec if we are asked to do so
		//it's important find gmec here, before the prunepmat step, since gmec finding
		//sets the ival to a level required to find the gmec
		List<EnergiedConf> minGMECConfs = computeMinGMEC() ? findMinGMEC(state) : null;

		if(settings.isReportingProgress) {
			System.out.println();
			System.out.println("pruning matrix...");
		}

		//first prune the pruning matrix
		boolean doPruning = isFinal() || settings.cfp.getParams().getBool("PRUNEPARTIALSEQCONFS");
		settings.search[state].prunePmat(doPruning, 
				settings.cfp.getParams().getInt("ALGOPTION")>=3,
				settings.isReportingProgress);

		if(settings.isReportingProgress) {
			System.out.println("...finished pruning matrix");
		}

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
		pf.setComputeMaxNumConfs(computeMaxNumConfs(state));

		pf.setComputeGMECRatio(computeGMECRatio());
		if(minGMECConfs != null) pf.setMinGMECConfs(minGMECConfs);
		
		//create priority queue for top confs if requested
		if(settings.isFinal && settings.numTopConfsToSave > 0) {

			pf.topConfs = new PriorityQueue<ScoredConf>(
					settings.numTopConfsToSave, 
					new ConfComparator()
					);

			pf.maxNumTopConfs = settings.numTopConfsToSave;

			pf.setConfListener((ScoredConf conf) -> {
				pf.saveConf(conf);
			});

		}

		//init partition function
		if(settings.isReportingProgress) {
			System.out.println();
			System.out.print("initializing partition function...");
		}

		pf.init(settings.targetEpsilon);
		
		//special case: if init results in epsilon due
		//to gmec computation then print top confs
		if(pf.getStatus() == Status.Estimated && settings.numTopConfsToSave > 0) {
			pf.writeTopConfs(settings.state, 
					settings.search[state], 
					settings.cfp.getParams().getValue("TopConfsDir"));
		}

		if(settings.isReportingProgress) {
			System.out.println(" done");
		}

		return true;
	}

	/**
	 * compute until maxNumConfs conformations have been processed
	 * @param maxNumConfs
	 */
	public void compute(long maxNumConfs) {

		for(int state=0;state<numStates;++state) {

			if(!constrSatisfied) {//state-specific constraints 
				return;
			}

			if(!initialized[state]) {
				initialized[state] = init(state);
			}

			if(partitionFunctions[state].getStatus() != Status.Estimated) {
				compute(state, maxNumConfs);
			}
		}

		//check all constraints now. technically, we should only check constraints
		//that don't pertain to only one state, which we have already checked in compute(state)
		if(settings.isFinal && constrSatisfied) {
			constrSatisfied = checkConstraints();
		}

		if(isComputed()) {
			cleanup();
		}
	}

	/**
	 * compute only unbound states
	 * @param maxNumConfs
	 */
	@Override
	public void computeUnboundStates(long maxNumConfs) {
		for(int state=0;state<numStates-1;++state){

			if(!constrSatisfied)//state-specific constraints
				return;

			if(!initialized[state]) {
				initialized[state] = init(state);
			}

			if(partitionFunctions[state].getStatus() != Status.Estimated) {
				compute(state, maxNumConfs);
			}

			//don't check all constraints, because we are not computing 
			//the bound state partition function
			if(settings.isFinal && constrSatisfied) {
				constrSatisfied = checkConstraints(state);
			}

			partitionFunctions[state].cleanup();
		}
	}

	public void computeBoundState(long maxNumConfs) {
		if(!constrSatisfied)
			return;

		int state = numStates-1;
		if(!initialized[state]) {
			initialized[state] = init(state);
		}

		if(partitionFunctions[state].getStatus() != Status.Estimated) {
			compute(state, maxNumConfs);
		}

		if(partitionFunctions[state].getStatus()==Status.Estimated) {//assumption: unbound states are complete
			if(settings.isFinal && constrSatisfied) 
				constrSatisfied = checkConstraints();
		}

		if(isComputed()) {
			cleanup();
		}
	}

	private void cleanup() {
		for(PartitionFunctionMinimized pf : partitionFunctions) {
			if(pf==null) continue;
			pf.cleanup();
		}
	}

	/**
	 * Unprune the pruning matrix until the new p* is small enough to give an
	 * epsilon approx
	 */
	private PartitionFunction phase2(int state, long maxNumConfs) {

		PartitionFunctionMinimized pf = partitionFunctions[state];
		BigDecimal qstar = pf.getValues().qstar;
		BigDecimal pstar;

		PartitionFunctionMinimized p2pf = null;
		MSSearchProblem search = settings.search[state];
		//double oldPruningWindow = search.settings.pruningWindow;
		//double oldStericThreshold = search.settings.stericThreshold;
		double effectiveEpsilon = 1.0;

		do {
			// unprune
			search.settings.pruningWindow += 0.15;

			boolean doPruning = isFinal() || settings.cfp.getParams().getBool("PRUNEPARTIALSEQCONFS");
			search.prunePmat(doPruning, 
					settings.cfp.getParams().getInt("ALGOPTION")>=3,
					settings.isReportingProgress);

			// else, we got the additional number of confs, so we proceed
			//make conf search factory (i.e. A* tree)
			ConfSearchFactory confSearchFactory = MSKStarFactory.makeConfSearchFactory(search, settings.cfp);

			//create partition function
			p2pf = (PartitionFunctionMinimized) MSKStarFactory.makePartitionFunction( 
					settings.pfTypes[state],
					search.emat, 
					search.pruneMat,
					new PruningMatrixInverted(search, search.pruneMat),
					confSearchFactory,
					settings.ecalcs[state]
					);

			p2pf.setReportProgress(settings.isReportingProgress);
			p2pf.setMinGMECConfs(pf.getMinGMECConfs());
			p2pf.init(settings.targetEpsilon);

			pstar = p2pf.getValues().pstar;

			// NaN means qstar and pstar are both 0
			effectiveEpsilon = Values.getEffectiveEpsilon(qstar, BigDecimal.ZERO, pstar);

		} while(!Double.isNaN(effectiveEpsilon) && 
				effectiveEpsilon > settings.targetEpsilon &&
				search.settings.pruningWindow < search.settings.stericThreshold);

		/*
		// we can optionally prune to a steric threshold of 100
		final double maxStericThresh = 100.0;
		if(effectiveEpsilon > settings.targetEpsilon && search.settings.stericThreshold < maxStericThresh) {
			// we have not been successful within the limits of the steric threshold, 
			// so unprune to max steric threshold and proceed
			search.settings.pruningWindow = maxStericThresh;
			search.settings.stericThreshold = maxStericThresh;

			boolean doPruning = isFinal() || settings.cfp.getParams().getBool("PRUNEPARTIALSEQCONFS");
			search.prunePmat(doPruning, 
					settings.cfp.getParams().getInt("ALGOPTION")>=3,
					settings.isReportingProgress);

			// else, we got the additional number of confs, so we proceed
			//make conf search factory (i.e. A* tree)
			ConfSearchFactory confSearchFactory = MSKStarFactory.makeConfSearchFactory(search, settings.cfp);

			//create partition function
			p2pf = (PartitionFunctionMinimized) MSKStarFactory.makePartitionFunction( 
					settings.pfTypes[state],
					search.emat, 
					search.pruneMat,
					new PruningMatrixInverted(search, search.pruneMat),
					confSearchFactory,
					settings.ecalcs[state]
					);

			p2pf.setReportProgress(settings.isReportingProgress);
			p2pf.setMinGMEC(pf.getMinGMEC());
			p2pf.init(settings.targetEpsilon);
		}
		 */

		p2pf.compute(maxNumConfs);
		p2pf.setStatus(Status.Estimated);

		return p2pf;
	}

	//in the bound state, can override maxNumConfs with value from config
	private boolean computeMaxNumConfs(int state) {
		return settings.isFinal && state==numStates-1 && settings.cfp.getParams().getBool("ComputeMaxNumConfs");
	}

	protected void compute(int state, long maxNumConfs) {			
		PartitionFunctionMinimized pf = partitionFunctions[state];

		//in the bound state, can override maxNumConfs with value from config
		boolean computeMaxNumConfs = computeMaxNumConfs(state);
		if(computeMaxNumConfs) maxNumConfs = settings.cfp.getParams().getInt("MaxNumConfs");

		pf.compute(maxNumConfs);

		if(pf.getStatus() == Status.Estimated) {
			//yay!
		}

		else if(pf.getStatus() == Status.NotEnoughFiniteEnergies) {
			pf.setStatus(Status.Estimated);
		}

		else if(computeMaxNumConfs && pf.getNumConfsEvaluated()>=maxNumConfs) {
			pf.setStatus(Status.Estimated);
		}

		//we are not trying to compute the partition function to completion
		else if(!computeMaxNumConfs && pf.getStatus() == Status.Estimating) {
			return;
		}

		//no more q conformations, and we have not reached epsilon
		else if(pf.getStatus() == Status.NotEnoughConformations) {

			if(settings.isReportingProgress) System.out.print("Attempting phase 2 ... ");

			PartitionFunctionMinimized p2pf = (PartitionFunctionMinimized) phase2(state, maxNumConfs);
			partitionFunctions[state] = p2pf;

			if(settings.isReportingProgress) System.out.println("success");
		}

		if(isFinal()) {//final is a superset of fully defined
			if(constrSatisfied) constrSatisfied = checkConstraints(state);
			if(constrSatisfied && settings.numTopConfsToSave > 0) {
				pf.writeTopConfs(settings.state, settings.search[state], settings.cfp.getParams().getValue("TopConfsDir"));
			}
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
				stateVals[s] = toLog10(stateVals[s]);
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
				stateVals[s] = toLog10(stateVals[s]);
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

	@Override
	public ArrayList<MSSearchProblem> getUpperBoundSearch() {
		ArrayList<MSSearchProblem> ans = new ArrayList<>();
		for(int i=0; i<numStates; ++i) {
			PartitionFunctionMinimized pf = partitionFunctions[i];
			if(pf.getClass().getSimpleName().toLowerCase().contains("upperbound")) {
				ans.add(settings.search[i]);
			}
		}
		return ans;
	}

}
