/*
** This file is part of OSPREY 3.0
** 
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
** 
** OSPREY relies on grants for its development, and since visibility
** in the scientific literature is essential for our success, we
** ask that users of OSPREY cite our papers. See the CITING_OSPREY
** document in this distribution for more information.
** 
** Contact Info:
**    Bruce Donald
**    Duke University
**    Department of Computer Science
**    Levine Science Research Center (LSRC)
**    Durham
**    NC 27708-0129
**    USA
**    e-mail: www.cs.duke.edu/brd/
** 
** <signature of Bruce Donald>, Mar 1, 2018
** Bruce Donald, Professor of Computer Science
*/

package edu.duke.cs.osprey.multistatekstar;

import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.PriorityQueue;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.gmec.ConfSearchFactory;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction.Status;
import edu.duke.cs.osprey.pruning.PruningMatrix;

/**
 * @author Adegoke Ojewole (ao68@duke.edu)
 * 
 */

public class KStarScoreMinimized implements KStarScore {

	public MSKStarSettings settings;
	public PartitionFunctionMinimized[] partitionFunctions;
	protected boolean[] initialized;
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

	@Override
	public MSKStarSettings getSettings() {
		return settings;
	}
	
	protected BigDecimal getDenom() {
		PartitionFunction pf;
		BigDecimal ans = BigDecimal.ONE.setScale(64, RoundingMode.HALF_UP);
		for(int state=0;state<numStates-1;++state) {
			pf = partitionFunctions[state];
			if(pf==null || pf.getValues().qstar.compareTo(BigDecimal.ZERO)==0)
				return BigDecimal.ZERO;
			ans = ans.multiply(pf.getValues().qstar);
		}
		return ans;
	}

	public BigDecimal getScore() {
		BigDecimal den = getDenom();
		if(den.compareTo(BigDecimal.ZERO) == 0) return BigDecimal.ZERO;
		PartitionFunction pf = partitionFunctions[numStates-1];
		return pf==null ? BigDecimal.ZERO : pf.getValues().qstar.divide(den, RoundingMode.HALF_UP);
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
		if(pf==null) return BigDecimal.ZERO;
		BigDecimal num = pf.getValues().qstar;
		if(pf.getStatus()!=Status.Estimated) num = num.add(pf.getValues().qprime).add(pf.getValues().pstar);
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

		partitionFunctions[state].setReportProgress(settings.isReportingProgress);

		//init partition function
		partitionFunctions[state].init(null, null, settings.targetEpsilon);

		//create priority queue for top confs if requested
		if(settings.search[state].isFullyAssigned() && settings.numTopConfsToSave > 0) {

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
		if(settings.isFinal && constrSatisfied) 
			constrSatisfied = checkConstraints();

		cleanup();
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

		p2pf.init(null, null, targetEpsilon);//enumerating over pstar, energies can be high
		p2pf.getValues().qstar = qstar;//keep old qstar
		p2pf.compute(targetScoreWeights);
		return p2pf;
	}

	protected void compute(int state, int maxNumConfs) {
		if(settings.isReportingProgress) 
			System.out.println("state"+state+": "+settings.search[state].settings.getFormattedSequence());
		PartitionFunctionMinimized pf = partitionFunctions[state];
		pf.compute(maxNumConfs);	

		//no more q conformations, and we have not reached epsilon
		double effectiveEpsilon = pf.getValues().getEffectiveEpsilon();
		if(!Double.isNaN(effectiveEpsilon) && effectiveEpsilon > settings.targetEpsilon) {
			PartitionFunctionMinimized p2pf = (PartitionFunctionMinimized) phase2(state);
			pf.getValues().qstar = p2pf.getValues().qstar;
			if(settings.search[state].isFullyAssigned() && settings.numTopConfsToSave > 0)
				pf.saveEConfs(p2pf.topConfs);
		}

		pf.setStatus(Status.Estimated);

		if(settings.isFinal) {//final is a superset of fully defined
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
	 * returns constraints that ONLY involve the spacified state
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
			if(constr.eval(stateVals).compareTo(BigDecimal.ZERO) > 0)
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
	public boolean isFullyProcessed() {
		if(!settings.isFinal) return false;
		int nulls = 0;
		for(PartitionFunctionMinimized pf : partitionFunctions) {
			if(pf==null) nulls++;
			else if(pf.getStatus() != Status.Estimated) return false;
		}
		//all non-null pfs are estimated; the reason why we skipped a pf must
		//be that a constraint is not satified
		if(nulls>0 && !constrSatisfied) return true;
		//otherwise, we erroneously skipped a partition function
		throw new RuntimeException("ERROR: illegally skipped a partition function computation");
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
