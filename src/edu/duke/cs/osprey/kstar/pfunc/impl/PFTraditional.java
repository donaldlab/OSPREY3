package edu.duke.cs.osprey.kstar.pfunc.impl;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;

import cern.colt.Arrays;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.energy.MultiTermEnergyFunction;
import edu.duke.cs.osprey.kstar.KSAllowedSeqs;
import edu.duke.cs.osprey.kstar.KSConf;
import edu.duke.cs.osprey.kstar.KSConfigFileParser;
import edu.duke.cs.osprey.kstar.KSSearchProblem;
import edu.duke.cs.osprey.kstar.RCEnergyContribs;
import edu.duke.cs.osprey.kstar.pfunc.PFAbstract;

/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 *
 */
@SuppressWarnings("serial")
public class PFTraditional extends PFAbstract implements Serializable {
	
	protected ConfSearch confSearch = null;

	public PFTraditional() { 
		super();
	}

	public PFTraditional( int strand, ArrayList<String> sequence, 
			ArrayList<Integer> absolutePos, 
			String checkPointPath, String reducedSPName, 
			KSConfigFileParser cfp, KSSearchProblem panSP ) {

		super( strand, sequence, absolutePos, checkPointPath, reducedSPName, cfp, panSP );
	}


	public void start() {

		setRunState(RunState.STARTED);
		
		if(canUseHotByManualSelection())
			createHotsFromCFG();
		
		initTradPStar();

		confSearch = getConfTree(false);

		startTime = System.currentTimeMillis();
	}


	protected void iterate() {

		ScoredConf conf;

		if( (conf = confSearch.nextConf()) != null ) {

			if( processedConfsSet.contains(conf.getAssignments()) ) return;

			KSConf ksConf = new KSConf(conf.getAssignments(), getConfBound(confSearch, conf.getAssignments()));

			accumulate(ksConf);
		}

		else {
			// no more conformations, and we have not reached target epsilon
			if( eAppx == EApproxReached.FALSE )
				eAppx = EApproxReached.NOT_POSSIBLE;
		}
		
		exitIfTimeOut();
	}


	@Override
	protected void computeSlice() {
		iterate();
	}


	protected void compute() {
		while( eAppx == EApproxReached.FALSE ) {
			iterate();
		}
	}

	
	protected void createHotsFromCFG() {
		
		if(HOTs == null) HOTs = new ArrayList<>();
		
		ArrayList<String> flexRes = KSAllowedSeqs.getFlexResFromSeq(getSequence());
		ArrayList<ArrayList<String>> hots = cfp.getHighOrderTuplesByStrand(strand);
		
		for( ArrayList<String> hot : hots ) {
			
			ArrayList<Integer> hotIndexes = new ArrayList<>();
			
			for(String res : hot ) {
				
				int pos = flexRes.indexOf(res);
				
				if( HOTsContains(pos) ) break;
				
				hotIndexes.add(pos);
			}
			
			if(hotIndexes.size() > 2) {
				Collections.sort(hotIndexes);
				combineResidues(KSConf.list2Array(hotIndexes));
			}
		}
	}
	
	
	protected void combineResidues(int[] pos) {
		
		System.out.print("Combining residues: "); for(int i : pos) System.out.print(getSequence().get(i) + " ");
		System.out.print("... ");
		
		long start = System.currentTimeMillis();
		reducedSP.mergeResiduePositions(pos);
		memoizePosInHot(pos);
		long duration = (System.currentTimeMillis()-start)/1000;

		System.out.println("done in " + duration + "s");
	}
	

	protected void combineResidues(KSConf conf, double pbe, double tpbe, int[] tpce) {

		System.out.print("% bound error: " + pbe + ". ");
		System.out.print("% bound error from top "+ getHotNumRes() +" RCs: " + tpbe + ". positions: " + Arrays.toString(tpce));		
		System.out.print(". ");
		
		combineResidues(tpce);
		
		//BigDecimal oldPartialQPLB = new BigDecimal(partialQLB.toString());
		partialQLB = reComputePartialQLB(null);
		//if( oldPartialQPLB.compareTo(partialQLB) < 0 )
		//	throw new RuntimeException("ERROR: old partial q' - new partial q': " + oldPartialQPLB.subtract(partialQLB) + " must be >= 0");
		
		confSearch = getConfTree(false);
		conf.setEnergyBound(getConfBound(null, conf.getConfArray()));
	}


	protected void tryHotForConf(KSConf conf, MultiTermEnergyFunction mef) {

		double pbe = 0, tpbe = 0; int[] tpce = null;

		RCEnergyContribs rce = new RCEnergyContribs(this, mef, conf.getConfArray());
		pbe = rce.getPercentBoundError();

		if(pbe >= getHotBoundPct()) {
			tpbe = rce.getPercentErrorForTopPos(getHotNumRes());
			if(tpbe >= getHotTopRotsPct()) {
				tpce = rce.getTopPosCausingError(getHotNumRes());
				combineResidues(conf, pbe, tpbe, tpce);
			}
		}
	}


	protected void accumulate( KSConf conf ) {

		double energy = 0, boundError = 0;
		MultiTermEnergyFunction mef = null;

		if( isContinuous() && isFullyDefined() ) {
			// we do not have a lock when minimizing
			mef = reducedSP.decompMinimizedEnergy(conf.getConfArray());
			energy = mef.getPreCompE();
		}

		else energy = conf.getEnergyBound();

		conf.setEnergy(energy);
		boundError = conf.getEnergyBound() - conf.getEnergy();
		
		updateQStar( conf );

		Et = conf.getEnergyBound();
		updateQPrime();

		// negative values of effective esilon are disallowed
		if( (effectiveEpsilon = computeEffectiveEpsilon()) < 0 ) {

			eAppx = EApproxReached.NOT_POSSIBLE;

			return;
		}

		long currentTime = System.currentTimeMillis();

		if( !PFAbstract.suppressOutput ) {
			if( !printedHeader ) printHeader();

			System.out.println(numberFormat.format(boundError) + "\t" + numberFormat.format(energy) + "\t" 
					+ numberFormat.format(effectiveEpsilon) + "\t" + getNumProcessed() + "\t" 
					+ getNumUnEnumerated() + "\t"+ (currentTime-startTime)/1000);
		}

		eAppx = effectiveEpsilon <= targetEpsilon || maxKSConfsReached() ? EApproxReached.TRUE: EApproxReached.FALSE;

		// hot
		double peb = (conf.getEnergyBound()-conf.getEnergy())/conf.getEnergy();
		if(canUseHotByConfError(peb)) 
			tryHotForConf(conf, mef);

		// for partial sequences when doing KAstar
		if( !isFullyDefined() && eAppx == EApproxReached.TRUE ) adjustQStar();
	}


	protected void printHeader() {

		System.out.println("error" + "\t" + "energy" + "\t" + "epsilon" + "\t" + "#processed" +
				"\t" + "#un-enum" + "\t" + "time(sec)");

		printedHeader = true;
	}


	public String getImpl() {
		return "Traditional";
	}


	public void cleanup() {
		super.cleanup();
		confSearch = null;
	}

}