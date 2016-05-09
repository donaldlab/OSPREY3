package edu.duke.cs.osprey.kstar.pfunc.impl;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;

import cern.colt.Arrays;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.energy.MultiTermEnergyFunction;
import edu.duke.cs.osprey.kstar.AllowedSeqs;
import edu.duke.cs.osprey.kstar.KSConf;
import edu.duke.cs.osprey.kstar.RCEnergyContribs;
import edu.duke.cs.osprey.kstar.pfunc.PFAbstract;

/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 *
 */
@SuppressWarnings("serial")
public class PFTrad extends PFAbstract implements Serializable {

	protected ConfSearch confSearch = null;

	// temp for benchmarking
	protected long startTime;

	public PFTrad() { 
		super();
	}

	public PFTrad( int strand, ArrayList<String> sequence, ArrayList<Integer> flexResIndexes, 
			String checkPointPath, String searchProblemName, 
			ConfigFileParser cfp, SearchProblem panSeqSP ) {

		super( strand, sequence, flexResIndexes, checkPointPath, searchProblemName, cfp, panSeqSP );
	}


	public void start() {

		setRunState(RunState.STARTED);

		if(canUseHotByManualSelection()) 
			createHotsFromCFG();
		
		initPStar();

		// first conf was merely to set p*
		confSearch = getConfTree(false);

		startTime = System.currentTimeMillis();
	}


	protected void iterate() throws Exception {

		int conf[];

		if( (conf = confSearch.nextConf()) != null ) {

			if( minimizedConfsSet.contains(conf) ) return;

			KSConf ksConf = new KSConf(conf, getConfBound(confSearch, conf, false));

			accumulate(ksConf);
		}

		else {
			// no more conformations, and we have not reached target epsilon
			if( eAppx == EApproxReached.FALSE )
				eAppx = EApproxReached.NOT_POSSIBLE;
		}
	}


	@Override
	protected void computeSlice() {

		try {

			iterate();

		} catch(Exception e) {
			System.out.println(e.getMessage());
			e.printStackTrace();
			System.exit(1);
		}
	}


	protected void compute() {

		try {

			while( eAppx == EApproxReached.FALSE ) {

				iterate();

			}

		} catch(Exception e) {
			System.out.println(e.getMessage());
			e.printStackTrace();
			System.exit(1);
		}
	}

	
	protected void createHotsFromCFG() {
		
		ArrayList<String> flexRes = AllowedSeqs.getFlexResFromSeq(getSequence());
		ArrayList<ArrayList<String>> hots = cfp.getHighOrderTuplesByStrand(strand);
		
		for( ArrayList<String> hot : hots ) {
			
			ArrayList<Integer> hotIndexes = new ArrayList<>();
			
			for(String res : hot ) {
				
				int pos = flexRes.indexOf(res);
				
				if( posInHot.contains(pos) ) break;
				
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
		sp.mergeResiduePositions(pos);
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
		conf.setEnergyBound(getConfBound(null, conf.getConfArray(), false));
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

		if( isFullyDefined() ) {
			// we do not have a lock when minimizing
			mef = sp.decompMinimizedEnergy(conf.getConfArray());
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

			System.out.println(boundError + "\t" + energy + "\t" + effectiveEpsilon + "\t" + getNumMinimized4Output() + 
					"\t" + getNumUnEnumerated() + "\t"+ (currentTime-startTime)/1000);
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

		System.out.println("boundError" + "\t" + "minE" + "\t" + "epsilon" + "\t" + "#min" +
				"\t" + "#un-enum" + "\t" + "time(sec)");

		printedHeader = true;
	}


	public String getImpl() {
		return "trad";
	}


	public void cleanup() {
		super.cleanup();
		confSearch = null;
	}

}