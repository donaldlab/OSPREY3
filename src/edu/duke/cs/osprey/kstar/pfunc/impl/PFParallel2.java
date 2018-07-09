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

package edu.duke.cs.osprey.kstar.pfunc.impl;

import java.io.Serializable;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.HashMap;

import edu.duke.cs.osprey.energy.MultiTermEnergyFunction;
import edu.duke.cs.osprey.kstar.KSConf;
import edu.duke.cs.osprey.kstar.KSConfQ;
import edu.duke.cs.osprey.kstar.KSConfigFileParser;
import edu.duke.cs.osprey.kstar.KSSearchProblem;
import edu.duke.cs.osprey.kstar.KSPStarCalculator;
import edu.duke.cs.osprey.kstar.KSQPrimeCalculator;
import edu.duke.cs.osprey.kstar.pfunc.PFAbstract;
import edu.duke.cs.osprey.tools.ObjectIO;

/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 *
 */
@SuppressWarnings("serial")
public class PFParallel2 extends PFParallel1 implements Serializable {

	protected static HashMap<Integer, ArrayList<KSSearchProblem>> strand2SPs = new HashMap<>();
	private ArrayList<Integer> indexes = new ArrayList<>(PFAbstract.getNumThreads());
	private ArrayList<KSSearchProblem> sps = null;
	private ArrayList<KSConf> partialQConfs = new ArrayList<>(PFAbstract.getNumThreads());

	public PFParallel2() {
		super();
	}

	public PFParallel2( int strand, ArrayList<String> sequence, 
			ArrayList<Integer> absolutePos, 
			String checkPointPath, String reducedSPName, 
			KSConfigFileParser cfp, KSSearchProblem panSP ) {

		super( strand, sequence, absolutePos, checkPointPath, reducedSPName, cfp, panSP );
	}


	public void cleanup() {
		super.cleanup();
		resetSPs();
	}


	protected ArrayList<KSSearchProblem> parallelCreateSPs( KSSearchProblem sp, int replicates ) {
		ArrayList<KSSearchProblem> ans = new ArrayList<>(replicates);
		ArrayList<Integer> indexes = new ArrayList<>(replicates);
		for(int i = 0; i < replicates; ++i) {
			indexes.add(i);
			ans.add(null);
		}

		indexes.parallelStream().forEach(i -> {
			ans.set(i, (KSSearchProblem) ObjectIO.deepCopy(sp));
			ans.get(i).emat = sp.emat;
			ans.get(i).pruneMat = sp.pruneMat;
			ans.get(i).inverseMat = sp.inverseMat;
		});

		return ans;
	}
	
	
	protected void appendSPs( ArrayList<KSSearchProblem> sps, KSSearchProblem sp, int target ) {
		
		while(sps.size() < target) {
			KSSearchProblem newSP = (KSSearchProblem) ObjectIO.deepCopy(sp);
			newSP.emat = sp.emat;
			newSP.pruneMat = sp.pruneMat;
			newSP.inverseMat = sp.inverseMat;
			sps.add(newSP);
		}
		
		sps.trimToSize();
	}


	protected void getSPs() {
		
		if(!strand2SPs.containsKey(strand))
			strand2SPs.put(strand, parallelCreateSPs(panSP, PFAbstract.getNumThreads()));
		
		sps = strand2SPs.get(strand);
		
		// make sure sps are of the right size
		if(sps.size() < PFAbstract.getNumThreads()) {
			appendSPs(sps, panSP, PFAbstract.getNumThreads());
			strand2SPs.put(strand, sps);
		}
	}


	protected void resetSPs() {
		sps = null;
	}


	public void start() {

		try {

			setRunState(RunState.STARTED);

			if(canUseHotByManualSelection()) 
				createHotsFromCFG();

			// set pstar
			pStarCalculator = null;
			if(prunedConfs.compareTo(BigInteger.ZERO) == 0) pStar = BigDecimal.ZERO;
			else {
				System.out.println("using p* calculator");
				pStarCalculator = new KSPStarCalculator( this );
				pStarCalculator.setPriority(Thread.MAX_PRIORITY);
			}

			getSPs();

			// initialize parallel data structures
			indexes.clear();
			for( int it = 0; it < PFAbstract.getNumThreads(); ++it ) indexes.add(it);
			indexes.trimToSize();

			partialQConfs.clear();
			for( int it = 0; it < indexes.size(); ++it ) partialQConfs.add(null);
			partialQConfs.trimToSize();

			confsQ = new KSConfQ( this, indexes.size(), partialQLB );
			qPrimeCalculator = new KSQPrimeCalculator( this );
			qPrimeCalculator.setPriority(Thread.MAX_PRIORITY);

			if(pStarCalculator != null) pStarCalculator.start();
			qPrimeCalculator.start();
			confsQ.start();
			
			if(!isContinuous() && isFullyDefined()) Thread.sleep(initSleepTime);

		} catch(Exception ex) {
			throw new Error("can't compute partition function", ex);
		}

		startTime = System.currentTimeMillis();
	}


	protected void iterate() {
		try {

			getSPs();

			synchronized( confsQ.lock ) {

				int request = partialQConfs.size();
				int granted = 0;

				if( (granted = canSatisfy(request)) == 0 )
					return;

				// reduce the size of partialQconfs and indexes to match request
				while( partialQConfs.size() > granted ) {

					partialQConfs.remove(partialQConfs.size()-1);

					indexes.remove(indexes.size()-1);
				}

				for( int i = 0; i < Math.min(granted, partialQConfs.size()); ++i ) {
					partialQConfs.set(i, confsQ.deQueue());
				}

				processingConfs = processingConfs.add( BigInteger.valueOf(partialQConfs.size()) );

				if( confsQ.getState() == Thread.State.WAITING ) confsQ.lock.notify();
			}

			// minimization hapens here
			accumulate(partialQConfs, false); 

			if( eAppx != EApproxReached.FALSE ) {
				// we leave this function
				confsQ.cleanUp(true);
				qPrimeCalculator.cleanUp(true);
				if(pStarCalculator != null) pStarCalculator.cleanUp(true);
			}

			resetSPs();
			
			exitIfTimeOut();
			
		} catch (InterruptedException ex) {
			// something interrupted us because it wants us to stop,
			// so throw an exception that no one's supposed to catch
			// and hopefully bubble to the top of the current thread
			throw new Error(ex);
		}
	}


	protected void tryHotForConfs( ArrayList<MultiTermEnergyFunction> mefs ) {

		for( int index : indexes ) {
			KSConf conf = partialQConfs.get(index);
			// update energy bound if updated emat
			conf.setEnergyBound(getConfBound(null, conf.getConfArray()));
			double peb = (conf.getEnergyBound()-conf.getEnergy())/conf.getEnergy();
			if(canUseHotByConfError(peb)) tryHotForConf(partialQConfs.get(index), mefs.get(index));
		}
	}
	

	protected void accumulate( ArrayList<KSConf> partialQConfs, boolean energiesEvaluated ) {

		ArrayList<MultiTermEnergyFunction> mefs = new ArrayList<>(partialQConfs.size());
		for(int i = 0; i < partialQConfs.size(); ++i) mefs.add(null);

		if( !energiesEvaluated ) {

			// we do not have a lock when minimizing
			indexes.parallelStream().forEach( i -> {

				double energy = 0;
				KSConf conf = partialQConfs.get(i);

				if( isContinuous() && isFullyDefined() ) {
					MultiTermEnergyFunction mef = sps.get(i).decompMinimizedEnergy(conf.getConfArray());
					mefs.set(i, mef);
					energy = mef.getPreCompE();
				}

				else energy = conf.getEnergyBound();

				conf.setEnergy(energy);
			});
		}

		double energy = 0, boundError = 0;
		for( KSConf conf : partialQConfs ) {

			processingConfs = processingConfs.subtract( BigInteger.ONE );

			energy = conf.getEnergy();
			boundError = conf.getEnergyBound() - conf.getEnergy();

			updateQStar( conf );

			updateQPrime();

			updatePStar();

			// negative values of effective epsilon are disallowed
			if( (effectiveEpsilon = computeEffectiveEpsilon()) < 0 ) {
				eAppx = EApproxReached.NOT_POSSIBLE;
				return;
			}

			if( effectiveEpsilon <= targetEpsilon || maxKSConfsReached() ) break;
		}

		long currentTime = System.currentTimeMillis();

		synchronized( confsQ.lock ) {

			if( !PFAbstract.suppressOutput ) {
				if( !printedHeader ) printHeader();
				System.out.println(numberFormat.format(boundError) + "\t" + numberFormat.format(energy) + "\t" 
						+ numberFormat.format(effectiveEpsilon) + "\t" + getNumProcessed() + "\t" 
						+ getNumUnEnumerated() + "\t" + confsQ.size() + "\t" + ((currentTime-startTime)/1000));
			}
		}

		eAppx = effectiveEpsilon <= targetEpsilon || maxKSConfsReached() ? EApproxReached.TRUE: EApproxReached.FALSE;

		// hot
		if(canUseHotByConfError()) 
			tryHotForConfs(mefs);

		// for partial sequences when doing KAstar
		if( !isFullyDefined() && eAppx == EApproxReached.TRUE ) adjustQStar();
	}


	public String getImpl() {
		return "Parallel2";
	}

}
