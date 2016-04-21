package edu.duke.cs.osprey.kstar.pfunc.impl;

import java.io.Serializable;
import java.math.BigInteger;
import java.util.ArrayList;

import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.kstar.KSConf;
import edu.duke.cs.osprey.kstar.KSConfQ;
import edu.duke.cs.osprey.kstar.QPrimeConfTree;
import edu.duke.cs.osprey.kstar.pfunc.PFAbstract;
import edu.duke.cs.osprey.tools.ObjectIO;

/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 *
 */
@SuppressWarnings("serial")
public class PFNew02 extends PFNew01 implements Serializable {

	private ArrayList<Integer> indexes = new ArrayList<>();
	private ArrayList<SearchProblem> sps = new ArrayList<>();
	private ArrayList<KSConf> partialQConfs = new ArrayList<>();

	public PFNew02() {
		super();
	}
	
	public PFNew02( int strand, ArrayList<String> sequence, ArrayList<Integer> flexResIndexes, 
			String checkPointPath, String searchProblemName, 
			ConfigFileParser cfp, SearchProblem panSeqSP ) {

		super( strand, sequence, flexResIndexes, checkPointPath, searchProblemName, cfp, panSeqSP );
	}


	public void cleanup() {
		super.cleanup();
		sps.clear();
	}


	protected ArrayList<SearchProblem> parallelCreateSPs( SearchProblem sp, int replicates ) {
		ArrayList<SearchProblem> ans = new ArrayList<>();
		ArrayList<Integer> indexes = new ArrayList<>();
		for(int i = 0; i < replicates; ++i) {
			indexes.add(i);
			ans.add(null);
		}

		indexes.parallelStream().forEach(i -> {
			ans.set(i, (SearchProblem)ObjectIO.deepCopy(sp));
		});

		indexes.trimToSize();
		ans.trimToSize();

		return ans;
	}


	public void start() {

		try {

			setRunState(RunState.STARTED);

			// set pstar
			initPStar();

			// initialize parallel data structures
			indexes.clear();
			for( int it = 0; it < PFAbstract.getNumThreads(); ++it ) indexes.add(it);
			indexes.trimToSize();

			sps.clear();
			sps = parallelCreateSPs(sp, indexes.size());

			partialQConfs.clear();
			for( int it = 0; it < indexes.size(); ++it ) partialQConfs.add(null);
			partialQConfs.trimToSize();

			confsQ = new KSConfQ( this, indexes.size() );
			qPrimeCalculator = new QPrimeConfTree( this, unPrunedConfs );
			
			qPrimeCalculator.start();
			confsQ.start();

		} catch(Exception e) {
			System.out.println(e.getMessage());
			e.printStackTrace();
			System.exit(1);
		}

		startTime = System.currentTimeMillis();
	}


	protected void iterate() throws Exception {

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

			minimizingConfs = minimizingConfs.add( BigInteger.valueOf(partialQConfs.size()) );

			if( confsQ.getState() == Thread.State.WAITING ) confsQ.lock.notify();
		}

		// minimization hapens here
		accumulate(partialQConfs, false); 

		if( eAppx != EApproxReached.FALSE ) {
			// we leave this function
			confsQ.cleanUp(true);
			qPrimeCalculator.cleanUp(true);
		}	
	}


	protected void accumulate( ArrayList<KSConf> partialQConfs, boolean isMinimized ) throws Exception {
		
		if( !isMinimized ) {
			
			// we do not have a lock when minimizing
			indexes.parallelStream().forEach( i -> {
				
				KSConf conf = partialQConfs.get(i);
				
				double energy = isFullyDefined() ? 
						sps.get(i).minimizedEnergy( conf.getConfArray() ) : conf.getEnergyBound();
				
				conf.setEnergy( energy );
			});
		}

		double energy = 0, boundError = 0;

		// we need a current snapshot of qDagger, so we lock here
		synchronized( confsQ.lock ) {
			
			while( confsQ.getPartialQLB().compareTo(qPrimeCalculator.getTotalQLB()) > 0 ) {
			//while( BigInteger.valueOf(confsQ.size()).compareTo(qPrimeCalculator.getNumEnumerated()) > 0 ) {
				Thread.sleep(5);
			}
			
			for( KSConf conf : partialQConfs ) {

				minimizingConfs = minimizingConfs.subtract( BigInteger.ONE );

				energy = conf.getEnergy();
				updateQStar( conf );

				boundError = Math.abs(conf.getEnergyBound()-conf.getEnergy())/Math.abs(conf.getEnergy())*100;
				
				confsQ.accumulatePartialQLB(getBoltzmannWeight(conf.getEnergyBound()));
				updateQPrime();

				// negative values of effective epsilon are disallowed
				if( (effectiveEpsilon = computeEffectiveEpsilon()) < 0 ) {
					eAppx = EApproxReached.NOT_POSSIBLE;
					return;
				}

				if( effectiveEpsilon <= targetEpsilon || maxKSConfsReached() ) break;
			}

			long currentTime = System.currentTimeMillis();

			if( !PFAbstract.suppressOutput ) {
				if( !printedHeader ) printHeader();
				
				System.out.println(boundError + "\t" + energy + "\t" + effectiveEpsilon + "\t" + 
						getNumMinimized4Output() + "\t" + getNumUnEnumerated() + "\t" + confsQ.size() + "\t" + ((currentTime-startTime)/1000));
			}

			eAppx = effectiveEpsilon <= targetEpsilon || maxKSConfsReached() ? EApproxReached.TRUE: EApproxReached.FALSE;

			// for partial sequences when doing KAstar
			if( !isFullyDefined() && eAppx == EApproxReached.TRUE ) adjustQStar();
		}
	}


	public String getImpl() {
		return "new02";
	}

}
