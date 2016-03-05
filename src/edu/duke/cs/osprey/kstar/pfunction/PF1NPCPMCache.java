package edu.duke.cs.osprey.kstar.pfunction;

import java.math.BigDecimal;
import java.math.BigInteger;
import java.util.ArrayList;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.kstar.KSConf;
import edu.duke.cs.osprey.kstar.KSConfQ;
import edu.duke.cs.osprey.pruning.PruningControl;
import edu.duke.cs.osprey.tools.ObjectIO;

/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 *
 */
public class PF1NPCPMCache extends PF1NPMCache {

	private ArrayList<Integer> indexes = new ArrayList<>();
	private ArrayList<SearchProblem> sps = new ArrayList<>();
	private ArrayList<KSConf> partialQConfs = new ArrayList<>();

	protected PF1NPCPMCache(ArrayList<String> sequence, ConfigFileParser cfp, 
			SearchProblem sp, PruningControl pc, DEEPerSettings dset, 
			ArrayList<String[]> moveableStrands, ArrayList<String[]> freeBBZones, 
			double EW_I0) {

		super( sequence, cfp, sp, pc, dset, moveableStrands, freeBBZones, EW_I0 );
	}


	protected void clearSearchProblem() {
		super.clearSearchProblem();
		sps.clear();
	}


	public void start() {

		try {

			setRunState(RunState.STARTED);

			// initialize parallel data structures
			for( int it = 0; it < PFAbstract.getNumThreads(); ++it ) indexes.add(indexes.size());
			indexes.trimToSize();

			for( int i = 0; i < indexes.size(); ++i )
				sps.add((SearchProblem)ObjectIO.deepCopy(sp));

			sps.trimToSize();

			partialQConfs.clear();
			for( int it = 0; it < indexes.size(); ++it ) partialQConfs.add(null);
			partialQConfs.trimToSize();

			confs = new KSConfQ( this, sp, indexes.size() );

			// set pstar
			setPStar( confs.getNextConfELB() );

			confs.start();

			if( waitUntilCapacity )
				confs.waitUntilCapacity();

			startTime = System.currentTimeMillis();

		} catch(Exception e) {
			System.out.println(e.getMessage());
			e.printStackTrace();
			System.exit(1);
		}

	}


	protected void iterate() throws Exception {

		synchronized( confs.qLock ) {

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
				partialQConfs.set(i, confs.deQueue());
			}

			minimizingConfs = minimizingConfs.add( BigInteger.valueOf(partialQConfs.size()) );

			if( confs.getState() == Thread.State.WAITING ) confs.qLock.notify();
		}

		// minimization hapens here
		accumulate(partialQConfs, false); 

		if( eAppx != EApproxReached.FALSE ) {
			// we leave this function
			confs.cleanUp();
		}	
	}


	protected void computeSlice() {

		try {
			/*
			synchronized( confs.qLock ) {
				if( confs.getState() == Thread.State.WAITING ) 
					confs.qLock.notify();
			}
			 */
			iterate();

		} catch(Exception e) {
			e.printStackTrace();
			System.exit(1);
		}

	}


	protected void compute() {

		try {

			// process conformations until e-approximation reached
			while( eAppx == EApproxReached.FALSE ) {

				iterate();
			}

		} catch(Exception e) {
			System.out.println(e.getMessage());
			e.printStackTrace();
			System.exit(1);
		}
	}


	protected BigDecimal computePartialQDagger( ArrayList<KSConf> partialQConfs ) {
		BigDecimal partialQDagger = BigDecimal.ZERO;

		for( KSConf conf : partialQConfs )
			partialQDagger = partialQDagger.add( getBoltzmannWeight(conf.getMinEnergyLB()) );

		return partialQDagger;
	}


	protected void accumulate( ArrayList<KSConf> partialQConfs, boolean isMinimized ) throws Exception {

		if( !isMinimized ) {
			// we do not have a lock when minimizing
			indexes.parallelStream().forEach( i -> {
				partialQConfs.get(i).setMinEnergy( sps.get(i).minimizedEnergy(partialQConfs.get(i).getConf()) );
			});
		}

		double E = 0;

		// we need a current snapshot of qDagger, so we lock here
		synchronized( confs.qLock ) {
			// update q*, qDagger, minimizingConfs, and q' atomically
			Et = confs.peekTail() != null ? confs.peekTail().getMinEnergyLB() 
					: partialQConfs.get(partialQConfs.size()-1).getMinEnergyLB();

			for( KSConf conf : partialQConfs ) {

				minimizingConfs = minimizingConfs.subtract( BigInteger.ONE );

				E = conf.getMinEnergy();
				updateQStar( conf );

				confs.setQDagger( confs.getQDagger().subtract( getBoltzmannWeight(conf.getMinEnergyLB()) ) );
				if(PFAbstract.useRigEUB) {
					confs.setQDot( confs.getQDot().subtract( getBoltzmannWeight(conf.getMinEnergyUB()) ) );
				}

				updateQPrime();

				// negative values of effective epsilon are disallowed
				if( (effectiveEpsilon = computeEffectiveEpsilon()) < 0 ) {
					eAppx = EApproxReached.NOT_POSSIBLE;
					return;
				}

				/*
				if(PFAbstract.useRigEnergy) {
					System.out.println("lb: " + conf.getMinEnergyLB() + "\t minE: " 
							+ conf.getMinEnergy() + "\t ub: " + conf.getMinEnergyUB() + "\t qstar/q.: " 
							+ qStar.divide(confs.getQDot(), 4)
							+ "\t q.size: " + confs.size());
				}
				 */

				if( effectiveEpsilon <= targetEpsilon || maxKSConfsReached() ) break;
			}

			long currentTime = System.currentTimeMillis();

			if( !printedHeader ) printHeader();

			System.out.println(E + "\t" + effectiveEpsilon + "\t" + 
					getNumMinimizedConfs() + "\t" + getNumUnMinimizedConfs() + "\t" + confs.size() + "\t" + ((currentTime-startTime)/1000));

			eAppx = effectiveEpsilon <= targetEpsilon || maxKSConfsReached() ? EApproxReached.TRUE: EApproxReached.FALSE;
		}
	}

}
