package edu.duke.cs.osprey.kstar;

import java.io.Serializable;
import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.LinkedHashSet;

import edu.duke.cs.osprey.astar.ConfTree;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.kstar.pfunc.PFAbstract;
import edu.duke.cs.osprey.kstar.pfunc.PFAbstract.EApproxReached;

/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 *
 */
@SuppressWarnings("serial")
public class KSConfQ extends Thread implements Serializable {

	private PFAbstract pf;
	private SearchProblem sp;
	private ConfSearch search;
	private int minCapacity;
	private BigDecimal capacityThresh = new BigDecimal(0.000001);
	
	// lock for queue access
	public final String qLock = new String("LOCK");

	// upper bound partition function
	private BigDecimal qDagger = BigDecimal.ZERO;

	LinkedHashSet<ArrayList<Integer>> q = null;
	private int qCap = (int)Math.pow(2, 20);
	private int origQCap = 0;
	private boolean confsExhausted = false;
	private ArrayList<Integer> tail = null;

	/**
	 * 
	 * @param pf
	 * @param search
	 * @param notificationThreshold = notify queue owner when queue contains
	 * this number of conformations
	 */
	public KSConfQ( PFAbstract pf, SearchProblem sp, int minCapacity ) {

		this.pf = pf;
		this.sp = sp;
		search = new ConfTree(sp);

		this.minCapacity = minCapacity;
		qCap = Math.max( minCapacity, PFAbstract.qCapacity );
		origQCap = qCap;
		
		q = new LinkedHashSet<>(qCap);
	}


	public void restartConfTree() {
		search = new ConfTree(sp);
	}
	
	
	public void waitUntilCapacity() throws InterruptedException {

		while( !isExhausted() && size() < getQCapacity() )
			Thread.sleep(250);
	}
	

	public double getNextConfELB() {

		int c[] = null;

		if( (c = search.nextConf()) != null ) {
			return enQueue(c);
		}

		// should never get here
		throw new RuntimeException("ERROR: all the conformations of this sequence were pruned");
		
		// return Double.MAX_VALUE;
	}


	public int size() {
		return q.size();
	}


	public boolean isExhausted() {
		return confsExhausted;
	}


	public boolean canSatisfy(int requested) {
		return size() >= requested;
	}


	public KSConf peekTail() {
		return size() > 0 ? new KSConf(tail, sp.lowerBound(KSConf.list2Array(tail))) : null;
	}


	public KSConf peekHead() {
		if(size() == 0) return null;
		
		ArrayList<Integer> value = q.iterator().next();
		return new KSConf(value, sp.lowerBound(KSConf.list2Array(value)));
	}
	

	protected double enQueue( int[] conf ) {
		
		double minELB = sp.lowerBound(conf);
		ArrayList<Integer> list = KSConf.array2List(conf);
		
		if(KSAbstract.doCheckpoint && size() > 0 && minELB < peekTail().getMinEnergyLB() ) return minELB;
		
		if( pf.getMinimizedConfsSet().contains(list) || q.contains(list) ) return minELB;

		qDagger = qDagger.add( pf.getBoltzmannWeight(minELB) );
		
		q.add(list);
		
		tail = list;
		
		return minELB;
	}


	public KSConf deQueue() {
		// assuming locks are in place
		KSConf conf = size() > 0 ? peekHead() : null;

		if(conf == null) 
			throw new RuntimeException("ERROR: attempting to dequeue from an empty list");
		
		q.remove(conf.getConf());
		
		if(size() == 0) tail = null;
		
		return conf;
	}


	public BigDecimal getCapacityThresh() {
		return capacityThresh;
	}


	public BigDecimal getQDagger() {
		return qDagger;
	}

	
	public void setQDagger( BigDecimal qDagger ) {
		this.qDagger = qDagger;
	}


	public int getQCapacity() {
		return qCap;
	}


	public void restoreQCapacity() {
		setQCapacity(origQCap);
	}


	public void setQCapacity( int newCap ) {
		qCap = Math.max(newCap, minCapacity);
	}


	public void cleanUp( boolean nullify ) throws InterruptedException {

		synchronized( qLock ) {
			qLock.notify();
		}
		
		this.join();
		
		if(nullify) {
			search = null;
			sp = null;
			q.clear();
			tail = null;
		}
	}


	public void run() {

		int conf[];

		while( true ) {

			conf = search.nextConf();

			synchronized( qLock ) {

				if( conf == null ) {
					confsExhausted = true;
					qLock.notify();
					return;
				}

				if( size() >= qCap ) {
					try {

						qLock.notify();

						if( pf.getEpsilonStatus() != EApproxReached.FALSE ) 
							return;

						else
							qLock.wait();

					} catch (InterruptedException e) {
						System.out.println(e.getMessage());
						e.printStackTrace();
						System.exit(1);
					}
				}

				// exit thread if we have an e-approximation
				if( pf.getEpsilonStatus() != EApproxReached.FALSE ) 
					return;

				enQueue(conf);

				// notify queue consumer ONLY if queue was empty before
				// i added latest conformation. this condition means that
				// the partition function is waiting for this signal to process
				// conformations.
				// it's wasteful to call notify for every insertion
				if( size() == minCapacity ) qLock.notify();
			}
		}
	}

}
