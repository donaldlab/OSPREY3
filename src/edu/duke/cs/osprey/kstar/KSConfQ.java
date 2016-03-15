package edu.duke.cs.osprey.kstar;

import java.io.Serializable;
import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Arrays;

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
	private BigDecimal capacityThresh = new BigDecimal(0.0001);
	private int size = 0;
	
	// lock for queue access
	public final String qLock = new String("LOCK");

	// upper bound partition function
	private BigDecimal qDagger = BigDecimal.ZERO;

	private final ArrayList<KSConf> q = new ArrayList<>();
	private int qCap = (int)Math.pow(2, 20);
	private int origQCap = 0;
	private boolean confsExhausted = false;

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
		return Double.MAX_VALUE;
	}


	public int size() {
		return size;
	}


	public boolean isExhausted() {
		return confsExhausted;
	}


	public boolean canSatisfy(int requested) {
		return size() >= requested;
	}


	public KSConf peekTail() {
		return q.size() > 0 ? q.get(q.size()-1) : null;
	}


	public KSConf peek() {
		return q.size() > 0 ? q.get(0) : null;
	}


	public double enQueue( int conf[] ) {
		
		double minELB = sp.lowerBound(conf);
		
		if( pf.getMinimizedConfsSet().contains(Arrays.toString(conf)) ) return minELB;

		KSConf ksc = new KSConf(conf, minELB);

		qDagger = qDagger.add( pf.getBoltzmannWeight(minELB) );

		q.add(ksc);
		
		size++;
		
		return minELB;
	}


	public KSConf deQueue() {
		// assuming locks are in place
		KSConf ksc = q.size() > 0 ? q.remove(0) : null;
		if(ksc == null) throw new RuntimeException("Error: attempting to dequeue from an empty list");

		// qDagger = qDagger.subtract( pf.getBoltzmannWeight(ksc.getELowerBound()) );
		
		size--;
		
		return ksc;
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


	public KSConf get( int index ) {
		return index >= 0 && index < q.size() ? q.get(index) : null;
	}


	public void cleanUp() throws InterruptedException {

		synchronized( this.qLock ) {
			this.qLock.notify();
		}
		
		this.join();
		
		this.sp = null;
		this.search = null;
		this.q.clear();
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

				if( q.size() >= qCap ) {
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
				if( q.size() == minCapacity ) qLock.notify();
			}
		}
	}

}
