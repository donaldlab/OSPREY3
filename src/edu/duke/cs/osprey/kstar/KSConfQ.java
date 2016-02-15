package edu.duke.cs.osprey.kstar;

import java.math.BigDecimal;
import java.util.ArrayList;

import edu.duke.cs.osprey.astar.ConfTree;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.kstar.PFAbstract.EApproxReached;

public class KSConfQ extends Thread {

	private PFAbstract pf;
	private SearchProblem sp;
	private ConfSearch search;
	private int minCapacity;

	// lock for queue access
	public final Object qLock = new Object();

	// upper bound partition function
	private BigDecimal qDagger = BigDecimal.ZERO;

	private final ArrayList<KSConf> q = new ArrayList<>();
	private int qCap = (int)Math.pow(2, 20);
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
	}


	public void waitUntilCapacity() throws InterruptedException {
		
		while( !exhausted() && size() < getQCapacity() )
			Thread.sleep(250);
	}

	
	public boolean getNextConf() {

		int c[] = null;

		if( (c = search.nextConf()) != null ) {
			enQueue(c);
			return true;
		}

		return false;
	}


	public int size() {
		return q.size();
	}


	public boolean exhausted() {
		return confsExhausted;
	}

	
	public KSConf peekTail() {
		return q.size() > 0 ? q.get(q.size()-1) : null;
	}

	
	public KSConf peek() {
		return q.size() > 0 ? q.get(0) : null;
	}


	public void enQueue( int conf[] ) {
		double minELowerBound = Double.POSITIVE_INFINITY;

		minELowerBound = sp.lowerBound(conf);
		qDagger = qDagger.add( pf.getBoltzmannWeight(minELowerBound) );

		KSConf ksc = new KSConf(conf, minELowerBound);
		q.add(ksc);
	}


	public KSConf deQueue() {
		// assuming locks are in place
		KSConf ksc = q.size() > 0 ? q.remove(0) : null;
		if(ksc == null) throw new RuntimeException("Error: attempting to dequeue from an empty list");

		// qDagger = qDagger.subtract( pf.getBoltzmannWeight(ksc.getELowerBound()) );
		return ksc;
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


	public void setQCapacity( int newCap ) {
		qCap = newCap;
	}


	public KSConf get( int index ) {
		return index >= 0 && index < q.size() ? q.get(index) : null;
	}


	public void cleanUp() throws InterruptedException {

		synchronized( this.qLock ) {
			this.qLock.notify();
		}

		this.join();
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

						qLock.wait();

					} catch (InterruptedException e) {
						System.out.println(e.getMessage());
						e.printStackTrace();
						System.exit(1);
					}
				}

				// exit thread if we have an e-approximation
				if( pf.eAppx != EApproxReached.FALSE ) 
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
