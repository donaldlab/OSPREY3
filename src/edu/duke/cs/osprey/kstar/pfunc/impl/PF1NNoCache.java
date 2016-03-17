package edu.duke.cs.osprey.kstar.pfunc.impl;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;

import edu.duke.cs.osprey.astar.ConfTree;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.kstar.KSConf;
import edu.duke.cs.osprey.kstar.pfunc.PFAbstract;

/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 *
 */
@SuppressWarnings("serial")
public class PF1NNoCache extends PFAbstract implements Serializable {

	protected ConfSearch search;

	// temp for benchmarking
	protected long startTime;

	public PF1NNoCache( ArrayList<String> sequence, String checkPointPath, 
			ConfigFileParser cfp, SearchProblem sp, double EW_I0 ) {

		super( sequence, checkPointPath, cfp, sp, EW_I0 );
	}


	public void start() {

		setRunState(RunState.STARTED);

		// replace new confrtree with a conf tree factory call 
		// to a function in the abstract base class
		search = new ConfTree(sp);
		int rawConf[];

		startTime = System.currentTimeMillis();

		if( (rawConf = search.nextConf()) != null ) {

			KSConf ksConf = new KSConf(rawConf, sp.lowerBound(rawConf));

			// get approx gmec LB to compute p*
			Et = ksConf.getMinEnergyLB();
			setPStar( Et );

			printedHeader = true;
			if( !minimizedConfsSet.contains(Arrays.toString(ksConf.getConfArray())) )
				accumulate( ksConf );
			printedHeader = false;
		}

		else
			eAppx = EApproxReached.NOT_POSSIBLE;
	}


	protected void iterate() throws Exception {

		int rawConf[];

		if( (rawConf = search.nextConf()) != null ) {

			if( minimizedConfsSet.contains(rawConf) ) return;

			KSConf conf = new KSConf(rawConf, sp.lowerBound(rawConf));

			accumulate(conf);
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


	/**
	 * Synchronous version evaluates conf energy immediately and adds to value.
	 */
	protected void accumulate( KSConf conf ) {

		conf.setMinEnergy(sp.minimizedEnergy(conf.getConfArray()));

		double E = conf.getMinEnergy();

		Et = conf.getMinEnergyLB();

		updateQStar( conf );

		// negative values of effective esilon are disallowed
		if( (effectiveEpsilon = computeEffectiveEpsilon()) < 0 ) {

			eAppx = EApproxReached.NOT_POSSIBLE;

			return;
		}

		long currentTime = System.currentTimeMillis();

		if( !printedHeader ) printHeader();

		System.out.println(E + "\t" + effectiveEpsilon + "\t" 
				+ getNumMinimized() + "\t" + getNumUnEnumerated() + "\t"+ (currentTime-startTime)/1000);

		eAppx = effectiveEpsilon <= targetEpsilon || maxKSConfsReached() ? EApproxReached.TRUE: EApproxReached.FALSE;
	}


	protected void printHeader() {

		System.out.println("minE" + "\t" + "epsilon" + "\t" + "#min" +
				"\t" + "#un-enum" + "\t" + "time(sec)");

		printedHeader = true;
	}

}