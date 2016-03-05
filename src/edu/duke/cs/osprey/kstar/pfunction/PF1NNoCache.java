package edu.duke.cs.osprey.kstar.pfunction;

import java.util.ArrayList;

import edu.duke.cs.osprey.astar.ConfTree;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.kstar.KSConf;
import edu.duke.cs.osprey.pruning.PruningControl;

/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 *
 */
public class PF1NNoCache extends PFAbstract {

	protected ConfSearch search;

	// temp for benchmarking
	protected long startTime;

	public PF1NNoCache(ArrayList<String> sequence, ConfigFileParser cfp, 
			SearchProblem sp, PruningControl pc, DEEPerSettings dset, 
			ArrayList<String[]> moveableStrands, ArrayList<String[]> freeBBZones, 
			double EW_I0) {

		super( sequence, cfp, sp, pc, dset, moveableStrands, freeBBZones, EW_I0 );
	}


	public void start() {

		setRunState(RunState.STARTED);

		// replace new confrtree with a conf tree factory call 
		// to a function in the abstract base class
		search = new ConfTree(sp);
		int rawConf[];

		startTime = System.currentTimeMillis();

		if( (rawConf = search.nextConf()) != null ) {

			KSConf ksConf = new KSConf(rawConf, sp.lowerBound(rawConf), Double.MAX_VALUE);
			
			// get approx gmec LB to compute p*
			Et = ksConf.getMinEnergyLB();
			setPStar( Et );
			
			printedHeader = true;
			if( !minimizedConfsSet.contains(ksConf.getConf()) )
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
			
			KSConf conf = new KSConf(rawConf, sp.lowerBound(rawConf), Double.MAX_VALUE);

			accumulate(conf);

			if( eAppx == EApproxReached.NOT_POSSIBLE )
				System.out.println("Cannot reach epsilon");
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
		
		conf.setMinEnergy(sp.minimizedEnergy(conf.getConf()));

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
				+ getNumMinimizedConfs() + "\t" + getNumUnMinimizedConfs() + "\t"+ (currentTime-startTime)/1000);

		eAppx = effectiveEpsilon <= targetEpsilon || maxKSConfsReached() ? EApproxReached.TRUE: EApproxReached.FALSE;
	}


	protected void printHeader() {

		System.out.println("minE" + "\t\t\t" + "epsilon" + "\t\t" + "#min" +
				"\t" + "#un-min" + "\t" + "time(sec)");

		printedHeader = true;
	}

}