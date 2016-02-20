package edu.duke.cs.osprey.kstar;

import java.util.ArrayList;

import edu.duke.cs.osprey.astar.ConfTree;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.pruning.PruningControl;

/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 *
 */
public class PF1NNoCache extends PFAbstract {

	ConfSearch search;

	// temp for benchmarking
	long startTime;

	public PF1NNoCache(ArrayList<String> sequence, ConfigFileParser cfp, 
			SearchProblem sp, PruningControl pc, DEEPerSettings dset, 
			ArrayList<String[]> moveableStrands, ArrayList<String[]> freeBBZones, 
			double EW_I0) {

		super( sequence, cfp, sp, pc, dset, moveableStrands, freeBBZones, EW_I0);
	}


	protected void start() {

		setRunState(RunState.STARTED);

		// replace new confrtree with a conf tree factory call 
		// to a function in the abstract base class
		search = new ConfTree(sp);
		int conf[];

		startTime = System.currentTimeMillis();

		if( (conf = search.nextConf()) != null ) {

			KSConf ksConf = new KSConf(conf, sp.lowerBound(conf));
			ksConf.setMinEnergy(sp.minimizedEnergy(conf));
			Et = ksConf.getMinEnergyLowerBound();

			// get approx gmec LB to compute p*
			updateQStar( ksConf );

			Et = ksConf.getMinEnergyLowerBound();

			setPStar( Et );
		}

	}


	protected void iterate() throws Exception {

		int rawConf[];

		if( (rawConf = search.nextConf()) != null ) {

			KSConf conf = new KSConf(rawConf, sp.lowerBound(rawConf));

			if( (eAppx = accumulate(conf)) == EApproxReached.NOT_POSSIBLE )
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
	protected EApproxReached accumulate( KSConf conf ) {

		conf.setMinEnergy(sp.minimizedEnergy(conf.getConf()));
		double E = conf.getMinEnergy();
		Et = conf.getMinEnergyLowerBound();

		updateQStar( conf );

		// negative values of effective esilon are disallowed
		if( (effectiveEpsilon = computeEffectiveEpsilon()) < 0 ) return EApproxReached.NOT_POSSIBLE;

		long currentTime = System.currentTimeMillis();

		System.out.println(Et + "\t" + E + "\t" + effectiveEpsilon + "\t" 
				+ getNumMinimizedConfs() + "\t" + getNumUnMinimizedConfs() + "\t "+ ((currentTime-startTime)/1000));

		// return Et < threshold ? true : false;
		return effectiveEpsilon > targetEpsilon ? EApproxReached.FALSE : EApproxReached.TRUE;
	}

}