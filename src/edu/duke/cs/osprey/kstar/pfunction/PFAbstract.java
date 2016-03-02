package edu.duke.cs.osprey.kstar.pfunction;

import java.io.File;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.PriorityQueue;

import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.kstar.KSAbstract;
import edu.duke.cs.osprey.kstar.KSConf;
import edu.duke.cs.osprey.pruning.PruningControl;
import edu.duke.cs.osprey.tools.ExpFunction;
import edu.duke.cs.osprey.tools.ObjectIO;

/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 *
 */
public abstract class PFAbstract {

	protected ArrayList<String> sequence;
	protected static String pFuncImplementation = "1npmcache";
	public static String eMinMethod = "ccd";
	protected static ArrayList<String> serverList = new ArrayList<>();
	protected static int threadConfsBuffer = 8;
	public static int qCapacity = (int)Math.pow(2, 17);
	public static boolean useRigEnergy = false;
	public static boolean waitUntilCapacity = false;
	protected static int numThreads = 1;
	protected static int numFibers = 1;
	protected static int numRemoteClients = 1;
	
	protected boolean printedHeader = false;

	public static boolean saveTopConfsAsPDB = false;
	protected static int numTopConfsToSave = 10;

	protected static BigDecimal stabilityThresh = BigDecimal.ONE;
	protected static double interval = getMaxInterval();

	protected static final double RT = 1.9891/1000.0 * 298.15;
	public static double targetEpsilon = 0.03;
	public static double rho = targetEpsilon / (1.0 - targetEpsilon);
	protected double effectiveEpsilon = 1.0;

	public static enum EApproxReached { TRUE, FALSE, NOT_POSSIBLE, NOT_STABLE, ABORTED }
	protected EApproxReached eAppx = EApproxReached.FALSE;

	public static enum RunState { NOTSTARTED, STARTED, SUSPENDED, TERMINATED }
	protected RunState runState = RunState.NOTSTARTED;

	protected ConfigFileParser cfp = null;
	protected SearchProblem sp = null;
	protected PruningControl pc = null;
	protected DEEPerSettings dset = null;
	protected ArrayList<String[]> moveableStrands = null;
	protected ArrayList<String[]> freeBBZones = null;

	protected BigDecimal qStar = BigDecimal.ZERO;
	protected BigDecimal qPrime = BigDecimal.ZERO;
	protected BigDecimal pStar = BigDecimal.ZERO;

	protected ExpFunction e = new ExpFunction();
	protected double EW_I0 = 5;
	protected double Et = 0;

	protected BigInteger prunedConfs = BigInteger.ZERO;
	protected BigInteger initialUnPrunedConfs = BigInteger.ZERO;
	protected BigInteger phaseOneMinimizedConfs = BigInteger.ZERO;
	protected BigInteger minimizedConfs = BigInteger.ZERO;
	protected BigInteger minimizedConfsDuringInterval = BigInteger.ZERO;
	protected BigInteger minimizingConfs = BigInteger.ZERO; // # confs being minimized at this instant

	protected PriorityQueue<KSConf> topConfsPQ = null;

	protected PFAbstract( ArrayList<String> sequence, ConfigFileParser cfp, 
			SearchProblem sp, PruningControl pc,
			DEEPerSettings dset, ArrayList<String[]> moveableStrands, ArrayList<String[]> freeBBZones,
			double EW_I0 ) {

		this.sequence = sequence;
		this.sp = sp;
		this.pc = pc;
		this.dset = dset;
		this.moveableStrands = moveableStrands;
		this.freeBBZones = freeBBZones;
		this.cfp = cfp;
		this.initialUnPrunedConfs = countUnPrunedConfs();
		this.prunedConfs = countPrunedConfsByDEE();
		this.EW_I0 = EW_I0;

		topConfsPQ = new PriorityQueue<KSConf>(getNumTopConfsToSave(), 
				(new KSConf(null, 0.0, Double.MAX_VALUE)).new KSConfMinEComparator());
	}


	protected abstract void printHeader();
	
	
	public static void setNumTopConfsToSave( int n ) {
		numTopConfsToSave = n;
	}


	public static int getNumTopConfsToSave() {
		return numTopConfsToSave;
	}


	protected int getNumTopSavedConfs() {
		return topConfsPQ.size();
	}


	protected void saveTopConf(KSConf conf) {

		if(getNumTopSavedConfs() >= getNumTopConfsToSave()) {

			if(topConfsPQ.peek().getMinEnergy() > conf.getMinEnergy()) {
				topConfsPQ.poll();
			}

			else
				return;
		}

		topConfsPQ.add(conf);
	}


	protected void printTopConfs() {

		System.out.println("\nWriting top " + getNumTopSavedConfs() + 
				" conformation(s) for sequence: " + KSAbstract.arrayList1D2String(sequence, " "));
		System.out.println();

		// create dir if it does not already exist
		String dir = "topConfs" + File.separator + KSAbstract.arrayList1D2String(sequence, ".");
		ObjectIO.makeDir(dir, false);

		String pdbName = null;
		for( int i = getNumTopSavedConfs()-1; i > -1; i-- ) {
			System.out.println("Saving: " + i +".pdb" + "\tminE:" + topConfsPQ.peek().getMinEnergy());
			pdbName = dir + File.separator + String.valueOf(i) +".pdb";
			sp.outputMinimizedStruct(topConfsPQ.poll().getConf(), pdbName);
		}
		
		System.out.println();
	}


	public ArrayList<String> getSequence() {
		return sequence;
	}


	public EApproxReached getEpsilonStatus() {
		return eAppx;
	}
	
	
	public void setEpsilonStatus( EApproxReached in ) {
		eAppx = in;
	}


	protected boolean epsilonPossible() {

		if( eAppx == EApproxReached.FALSE || eAppx == EApproxReached.NOT_POSSIBLE ) {
			eAppx = EApproxReached.NOT_POSSIBLE;
			return false;
		}

		return true;
	}


	public void abort() {}


	public BigDecimal getQStar() {
		return qStar;
	}


	public BigDecimal getQPrime() {
		return qPrime;
	}


	public BigDecimal getPStar() {
		return pStar;
	}


	public BigDecimal getQDagger() {
		return null;
	}


	public BigDecimal getLowerBound() {
		return getQStar();
	}


	public BigDecimal getUpperBound() {
		if( eAppx == EApproxReached.TRUE ) return getQStar();

		updateQPrime();
		return qStar.add(qPrime);
	}


	public BigDecimal getUpperBoundAtEpsilon() {
		return getUpperBound();
	}


	public void setNumUnPrunedConfs() {
		initialUnPrunedConfs = countUnPrunedConfs();
	}


	private BigInteger countUnPrunedConfs() {
		if(pc == null) return BigInteger.ZERO;

		BigInteger numUPConfs = BigInteger.ONE;

		for( int pos = 0; pos < sp.confSpace.numPos; ++pos ) {
			numUPConfs = numUPConfs.multiply( BigInteger.valueOf( sp.pruneMat.unprunedRCsAtPos(pos).size() ) );
		}
		return numUPConfs;
	}


	public void setNumPrunedConfsByDEE() {
		prunedConfs = countPrunedConfsByDEE();
	}


	private BigInteger countPrunedConfsByDEE() {
		if(pc == null) return BigInteger.ZERO;

		BigInteger numPConfs = BigInteger.ONE;

		for( int pos = 0; pos < sp.confSpace.numPos; ++pos ) {
			numPConfs = numPConfs.multiply( BigInteger.valueOf( sp.pruneMat.prunedRCsAtPos(pos).size() ) );
		}

		return numPConfs;
	}


	protected double computeEffectiveEpsilon() {
		updateQPrime();

		BigDecimal divisor = (qStar.add(qPrime)).add(pStar);

		// divisor is 0 iff qstar = 0. this means the energies are too high
		// so epsilon can never be reached
		if( divisor.compareTo(BigDecimal.ZERO) == 0 ) return -1.0;

		BigDecimal maxQStar = qStar.add(qPrime);

		double minEpsilon = BigDecimal.ONE.subtract( maxQStar.divide(divisor, 4) ).doubleValue();

		if( minEpsilon > targetEpsilon ) return -1.0;

		return BigDecimal.ONE.subtract( qStar.divide(divisor, 4) ).doubleValue();
	}


	public double getEffectiveEpsilon() {
		return effectiveEpsilon;
	}


	protected void setPStar( double eLB ) {
		pStar = ( getBoltzmannWeight( eLB + EW_I0 )).multiply( new BigDecimal(prunedConfs) );
	}


	protected void updateQPrime() {
		qPrime = ( getBoltzmannWeight( Et )).multiply( new BigDecimal(getNumUnMinimizedConfs()) );
	}


	protected void updateQStar( KSConf conf ) {

		if(saveTopConfsAsPDB) {
			saveTopConf(conf);
		}

		qStar = qStar.add( getBoltzmannWeight( conf.getMinEnergy() ) );
		minimizedConfs = minimizedConfs.add(BigInteger.ONE);
		minimizedConfsDuringInterval = minimizedConfsDuringInterval.add(BigInteger.ONE);
	}


	public BigDecimal getBoltzmannWeight( double E ) {
		return e.exp(-E / RT);
	}


	public abstract void start();

	public void resume() {

		long startTime = System.currentTimeMillis();

		// setRunState( RunState.STARTED );

		while( eAppx == EApproxReached.FALSE && (double)(System.currentTimeMillis() - startTime) < interval ) {
			computeSlice();

			if( eAppx == EApproxReached.NOT_POSSIBLE ) {
				restart(100);
				computeSlice();
			}
		}

		if( eAppx == EApproxReached.TRUE && saveTopConfsAsPDB ) {
			printTopConfs();
		}

		if( eAppx != EApproxReached.FALSE ) {
			clearSearchProblem();
		}

		/*
		if( eAppx == EApproxReached.FALSE ) 
			setRunState( RunState.SUSPENDED );

		else 
			setRunState(RunState.TERMINATED);
		 */
	}

	public void runToCompletion() {
		compute();

		if( eAppx == EApproxReached.NOT_POSSIBLE ) {
			restart(100);
			compute();
		}
		//setRunState(RunState.TERMINATED);

		if( eAppx == EApproxReached.TRUE && saveTopConfsAsPDB ) {
			printTopConfs();
		}

		clearSearchProblem();
	}

	/**
	 * computes partition function for an interval of time then yields
	 */
	abstract protected void computeSlice();

	abstract protected void compute();

	abstract protected void iterate() throws Exception;

	protected void restart( double EW_I0 ) {

		this.EW_I0 = EW_I0;

		System.out.println("Could not reach target epsilon approximation of " + targetEpsilon + " for sequence: " +
				KSAbstract.arrayList1D2String(sequence, " "));

		System.out.println("Restarting with EW+I0 = " + EW_I0);

		pc = getPruningControl();

		pc.prune();

		initialUnPrunedConfs = countUnPrunedConfs();

		prunedConfs = countPrunedConfsByDEE();

		phaseOneMinimizedConfs = phaseOneMinimizedConfs.add(getNumMinimizedConfs());

		minimizedConfs = BigInteger.ZERO;

		qStar = BigDecimal.ZERO;

		eAppx = EApproxReached.FALSE;

		start();
	}


	protected void clearSearchProblem() {
		sp = null;
		pc = null;
	}


	public SearchProblem getSearchProblem() {
		return sp;
	}


	public PruningControl getPruningControl() {
		if(pc == null) pc = cfp.getPruningControl(sp, EW_I0, false, false);
		return pc;
	}


	protected BigInteger getNumPrunedConfs() {
		return prunedConfs;
	}


	protected BigInteger getNumUnMinimizedConfs() {
		return initialUnPrunedConfs.subtract(minimizedConfs);
	}


	public BigInteger getNumInitialUnPrunedConfs() {
		return initialUnPrunedConfs;
	}


	public BigInteger getNumMinimizedConfs() {
		return minimizedConfs;
	}


	public BigInteger getPhaseOneMinimizedConfs() {
		return phaseOneMinimizedConfs;
	}


	protected BigInteger getNumMinimizedConfsDuringInterval() {
		return minimizedConfsDuringInterval;
	}


	protected void resetNumMinimizedConfsDuringInterval() {
		minimizedConfsDuringInterval = BigInteger.ZERO;
	}


	private static int setNumPEs( int requested ) { 

		if(requested < 1) requested = 1;

		else if(requested > Runtime.getRuntime().availableProcessors()) 
			requested = Runtime.getRuntime().availableProcessors();

		return requested;
	}


	public static void setNumThreads( int threads ) { 
		numThreads = setNumPEs(threads);
	}


	public static int getNumThreads() {
		return numThreads;
	}


	public static void setNumFibers( int fibers ) {
		numFibers = setNumPEs(fibers);
	}


	public static int getNumFibers() {
		return numFibers;
	}


	public static int getConfsThreadBuffer() {
		return threadConfsBuffer;
	}


	public static void setConfsThreadBuffer( int confsBuffer ) {
		if( confsBuffer > 0 ) threadConfsBuffer = confsBuffer;
	}


	public static void setServerList( String[] servers ) {
		for( String server : servers ) 
			serverList.add( server.trim() );
	}


	public static void setNumRemoteClients( int clients ) {
		int nServers = serverList.size() == 0 ? 1 : serverList.size();

		numRemoteClients = clients > nServers ? nServers : clients;
	}


	public static int getNumRemoteClients() {
		return numRemoteClients;
	}


	public static void setInterval( String val ) {
		double dVal = -1;

		if(val.toLowerCase().equals("max")) {
			interval = Double.MAX_VALUE;
		} else {
			try{ dVal = Double.parseDouble(val); } catch( Exception e ) {}
			finally { if( dVal > 0 ) interval = dVal; }
		}
	}


	public static double getInterval() {
		return interval;
	}


	public static double getMaxInterval() {
		return Double.MAX_VALUE;
	}


	public void setRunState( RunState newState ) {
		runState = newState;
	}


	public RunState getRunState() {
		return runState;
	}


	public static String getImplementation() {
		return pFuncImplementation;
	}


	public static void setStabilityThreshold(double threshold) {
		threshold = threshold < 0.01 ? 1.0 : threshold;
		stabilityThresh = new BigDecimal(threshold);
	}


	public static BigDecimal getStabilityThreshold() {
		return stabilityThresh;
	}


	public static void setImplementation( String implementation ) {

		switch( implementation.toLowerCase() ) {

		case "1nnocache":
		case "1npmcache":
		case "1npcpmcache":
		case "1nmtpcpmcache":
		case "mnpcpmcache":
			pFuncImplementation = implementation;
			break;

		default:
			throw new RuntimeException("ERROR: specified value of parameter pFuncMethod is invalid");
		}
	}

}