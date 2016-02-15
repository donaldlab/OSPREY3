package edu.duke.cs.osprey.kstar;

import java.math.BigDecimal;
import java.math.BigInteger;
import java.util.ArrayList;

import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.ematrix.epic.EPICSettings;
import edu.duke.cs.osprey.pruning.PruningControl;
import edu.duke.cs.osprey.tools.ExpFunction;

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
	public static boolean waitUntilCapacity = false;
	protected static int numThreads = 1;
	protected static int numFibers = 1;
	protected static int numRemoteClients = 1;

	protected static BigDecimal stabilityThresh = BigDecimal.ONE;
	protected static double interval = 0.5;

	protected static final double RT = 1.9891/1000.0 * 298.15;
	public static double targetEpsilon = 0.03;
	public static double rho = targetEpsilon / (1.0 - targetEpsilon);
	protected double effectiveEpsilon = 1.0;

	public static enum EApproxReached { TRUE, FALSE, NOT_POSSIBLE, NOT_STABLE, ABORTED }
	protected EApproxReached eAppx = EApproxReached.FALSE;

	public static enum RunState { NOTSTARTED, STARTED, SUSPENDED, TERMINATED }
	protected RunState runState = RunState.NOTSTARTED;

	protected ConfigFileParser cfp;
	protected SearchProblem sp;
	protected PruningControl pc;
	protected DEEPerSettings dset;
	protected ArrayList<String[]> moveableStrands;
	protected ArrayList<String[]> freeBBZones;

	protected BigDecimal qStar = BigDecimal.ZERO;
	protected BigDecimal qPrime = BigDecimal.ZERO;
	protected BigDecimal pStar = BigDecimal.ZERO;

	protected ExpFunction e = new ExpFunction();
	protected double EW_I0 = 5;
	protected double Et = 0;

	protected BigInteger prunedConfs = BigInteger.ZERO;
	protected BigInteger initialUnPrunedConfs = BigInteger.ZERO;
	protected BigInteger minimizedConfs = BigInteger.ZERO;
	protected BigInteger minimizedConfsDuringInterval = BigInteger.ZERO;
	protected BigInteger minimizingConfs = BigInteger.ZERO; // # confs being minimized at this instant


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
	}


	public ArrayList<String> getSequence() {
		return sequence;
	}


	public EApproxReached getEApproxReached() {
		return eAppx;
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


	private BigInteger countUnPrunedConfs() {
		BigInteger numUPConfs = BigInteger.ONE;

		for( int pos = 0; pos < sp.confSpace.numPos; ++pos ) {
			numUPConfs = numUPConfs.multiply( BigInteger.valueOf( sp.pruneMat.unprunedRCsAtPos(pos).size() ) );
		}
		return numUPConfs;
	}


	private BigInteger countPrunedConfsByDEE() {
		BigInteger numPConfs = BigInteger.ONE;

		for( int pos = 0; pos < sp.confSpace.numPos; ++pos ) {
			numPConfs = numPConfs.multiply( BigInteger.valueOf( sp.pruneMat.prunedRCsAtPos(pos).size() ) );
		}

		return numPConfs;
	}


	protected double getStopThreshold() {
		return BigDecimal.valueOf(-RT).multiply
				( ( e.log
						( ( qStar.multiply
								( BigDecimal.valueOf(rho) ) ).subtract
								( pStar ) ) ).subtract
						( e.log
								( new BigDecimal(getNumUnMinimizedConfs()) ) ) ).doubleValue();
	}


	protected double computeEffectiveEpsilon() {
		updateQPrime();

		BigDecimal divisor = (qStar.add(qPrime)).add(pStar);

		// divisor is 0 iff qstar = 0. this means the energies are too high
		// so epsilon can never be reached
		if( divisor.compareTo(BigDecimal.ZERO) == 0 ) return -1.0;

		BigDecimal maxQStar = qStar.add(qPrime);

		BigDecimal minEpsilon = BigDecimal.ONE.subtract( maxQStar.divide(divisor, 4) );

		if( minEpsilon.compareTo(BigDecimal.valueOf(targetEpsilon)) > 0 ) return -1.0;

		return BigDecimal.ONE.subtract( qStar.divide(divisor, 4) ).doubleValue();
	}


	protected double getEffectiveEpsilon() {
		return effectiveEpsilon;
	}


	protected void setPStar( double eLB ) {
		pStar = ( getBoltzmannWeight( eLB + EW_I0 )).multiply( new BigDecimal(prunedConfs) );
	}


	protected void updateQPrime() {
		qPrime = ( getBoltzmannWeight( Et )).multiply( new BigDecimal(getNumUnMinimizedConfs()) );
	}


	protected void updateQStar( double E ) {
		qStar =  qStar.add( getBoltzmannWeight( E ) );
		minimizedConfs = minimizedConfs.add(BigInteger.ONE);
		minimizedConfsDuringInterval = minimizedConfsDuringInterval.add(BigInteger.ONE);
	}


	protected BigDecimal getBoltzmannWeight( double E ) {
		return e.exp(-E / RT);
	}


	abstract protected void start();

	protected void resume() {

		long startTime = System.currentTimeMillis();

		// setRunState( RunState.STARTED );

		while( eAppx == EApproxReached.FALSE && (double)(System.currentTimeMillis() - startTime) < interval ) {
			computeSlice();
			
			if( eAppx == EApproxReached.NOT_POSSIBLE ) {
				restart(100);
				computeSlice();
			}
		}
		/*
		if( eAppx == EApproxReached.FALSE ) 
			setRunState( RunState.SUSPENDED );

		else 
			setRunState(RunState.TERMINATED);
		 */
	}

	protected void runToCompletion() {
		compute();

		if( eAppx == EApproxReached.NOT_POSSIBLE ) {
			restart(100);
			compute();
		}
		//setRunState(RunState.TERMINATED);
	}

	/**
	 * computes partition function for an interval of time then yields
	 */
	abstract protected void computeSlice();

	abstract protected void compute();

	abstract protected void iterate() throws Exception;

	protected void restart( double EW_I0 ) {

		System.out.print("Could not reach target epsilon approximation of " + targetEpsilon + " for sequence: ");
		KSCalc.print(sequence, System.out);
		System.out.println();
		System.out.println("Restarting with EW+I0 = " + EW_I0);

		this.EW_I0 = EW_I0;
		
		sp = new SearchProblem( sp.name, sp.PDBFile, sp.flexibleRes, sp.allowedAAs, false, true,
				cfp.getParams().getBool("UseEPIC"),
                new EPICSettings(cfp.getParams()),
                cfp.getParams().getBool("UseTupExp"),
                dset, moveableStrands, freeBBZones,
                cfp.getParams().getBool("useEllipses"),
                cfp.getParams().getBool("useERef"),
                cfp.getParams().getBool("AddResEntropy"));

		sp.loadEnergyMatrix();

		pc = cfp.getPruningControl(sp, EW_I0, false, false);

		pc.prune();

		initialUnPrunedConfs = countUnPrunedConfs();

		prunedConfs = countPrunedConfsByDEE();
		
		qStar = BigDecimal.ZERO;

		eAppx = EApproxReached.FALSE;
		
		start();
	}

	protected EApproxReached accumulate( int conf[] ) {
		return EApproxReached.FALSE;
	}


	protected BigInteger getNumPrunedConfs() {
		return prunedConfs;
	}


	protected BigInteger getNumUnMinimizedConfs() {
		return initialUnPrunedConfs.subtract(minimizedConfs);
	}


	protected BigInteger getNumInitialUnPrunedConfs() {
		return initialUnPrunedConfs;
	}


	protected BigInteger getNumMinimizedConfs() {
		return minimizedConfs;
	}


	protected BigInteger getNumMinimizedConfsDuringInterval() {
		return minimizedConfsDuringInterval;
	}


	protected void resetNumMinimizedConfsDuringInterval() {
		minimizedConfsDuringInterval = BigInteger.ZERO;
	}


	private static int setNumPEs( int requested ) { 

		if(requested < 1) requested = 1;

		else if(requested > Runtime.getRuntime().availableProcessors() - 1) 
			requested = Runtime.getRuntime().availableProcessors() - 1;

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
		case "1nastar":
		case "1npmcache":
		case "1npcpmcache":
		case "1nmtpcpmcache":
		case "mnpcpmcache":
		case "mnmcpcpmcache":
			pFuncImplementation = implementation;
			break;

		default:
			throw new RuntimeException("ERROR: specified value of parameter pFuncMethod is invalid");
		}
	}

}