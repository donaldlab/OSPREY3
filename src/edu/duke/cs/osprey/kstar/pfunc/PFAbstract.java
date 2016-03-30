package edu.duke.cs.osprey.kstar.pfunc;

import java.io.File;
import java.io.Serializable;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashSet;
import java.util.PriorityQueue;

import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.ConfigFileParser;
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
@SuppressWarnings("serial")
public abstract class PFAbstract implements Serializable {

	public static enum EApproxReached { TRUE, FALSE, NOT_POSSIBLE, NOT_STABLE }
	protected EApproxReached eAppx = EApproxReached.FALSE;

	public static enum RunState { NOTSTARTED, STARTED }
	protected RunState runState = RunState.NOTSTARTED;
	
	protected ArrayList<String> sequence;
	protected static String pFuncCFGImpl = "trad";
	public static String eMinMethod = "ccd";
	protected static ArrayList<String> serverList = new ArrayList<>();
	protected static int threadConfsBuffer = 8;
	public static int qCapacity = 4194304;
	public static boolean waitUntilCapacity = false;
	protected static int numThreads = 1;
	protected static int numFibers = 1;
	protected static int numRemoteClients = 1;

	public static boolean suppressOutput = false;
	protected boolean printedHeader = false;
	protected boolean restarted = false;

	public static boolean saveTopConfsAsPDB = false;
	protected static int numTopConfsToSave = 10;

	public static boolean useMaxKSConfs = false;
	protected static long maxKSConfs = 100000;

	protected String checkPointPath = null;

	protected static BigDecimal stabilityThresh = BigDecimal.ONE;

	protected static final double RT = 1.9891/1000.0 * 298.15;
	public static double targetEpsilon = 0.03;
	protected double effectiveEpsilon = 1.0;

	protected ConfigFileParser cfp = null;
	protected SearchProblem sp = null;
	protected PruningControl pc = null;

	protected BigDecimal qStar = BigDecimal.ZERO;
	protected BigDecimal qPrime = BigDecimal.ZERO;
	protected BigDecimal pStar = BigDecimal.ZERO;

	protected ExpFunction e = new ExpFunction();
	protected double EW_I0 = 5;
	protected double Et = 0;
	protected double E0 = 0;

	protected BigInteger prunedConfs = BigInteger.ZERO;
	protected BigInteger unPrunedConfs = BigInteger.ZERO;
	protected BigInteger minimizedConfs = BigInteger.ZERO;
	protected HashSet<ArrayList<Integer>> minimizedConfsSet = new HashSet<>();
	protected BigInteger minimizedConfsDuringInterval = BigInteger.ZERO;
	protected BigInteger minimizingConfs = BigInteger.ZERO; // # confs being minimized at this instant

	protected PriorityQueue<KSConf> topConfsPQ = null;

	protected PFAbstract( ArrayList<String> sequence, String checkPointPath, 
			ConfigFileParser cfp, SearchProblem sp, double EW_I0 ) {

		this.sequence = sequence;
		this.checkPointPath = checkPointPath;
		this.sp = sp;
		this.cfp = cfp;
		this.EW_I0 = EW_I0;
		
		Comparator<KSConf> comparator = new KSConf(new ArrayList<>(), 0.0).new KSConfComparator();
		topConfsPQ = new PriorityQueue<KSConf>(getNumTopConfsToSave(), comparator);
	}
	

	public HashSet<ArrayList<Integer>> getMinimizedConfsSet() {
		return minimizedConfsSet;
	}


	protected abstract void printHeader();


	public static void setNumTopConfsToSave( int n ) {
		numTopConfsToSave = n;
	}


	public static int getNumTopConfsToSave() {
		return numTopConfsToSave;
	}


	public int getNumTopSavedConfs() {
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


	public void writeTopConfs() {

		if( getNumTopSavedConfs() == 0 ) return;
		
		System.out.println("\nWriting top " + getNumTopSavedConfs() + 
				" conformation(s) for sequence: " + KSAbstract.list1D2String(sequence, " "));
		System.out.println();

		@SuppressWarnings("unchecked")
		PriorityQueue<KSConf> tmp = (PriorityQueue<KSConf>) ObjectIO.deepCopy(topConfsPQ);
		
		// create dir if it does not already exist
		String dir = "topConfs" + File.separator + KSAbstract.list1D2String(sequence, ".");
		ObjectIO.makeDir(dir, false);

		String pdbName = null;
		for( int i = getNumTopSavedConfs()-1; i > -1; i-- ) {
			System.out.println("Saving: " + i +".pdb" + "\tminE:" + tmp.peek().getMinEnergy());
			pdbName = dir + File.separator + String.valueOf(i) +".pdb";
			sp.outputMinimizedStruct(tmp.poll().getConfArray(), pdbName);
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


	public void abort(boolean nullify) {}


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


	public void setNumUnPruned() {
		unPrunedConfs = sp.numUnPruned();
	}


	public void setNumPruned() {
		prunedConfs = sp.numPruned();
	}


	protected double computeEffectiveEpsilon() {
		updateQPrime();

		BigDecimal divisor = (qStar.add(qPrime)).add(pStar);

		// energies are too high so epsilon can never be reached
		if( divisor.compareTo(BigDecimal.ZERO) == 0 ) return -1.0;

		return BigDecimal.ONE.subtract( qStar.divide(divisor, 4) ).doubleValue();
	}


	public double getEffectiveEpsilon() {
		return effectiveEpsilon;
	}


	protected void setPStar( double eLB ) {
		E0 = eLB + EW_I0;
		pStar = ( getBoltzmannWeight( E0 )).multiply( new BigDecimal(prunedConfs) );
	}


	protected void updateQPrime() {
		qPrime = getBoltzmannWeight( Et ).multiply(new BigDecimal(getNumUnEnumerated()));
	}


	protected void updateQStar( KSConf conf ) {

		if(saveTopConfsAsPDB) {
			saveTopConf(conf);
		}

		qStar = qStar.add( getBoltzmannWeight( conf.getMinEnergy() ) );

		minimizedConfs = minimizedConfs.add(BigInteger.ONE);
		minimizedConfsSet.add(conf.getConf());
		minimizedConfsDuringInterval = minimizedConfsDuringInterval.add(BigInteger.ONE);
	}


	public BigDecimal getBoltzmannWeight( double E ) {
		return e.exp(-E / RT);
	}


	public abstract void start();

	public void runSlice(long target) {

		while( eAppx == EApproxReached.FALSE && getMinDuringInterval().longValue() < target ) {
			computeSlice();

			if( eAppx == EApproxReached.NOT_POSSIBLE ) {
				restart();
				
				if( eAppx == EApproxReached.NOT_POSSIBLE ) break;
				
				computeSlice();
			}
		}

		if( eAppx == EApproxReached.TRUE && saveTopConfsAsPDB ) writeTopConfs();

		if( eAppx != EApproxReached.FALSE ) cleanup();
		
		resetMinDuringInterval();
	}
	

	public void runToCompletion() {
		compute();

		if( eAppx == EApproxReached.NOT_POSSIBLE ) {
			restart();
			
			if( eAppx == EApproxReached.FALSE ) compute();
		}

		if( eAppx == EApproxReached.TRUE && saveTopConfsAsPDB ) writeTopConfs();

		cleanup();
	}

	/**
	 * computes partition function for an interval of time then yields
	 */
	abstract protected void computeSlice();

	abstract protected void compute();

	abstract protected void iterate() throws Exception;

	protected void restart() {

		System.out.println("\nCould not reach target epsilon approximation of " + targetEpsilon + " for sequence: " +
				KSAbstract.list1D2String(sequence, " "));

		BigDecimal rho = BigDecimal.valueOf(targetEpsilon/(1-targetEpsilon));
		BigDecimal bE0 = getBoltzmannWeight(E0);

		double pruningInterval = sp.pruneMat.getPruningInterval();
		
		if( bE0.compareTo(BigDecimal.ZERO) == 0 ) {
			pruningInterval = 100.0;
			pc = getPruningControl(pruningInterval); pc.prune();
			sp.pruneMat.setPruningInterval(pruningInterval);
		}
		
		else {

			BigDecimal l = new BigDecimal(prunedConfs).subtract( (qStar.multiply(rho)).divide(bE0, 4) );
			BigInteger li = l.add(BigDecimal.ONE).toBigInteger();
			BigInteger pruningTarget = prunedConfs.subtract(li);

			System.out.println("\nOld pruning window: " + pruningInterval);
			System.out.println("Number of pruned confs.: " + prunedConfs);
			System.out.println("Pruning target: " + pruningTarget + " confs.");

			// prune until we reach pruning target
			BigInteger numPruned = BigInteger.ZERO;

			do {
				pruningInterval = Math.min(pruningInterval + 2.5, 100.0);

				pc = getPruningControl(pruningInterval);
				pc.prune();

				numPruned = sp.numPruned();
				sp.pruneMat.setPruningInterval(pruningInterval);

				System.out.println("New pruning window: " + pruningInterval);
				System.out.println("Pruning target: " + pruningTarget + " confs.");
				System.out.println("Number of pruned confs.: " + numPruned + "\n");

			} while( numPruned.compareTo(pruningTarget) > 0 && pruningInterval < 100.0 );

			if( numPruned.compareTo(pruningTarget) > 0 ) {
				eAppx = EApproxReached.NOT_POSSIBLE;
				return;
			}
		}

		unPrunedConfs = sp.numUnPruned();

		prunedConfs = sp.numPruned();

		eAppx = EApproxReached.FALSE;
		
		printedHeader = false;
		restarted = true;
		
		// i should not need to do the next three lines. ask mark for help
		//minimizedConfsSet.clear();
		//qStar = BigDecimal.ZERO;
		minimizedConfs = BigInteger.ZERO;

		start();
	}


	public void cleanup() {
		sp = null;
		pc = null;
		minimizedConfsSet.clear();
	}


	public SearchProblem getSearchProblem() {
		return sp;
	}


	public PruningControl getPruningControl(double pruningInterval) {
		return cfp.getPruningControl(sp, pruningInterval, sp.useEPIC, sp.useTupExpForSearch);
	}


	protected BigInteger getNumPruned() {
		return prunedConfs;
	}


	protected BigInteger getNumUnEnumerated() {
		return unPrunedConfs.subtract(minimizedConfs);
	}


	public BigInteger getNumUnPruned() {
		return unPrunedConfs;
	}


	public BigInteger getNumMinimized() {
		return minimizedConfs;
	}
	
	
	// ugly hack to count number of minimized confs after a re-start occurs
	// used only for output purposes
	public BigInteger getNumMinimized4Output() {
		BigInteger numMinimized = getNumMinimized();
		if(restarted) numMinimized = numMinimized.add(BigInteger.valueOf(getMinimizedConfsSet().size()));
		return numMinimized;
	}


	protected BigInteger getMinDuringInterval() {
		return minimizedConfsDuringInterval;
	}


	protected void resetMinDuringInterval() {
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


	public static double getMaxInterval() {
		return Double.MAX_VALUE;
	}


	public void setRunState( RunState newState ) {
		runState = newState;
	}


	public RunState getRunState() {
		return runState;
	}


	public static void setStabilityThresh(double threshold) {
		threshold = threshold < 0.0 ? 0.0 : threshold;
		stabilityThresh = new BigDecimal(threshold);
	}


	public static BigDecimal getStabilityThresh() {
		return stabilityThresh;
	}


	public static void setMaxKSconfs( long in ) {
		if( in < 1 ) in = 1;
		maxKSConfs = in;
	}


	public boolean maxKSConfsReached() {
		return useMaxKSConfs && minimizedConfs.longValue() >= maxKSConfs;
	}

	
	public String getCheckPointPath() {
		return checkPointPath;
	}
	

	public boolean checkPointExists() {
		return new File(getCheckPointPath()).exists();
	}
	
	
	public abstract String getImpl();
	
	
	public static String getCFGImpl() {
		return pFuncCFGImpl;
	}
	
	
	public static void setCFGImpl( String implementation ) {

		switch( implementation.toLowerCase() ) {

		case "trad":
		case "new00":
		case "new01":
		case "new02":
		case "new03":
		case "new04":
			pFuncCFGImpl = implementation;
			break;

		default:
			throw new RuntimeException("ERROR: specified value of parameter pFuncMethod is invalid");
		}
	}

}