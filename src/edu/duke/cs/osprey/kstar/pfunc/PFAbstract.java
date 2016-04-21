package edu.duke.cs.osprey.kstar.pfunc;

import java.io.File;
import java.io.Serializable;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashSet;
import java.util.PriorityQueue;

import edu.duke.cs.osprey.astar.ConfTree;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.kstar.AllowedSeqs;
import edu.duke.cs.osprey.kstar.KAStarConfTree;
import edu.duke.cs.osprey.kstar.KSAbstract;
import edu.duke.cs.osprey.kstar.KSConf;
import edu.duke.cs.osprey.kstar.PStarConfTree;
import edu.duke.cs.osprey.kstar.pfunc.impl.PFTrad;
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
	protected static String pFuncCFGImpl = new PFTrad().getImpl();
	public static String eMinMethod = "ccd";
	protected static ArrayList<String> serverList = new ArrayList<>();
	protected static int threadConfsBuffer = 8;
	public static int qCapacity = 1024000;
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
	protected String searchProblemName = null;

	protected static BigDecimal stabilityThresh = BigDecimal.ONE;

	protected static final double RT = 1.9891/1000.0 * 298.15;
	public static double targetEpsilon = 0.03;
	protected double effectiveEpsilon = 1.0;

	protected int strand = -1;
	protected ConfigFileParser cfp = null;
	protected SearchProblem sp = null;
	protected ArrayList<Integer> flexResIndexes;
	protected SearchProblem panSeqSP = null;
	protected PruningControl pc = null;
	private boolean isFullyDefined = true;

	protected BigDecimal qStar = BigDecimal.ZERO;
	protected BigDecimal qPrime = BigDecimal.ZERO;
	protected BigDecimal pStar = BigDecimal.ZERO;

	protected ExpFunction e = new ExpFunction();
	protected double Et = 0;
	protected double E0 = 0;

	protected BigInteger prunedConfs = BigInteger.ZERO;
	protected BigInteger unPrunedConfs = BigInteger.ZERO;
	protected BigInteger minimizedConfs = BigInteger.ZERO;
	protected HashSet<ArrayList<Integer>> minimizedConfsSet = new HashSet<>();
	protected BigInteger minimizedConfsDuringInterval = BigInteger.ZERO;
	protected BigInteger minimizingConfs = BigInteger.ZERO; // # confs being minimized at this instant

	protected PriorityQueue<KSConf> topConfsPQ = null;

	protected PFAbstract() {}

	protected PFAbstract( int strand, ArrayList<String> sequence, ArrayList<Integer> flexResIndexes, 
			String checkPointPath, String searchProblemName, 
			ConfigFileParser cfp, SearchProblem panSeqSP ) {

		this.sequence = sequence;
		this.flexResIndexes = flexResIndexes;
		this.checkPointPath = checkPointPath;
		this.searchProblemName = searchProblemName;
		this.panSeqSP = panSeqSP;
		this.strand = strand;
		this.sp = createSingleSeqSP(panSeqSP.contSCFlex, strand, sequence, flexResIndexes, true);
		
		this.isFullyDefined = sp.confSpace.numPos == panSeqSP.confSpace.numPos ? true : false;
		this.cfp = cfp;
		
		// re-prune, since we have fewer witnesses now that we have trimmed the emat
		double oldPI = sp.pruneMat.getPruningInterval();
		sp.pruneMat.setPruningInterval(0);
		pc = getPruningControl(oldPI); pc.prune(); pc = null;
		sp.pruneMat.setPruningInterval(oldPI);

		Comparator<KSConf> comparator = new KSConf(new ArrayList<>(), 0.0).new KSConfMinEComparator();
		topConfsPQ = new PriorityQueue<KSConf>(getNumTopConfsToSave(), comparator);
	}


	public boolean isFullyDefined() {
		return isFullyDefined;
	}


	protected BigDecimal productUndefinedRots() {

		BigDecimal ans = BigDecimal.ONE;

		// get unassigned residue positions
		int numPos = panSeqSP.confSpace.numPos;
		ArrayList<Integer> undefinedPos = new ArrayList<>(numPos);
		for(int pos = 0; pos < numPos; ++pos) undefinedPos.add(pos);

		undefinedPos.removeAll(flexResIndexes);

		boolean minimizeProduct = panSeqSP.contSCFlex ? false : true;

		for( int pos : undefinedPos ) {

			int rcNumAtPos = minimizeProduct ? Integer.MAX_VALUE : Integer.MIN_VALUE;

			for( String AAType : panSeqSP.allowedAAs.get(pos) ) {
				// get number of unpruned rcs at that level for that AAType
				int num = panSeqSP.pruneMat.getNumRCsAtPosForAA(panSeqSP.confSpace, pos, AAType, false);

				if(!minimizeProduct)
					num += panSeqSP.pruneMat.getNumRCsAtPosForAA(panSeqSP.confSpace, pos, AAType, true);

				rcNumAtPos = minimizeProduct ? Math.min(rcNumAtPos, num) : Math.max(rcNumAtPos, num);
			}

			ans = ans.multiply(BigDecimal.valueOf(rcNumAtPos));
		}

		return ans;
	}


	protected void adjustQStar() {
		if(eAppx != EApproxReached.TRUE) return;

		// no undefined positions
		if(isFullyDefined()) return;

		// don't want to zero out our qstar
		BigDecimal product = productUndefinedRots();
		if(product.compareTo(BigDecimal.ZERO) == 0)
			return;

		qStar = qStar.multiply(product);
	}


	public ConfSearch getConfTree( boolean usePrunedConfs ) {
		if(isFullyDefined()) {
			if( usePrunedConfs && canUseTightPStar() ) 
				return new PStarConfTree(sp);

			return new ConfTree(sp);
		}

		else {
			if( usePrunedConfs && !canUseTightPStar() ) usePrunedConfs = false;
			return new KAStarConfTree(sp, panSeqSP, flexResIndexes, usePrunedConfs);
		}
	}


	public double getConfBound( ConfSearch confSearch, int[] conf, boolean usePrunedConfs ) {
		double bound = 0;

		if( isFullyDefined() ) {
			if( usePrunedConfs && canUseTightPStar() )
				bound = ((PStarConfTree)confSearch).confBound(conf);

			else bound = sp.lowerBound(conf);
		}

		else
			bound = ((KAStarConfTree)confSearch).confBound(conf);

		return bound;
	}


	protected boolean canUseTightPStar() {

		boolean singleSeqPrunedConfsExist = sp.numConfs(true).compareTo(BigInteger.ZERO) > 0;

		if( isFullyDefined() )
			return singleSeqPrunedConfsExist;

		boolean panSeqPrunedConfsExist = panSeqSP.numConfs(true).compareTo(BigInteger.ZERO) > 0;

		return singleSeqPrunedConfsExist && panSeqPrunedConfsExist;
	}


	public String getSearchProblemName() {
		return searchProblemName;
	}


	public int getStrand() {
		return strand;
	}


	public void setPanSeqSP( SearchProblem in ) {
		panSeqSP = in;
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

			if(topConfsPQ.peek().getEnergy() > conf.getEnergy()) {
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
			System.out.println("Saving: " + i +".pdb" + "\tminE:" + tmp.peek().getEnergy());
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
		return BigDecimal.ZERO;
	}


	public BigDecimal getQStarLowerBound() {
		return getQStar();
	}


	public BigDecimal getQStarUpperBound() {
		if( eAppx == EApproxReached.TRUE ) return getQStar();

		updateQPrime();
		return (qStar.add(qPrime)).add(pStar);
	}


	public void setNumUnPruned() {
		unPrunedConfs = sp.numConfs(false);
	}


	public void setNumPruned() {
		prunedConfs = sp.numConfs(true);
	}


	protected double computeEffectiveEpsilon() {

		BigDecimal divisor = qStar.add(qPrime.add(pStar));

		// energies are too high so epsilon can never be reached
		if( divisor.compareTo(BigDecimal.ZERO) == 0 ) return -1.0;

		return BigDecimal.ONE.subtract( qStar.divide(divisor, 4) ).doubleValue();
	}


	public double getEffectiveEpsilon() {
		return effectiveEpsilon;
	}


	protected void initPStar() {

		ConfSearch confSearch = getConfTree(true);
		int conf[];

		if( (conf = confSearch.nextConf()) != null ) {
			setPStar( getConfBound(confSearch, conf, true), true );
		}

		else
			throw new RuntimeException("ERROR: cannot set P*");
	}


	private void setPStar( double eLB, boolean usePrunedConfs ) {

		if( usePrunedConfs && canUseTightPStar() ) {
			setPStarTight(eLB);
		}

		else
			setPStarLoose(eLB);
	}


	private void setPStarLoose( double eLB ) {
		E0 = eLB + getSearchProblem().pruneMat.getPruningInterval();
		pStar = ( getBoltzmannWeight( E0 )).multiply( new BigDecimal(prunedConfs) );
	}


	private void setPStarTight( double eLB ) {
		E0 = eLB;
		pStar = ( getBoltzmannWeight( E0 )).multiply( new BigDecimal(prunedConfs) );
	}


	protected void updateQPrime() {
		qPrime = getBoltzmannWeight( Et ).multiply(new BigDecimal(getNumUnEnumerated()));
	}


	protected void updateQStar( KSConf conf ) {

		if(saveTopConfsAsPDB) {
			saveTopConf(conf);
		}

		qStar = qStar.add( getBoltzmannWeight( conf.getEnergy() ) );

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
				phase2();

				if( eAppx == EApproxReached.NOT_POSSIBLE ) break;

				computeSlice();
			}
		}

		if( saveTopConfsAsPDB && eAppx == EApproxReached.TRUE ) writeTopConfs();

		if( eAppx != EApproxReached.FALSE ) cleanup();

		resetMinDuringInterval();
	}


	public void runToCompletion() {
		compute();

		if( eAppx == EApproxReached.NOT_POSSIBLE ) {
			phase2();

			if( eAppx == EApproxReached.FALSE ) 
				compute();
		}

		if( saveTopConfsAsPDB && eAppx == EApproxReached.TRUE ) writeTopConfs();

		cleanup();
	}

	/**
	 * computes partition function for an interval of time then yields
	 */
	abstract protected void computeSlice();

	abstract protected void compute();

	abstract protected void iterate() throws Exception;

	protected void restart() {
		eAppx = EApproxReached.FALSE;
		printedHeader = false;
		restarted = true;
		minimizedConfs = BigInteger.ZERO;
		minimizingConfs = BigInteger.ZERO;
		start();
	}

	protected void phase2() {

		/*
		 * TODO
		 * 1) what needs to happen for partially defined sequences?
		 */

		System.out.println("\nCould not reach target epsilon approximation of " + targetEpsilon + " for sequence: " +
				KSAbstract.list1D2String(sequence, " "));

		if( getEffectiveEpsilon() < 0 ) {
			// we can never reach epsilon because q* + q' + p* = 0
			System.out.println("\nCan never reach target epsilon approximation of " + targetEpsilon + " for sequence: " +
					KSAbstract.list1D2String(sequence, " "));
			return;
		}

		System.out.println("Attempting Phase 2...");

		// completely relax pruning
		double maxPruningInterval = 100;
		pc = getPruningControl(maxPruningInterval); pc.prune(); pc = null;
		sp.pruneMat.setPruningInterval(maxPruningInterval);

		setNumUnPruned();
		setNumPruned(); // needed for p*

		/*
		// conservative implementation
		restart();
		return;
		 */

		// shortcut implementation
		// get new value of epsilon with q' = 0 and reduced p*
		// MUST call abstract base class version of this method

		qPrime = BigDecimal.ZERO;

		initPStar(); // re-calculate p*

		double effectiveEpsilon = computeEffectiveEpsilon();

		if( getQStar().compareTo(BigDecimal.ZERO) > 0 
				&& effectiveEpsilon != -1.0 && effectiveEpsilon <= targetEpsilon ) {

			setEpsilonStatus(EApproxReached.TRUE);

			System.out.println("\nReached target epsilon approximation of " + targetEpsilon + " for sequence: " +
					KSAbstract.list1D2String(sequence, " "));
			return;
		}

		else {
			System.out.println("\nRe-starting K* for sequence" + KSAbstract.list1D2String(sequence, " "));
			restart();
		}
	}


	public void cleanup() {
		panSeqSP = null;
		sp = null;
		pc = null;
		minimizedConfsSet = null;
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


	protected SearchProblem createSingleSeqSP( boolean contSCFlex, int strand, 
			ArrayList<String> seq, ArrayList<Integer> flexResIndexes, boolean fast ) {

		if(fast) return createSingleSeqSPFast(contSCFlex, strand, seq, flexResIndexes);

		return createSingleSeqSPSlow(contSCFlex, strand, seq);
	}


	private SearchProblem createSingleSeqSPSlow( boolean contSCFlex, int strand, ArrayList<String> seq ) {

		ArrayList<ArrayList<String>> allowedAAs = KSAbstract.list1D2ListOfLists(AllowedSeqs.getAAsFromSeq(seq));
		ArrayList<String> flexibleRes = AllowedSeqs.getFlexResFromSeq(seq);

		// create searchproblem
		SearchProblem seqSP = new SearchProblem( 
				searchProblemName, 
				panSeqSP.PDBFile, 
				flexibleRes, 
				allowedAAs, 
				false, 
				contSCFlex,
				panSeqSP.useEPIC,
				panSeqSP.epicSettings,
				panSeqSP.useTupExpForSearch,
				panSeqSP.dset, 
				panSeqSP.moveableStrands, 
				panSeqSP.freeBBZones,
				panSeqSP.useEllipses,
				panSeqSP.useERef,
				panSeqSP.addResEntropy,
				panSeqSP.limits);

		return seqSP;
	}


	private SearchProblem createSingleSeqSPFast( boolean contSCFlex, int strand, 
			ArrayList<String> seq, ArrayList<Integer> flexResIndexes ) {

		SearchProblem seqSP = panSeqSP.singleSeqSearchProblem(searchProblemName, 
				AllowedSeqs.getAAsFromSeq(seq), AllowedSeqs.getFlexResFromSeq(seq), 
				flexResIndexes);

		return seqSP;
	}

}