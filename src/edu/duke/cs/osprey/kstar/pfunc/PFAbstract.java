package edu.duke.cs.osprey.kstar.pfunc;

import java.io.File;
import java.io.Serializable;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashSet;
import java.util.PriorityQueue;

import edu.duke.cs.osprey.astar.ConfTree;
import edu.duke.cs.osprey.astar.FullAStarNode;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.kstar.KAStarConfTree;
import edu.duke.cs.osprey.kstar.KSAbstract;
import edu.duke.cs.osprey.kstar.KSAllowedSeqs;
import edu.duke.cs.osprey.kstar.KSConf;
import edu.duke.cs.osprey.kstar.KSConfigFileParser;
import edu.duke.cs.osprey.kstar.KSSearchProblem;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tools.ExpFunction;
import edu.duke.cs.osprey.tools.ObjectIO;

/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 *
 */
@SuppressWarnings("serial")
public abstract class PFAbstract implements Serializable {

	protected long startTime;
	
	protected DecimalFormat numberFormat = new DecimalFormat("0.0000");
	
	public static enum EApproxReached { TRUE, FALSE, NOT_POSSIBLE, NOT_STABLE }
	protected EApproxReached eAppx = EApproxReached.FALSE;

	protected final double EPSILON_NEVER_POSSIBLE = -1.0;
	protected final double EPSILON_PHASE_2 = -2.0;

	public static enum RunState { NOTSTARTED, STARTED }
	protected RunState runState = RunState.NOTSTARTED;

	protected ArrayList<String> sequence;
	protected static String pFuncCFGImpl = null;
	protected static int threadConfsBuffer = 8;
	public static int qCapacity = 1024000;
	protected static int numThreads = 1;

	public static boolean suppressOutput = false;
	protected boolean printedHeader = false;
	protected boolean doingPhase2 = false;

	public static boolean saveTopConfsAsPDB = false;
	protected static int numTopConfsToSave = 10;

	public static boolean useMaxKSConfs = false;
	protected static long maxKSConfs = 100000;

	protected static String hotMethod = "none"; 
	protected static double hotBoundPct = 0.2;
	protected static int hotNumRes = 3;
	protected static double hotTopRotsPct = 0.2;
	protected ArrayList<ArrayList<Integer>> HOTs = null;

	protected String checkPointPath = null;
	protected String reducedSPName = null;

	protected static BigDecimal stabilityThresh = BigDecimal.ONE;

	public static final double RT = 1.9891/1000.0 * 298.15;
	public static double targetEpsilon = 0.03;
	protected double effectiveEpsilon = 1.0;
	protected static String phase2Method = "fast";

	protected int strand = -1;
	protected KSConfigFileParser cfp = null;
	protected KSSearchProblem reducedSP = null;
	private boolean isContinuous = true;
	protected ArrayList<Integer> absolutePos = null;
	protected KSSearchProblem panSP = null;
	private boolean isFullyDefined = true;

	protected BigDecimal qStar = BigDecimal.ZERO;
	protected BigDecimal qPrime = BigDecimal.ZERO;
	protected BigDecimal pStar = BigDecimal.ZERO;

	protected ExpFunction e = new ExpFunction();
	protected double Et = 0;
	protected double E0 = 0;

	protected BigInteger prunedConfs = BigInteger.ZERO;
	protected BigInteger unPrunedConfs = BigInteger.ZERO;
	protected HashSet<ArrayList<Integer>> processedConfsSet = new HashSet<>();
	protected BigDecimal partialQLB = BigDecimal.ZERO;
	protected BigInteger processedConfsDuringInterval = BigInteger.ZERO;
	protected BigInteger processedConfs = BigInteger.ZERO;
	protected BigInteger processingConfs = BigInteger.ZERO; // # confs being minimized at this instant

	protected PriorityQueue<KSConf> topConfsPQ = null;

	protected PFAbstract() {}

	protected PFAbstract( int strand, ArrayList<String> sequence, 
			ArrayList<Integer> absolutePos, 
			String checkPointPath, String reducedSPName, 
			KSConfigFileParser cfp, KSSearchProblem panSP ) {

		this.sequence = sequence; this.sequence.trimToSize();
		this.absolutePos = absolutePos; this.absolutePos.trimToSize();
		this.checkPointPath = checkPointPath;
		this.reducedSPName = reducedSPName;
		this.panSP = panSP;
		this.strand = strand;
		this.isFullyDefined = sequence.size() == panSP.confSpace.numPos ? true : false;
		this.isContinuous = panSP.contSCFlex;
		this.reducedSP = createReducedSP(panSP.contSCFlex, strand, sequence, absolutePos);
		this.cfp = cfp;

		Comparator<KSConf> comparator = new KSConf(new ArrayList<>(), 0.0).new KSConfMinEComparator();
		topConfsPQ = new PriorityQueue<KSConf>(getNumTopConfsToSave(), comparator);
	}


	public static void setPhase2Method( String in ) {
		ArrayList<String> allowedMethods = new ArrayList<String>(Arrays.asList("fast", "slow"));
		in = in.toLowerCase();

		if(!allowedMethods.contains(in) )
			throw new RuntimeException("ERROR: allowed values of parameter kStarPhase2Method are "
					+ "'fast' and 'slow'");

		phase2Method = in;
	}


	public static String getPhase2Method() {
		return phase2Method;
	}

	
	public PruningMatrix getReducedPruningMatrix() {
		return reducedSP.reducedMat;
	}
	
	
	public PruningMatrix getPanPruningMatrix() {
		return panSP.pruneMat;
	}
	

	public PruningMatrix getReducedInversePruningMatrix() {
		return reducedSP.inverseMat;
	}
	
	
	public PruningMatrix getPanInversePruningMatrix() {
		return panSP.inverseMat;
	}


	public ArrayList<Integer> getAbsolutePos() {
		return absolutePos;
	}


	public boolean isContinuous() {
		return isContinuous;
	}


	public String getFlexibility() {
		return isContinuous ? "(continuous)": "(discrete)";
	}


	protected BigDecimal reComputePartialQLB( ConfSearch confSearch ) {
		partialQLB = BigDecimal.ZERO;

		for(ArrayList<Integer> conf : processedConfsSet) {
			int[] confArray = KSConf.list2Array(conf);
			partialQLB = partialQLB.add( getBoltzmannWeight(getConfBound(confSearch, confArray)) );
		}

		return partialQLB;
	}


	protected boolean canUseHotByManualSelection() {
		if(!isContinuous()) return false;

		if(!getHotMethod().equalsIgnoreCase("manual")) return false;

		if(!isFullyDefined()) return false;

		if(!reducedSP.contSCFlex) return false;

		return true;
	}


	protected boolean canUseHotByConfError( double boundError ) {
		if(!canUseHotByConfError()) return false;

		if(boundError < getHotBoundPct()) return false;

		return true;
	}


	protected boolean canUseHotByConfError() {
		if(!getHotMethod().equalsIgnoreCase("error")) return false;

		if(reducedSP.confSpace.numPos - getNumPosInHOTs() < getHotNumRes()) return false;

		if(!isFullyDefined()) return false;

		if(!reducedSP.contSCFlex) return false;

		return true;
	}


	public boolean HOTsContains(Integer pos) {
		for(ArrayList<Integer> hot : HOTs)
			if(hot.contains(pos)) return true;

		return false;
	}


	protected int getMaxHOTSize( boolean max ) {
		int ans = max == true ? Integer.MIN_VALUE : Integer.MAX_VALUE;

		for(ArrayList<Integer> hot : HOTs) {
			if(max)
				ans = (int)Math.max(ans, hot.size());

			else
				ans = (int)Math.min(ans, hot.size());
		}

		return ans;
	}


	protected int getNumPosInHOTs() {
		int ans = 0;

		for(ArrayList<Integer> hot : HOTs)
			ans += hot.size();

		return ans;
	}


	protected void memoizePosInHot(int... pos) {
		ArrayList<Integer> hot = new ArrayList<>();
		for(int i : pos) hot.add(i);
		HOTs.add(hot);
	}


	public boolean isFullyDefined() {
		return isFullyDefined;
	}


	protected BigDecimal productUndefinedRots() {

		BigDecimal ans = BigDecimal.ONE;

		// get unassigned residue positions
		int numPos = panSP.confSpace.numPos;
		ArrayList<Integer> undefinedPos = new ArrayList<>(numPos);
		for(int pos = 0; pos < numPos; ++pos) undefinedPos.add(pos);
		undefinedPos.removeAll(absolutePos);

		boolean minimizeProduct = isContinuous() ? false : true;
		
		for( int pos : undefinedPos ) {
			
			long unprunedSize = minimizeProduct ? Integer.MAX_VALUE : Integer.MIN_VALUE;
			long prunedSize = minimizeProduct ? Integer.MAX_VALUE : Integer.MIN_VALUE;
			
			ArrayList<String> aasAtPos = panSP.getAAsAtPos(panSP.pruneMat, pos);
			for(String aa : aasAtPos) {
				ArrayList<Integer> aaUnprunedRCNums = panSP.rcsAtPos(panSP.pruneMat, pos, aa, false);
				unprunedSize = minimizeProduct ? Math.min(unprunedSize, aaUnprunedRCNums.size()) : Math.max(unprunedSize, aaUnprunedRCNums.size());
				
				if(!minimizeProduct) {
					ArrayList<Integer> aaPrunedRCNums = panSP.rcsAtPos(panSP.pruneMat, pos, aa, true);
					prunedSize = Math.max(prunedSize, aaPrunedRCNums.size());
				}
			}
			
			if(minimizeProduct) prunedSize = 0;
			ans = ans.multiply(BigDecimal.valueOf(unprunedSize + prunedSize));
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


	public ConfSearch getConfTree( boolean invertPruneMat ) {

		PruningMatrix reducedPmat = invertPruneMat ? reducedSP.inverseMat : reducedSP.reducedMat;
		
		if( isFullyDefined() ) {
			if(reducedPmat != null)
				return new ConfTree<FullAStarNode>(new FullAStarNode.Factory(reducedSP.posNums.size()), reducedSP, reducedPmat);

			else 
				return null;
		}

		else {
			
			PruningMatrix panPmat = invertPruneMat ? panSP.inverseMat : panSP.pruneMat;
			
			if(reducedPmat != null && panPmat != null)
				return new KAStarConfTree(reducedSP, reducedPmat, panPmat);

			else 
				return null;
		}
	}


	public double getConfBound( ConfSearch confSearch, int[] conf ) {
		double bound = 0;

		if( isFullyDefined() ) {
			bound = reducedSP.lowerBound(conf);
		}

		else
			bound = ((KAStarConfTree)confSearch).confEnergyBound(conf);

		return bound;
	}


	public String getReducedSearchProblemName() {
		return reducedSPName;
	}


	public int getStrand() {
		return strand;
	}


	public void setPanSeqSP( KSSearchProblem other ) {
		this.panSP = other;
	}


	public HashSet<ArrayList<Integer>> getProcessedConfsSet() {
		return processedConfsSet;
	}


	public BigInteger getNumProcessed() {
		return processedConfs;
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
			reducedSP.outputMinimizedStruct(tmp.poll().getConfArray(), pdbName);
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
		if( eAppx != EApproxReached.FALSE ) return getQStar();

		updateQPrime();
		return (qStar.add(qPrime)).add(pStar);
	}


	public void setNumUnPruned() {
		unPrunedConfs = reducedSP.numConfs(getReducedPruningMatrix());
	}


	public void setNumPruned() {
		prunedConfs = reducedSP.numConfs(getReducedInversePruningMatrix());
	}


	protected double computeEffectiveEpsilon() {
		
		BigDecimal dividend = qPrime.add(pStar);
		BigDecimal divisor = qStar.add(dividend);
		
		// energies are too high so epsilon can never be reached
		if( divisor.compareTo(BigDecimal.ZERO) == 0 ) return EPSILON_NEVER_POSSIBLE;

		/*
		BigDecimal maxSubtrahend = qPrime.max(pStar);
		BigDecimal minEpsilon = (dividend.subtract(maxSubtrahend)).divide(divisor.subtract(maxSubtrahend), 4);
		if(minEpsilon.doubleValue() > targetEpsilon) {
			System.out.println("WARNING: minEpsilon = " + minEpsilon + " > targetEpsilon = " + targetEpsilon);
			return EPSILON_PHASE_2;
		}
		*/
		
		else return dividend.divide(divisor, 4).doubleValue();
	}


	public double getEffectiveEpsilon() {
		return effectiveEpsilon;
	}


	protected void initTradPStar() {

		ConfSearch confSearch = getConfTree(true);

		if(confSearch == null) {
			setPStar( Double.POSITIVE_INFINITY );
			return;
		}

		ScoredConf conf;
		if( (conf = confSearch.nextConf()) != null ) {
			setPStar( getConfBound(confSearch, conf.getAssignments()) );
		}

		else 
			setPStar( Double.POSITIVE_INFINITY );
	}


	private void setPStar( double eLB ) {
		E0 = eLB;
		pStar = ( getBoltzmannWeight( E0 )).multiply( new BigDecimal(prunedConfs) );
	}


	protected void updateQPrime() {
		qPrime = getBoltzmannWeight( Et ).multiply(new BigDecimal(getNumUnEnumerated()));
	}


	protected void updateQStar( KSConf conf ) {

		if( processedConfsSet.contains(conf.getConf()) ) 
			return;
		
		processedConfsSet.add(conf.getConf());
		
		qStar = qStar.add( getBoltzmannWeight( conf.getEnergy() ) );

		if( getImpl().toLowerCase().contains("parallel") )
			partialQLB = partialQLB.add( getBoltzmannWeight( conf.getEnergyBound() ) );
		
		if(isFullyDefined() && saveTopConfsAsPDB)
			saveTopConf(conf);

		processedConfs = processedConfs.add(BigInteger.ONE);
		processedConfsDuringInterval = processedConfsDuringInterval.add(BigInteger.ONE);
	}


	public BigDecimal getBoltzmannWeight( double E ) {
		return e.exp(-E / RT);
	}

	protected void exitIfTimeOut() {
		if( KSAbstract.runTimeout != 0 && System.currentTimeMillis()-startTime > KSAbstract.runTimeout * 86400000 )
			throw new RuntimeException("ERROR: running time exceeds " + KSAbstract.runTimeout + " days!");
	}
	
	public abstract void start();

	public void runSlice(long target) {

		while( eAppx == EApproxReached.FALSE && getProcessedDuringInterval().longValue() < target ) {
			computeSlice();

			if( !doingPhase2 && eAppx == EApproxReached.NOT_POSSIBLE ) {
				phase2();

				//if( eAppx == EApproxReached.FALSE )
				//	computeSlice();
			}
			
			else if( doingPhase2 && eAppx == EApproxReached.NOT_POSSIBLE )
				System.out.println("\nCan never reach target epsilon approximation of " + targetEpsilon + " for sequence: " + KSAbstract.list1D2String(sequence, " ") + " " + getFlexibility());
		}

		if( isFullyDefined() && saveTopConfsAsPDB && eAppx == EApproxReached.TRUE ) writeTopConfs();

		if( eAppx != EApproxReached.FALSE ) 
			cleanup();

		resetProcessedDuringInterval();
	}


	public void runToCompletion() {
		compute();

		if( eAppx == EApproxReached.NOT_POSSIBLE ) {
			phase2();

			if( eAppx == EApproxReached.FALSE ) {
				compute();
				
				if( eAppx != EApproxReached.TRUE )
					System.out.println("\nCan never reach target epsilon approximation of " + targetEpsilon + " for sequence: " + KSAbstract.list1D2String(sequence, " ") + " " + getFlexibility());
			}
		}

		if( isFullyDefined() && saveTopConfsAsPDB ) writeTopConfs();

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
		processingConfs = BigInteger.ZERO;
		resetProcessedDuringInterval();
		start();
	}

	protected void phase2() {

		System.out.println("\nCould not reach target epsilon approximation of " + targetEpsilon + " for sequence: " + KSAbstract.list1D2String(sequence, " ") + " " + getFlexibility());
                
		doingPhase2 = true;
		
		if( getEffectiveEpsilon() == EPSILON_NEVER_POSSIBLE ) {
			// we can never reach epsilon because q* + q' + p* = 0
			System.out.println("\nCan never reach target epsilon approximation of " + targetEpsilon + " for sequence: " + KSAbstract.list1D2String(sequence, " ") + " " + getFlexibility());
			return;
		}

		System.out.println("Attempting Phase 2...\n");

		if(HOTs != null && HOTs.size() > 0) {
			// considerations for HOT
			HOTs.clear(); // clear HOTs
		}

		// completely relax pruning
		double maxPruningInterval = cfp.params.getDouble("StericThresh");
		rePruneReducedSP(maxPruningInterval);
		
		setNumUnPruned();
		setNumPruned(); // needed for p*

		if( !phase2Method.equalsIgnoreCase("fast") || !isFullyDefined() ) {
			// conservative implementation that re-enumerates all
			restart();
			return;
		}

		else {
			// shortcut implementation that follows the protocol described in the paper:
			// get new value of epsilon with q' = 0 and the reduced (unpruned) p*

			qPrime = BigDecimal.ZERO;

			// MUST call abstract base class version of this method
			initTradPStar(); // re-calculate p*

			double effectiveEpsilon = computeEffectiveEpsilon();

			if( getQStar().compareTo(BigDecimal.ZERO) > 0 
					&& effectiveEpsilon != EPSILON_NEVER_POSSIBLE && effectiveEpsilon <= targetEpsilon ) {

				setEpsilonStatus(EApproxReached.TRUE);

				System.out.println("\nReached target epsilon approximation of " + targetEpsilon + " for sequence: " +
						KSAbstract.list1D2String(sequence, " ") + " " + getFlexibility());
				return;
			}

			else {
				System.out.println("\nRe-starting K* for sequence: " + KSAbstract.list1D2String(sequence, " ") + " " + getFlexibility());
				restart();
			}
		}
	}


	public void cleanup() {
		e = null;
		HOTs = null;
		cfp = null;
		panSP = null;
		reducedSP = null;
		processedConfsSet = null;
	}


	public KSSearchProblem getReducedSearchProblem() {
		return reducedSP;
	}


	public KSSearchProblem getPanSeqSearchProblem() {
		return panSP;
	}


	public BigInteger getNumPruned() {
		return prunedConfs;
	}


	protected BigInteger getNumUnEnumerated() {
		return unPrunedConfs.subtract(getNumProcessed());
	}


	public BigInteger getNumUnPruned() {
		return unPrunedConfs;
	}


	protected BigInteger getProcessedDuringInterval() {
		return processedConfsDuringInterval;
	}


	protected void resetProcessedDuringInterval() {
		processedConfsDuringInterval = BigInteger.ZERO;
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


	public static int getConfsThreadBuffer() {
		return threadConfsBuffer;
	}


	public static void setConfsThreadBuffer( int confsBuffer ) {
		if( confsBuffer > 0 ) threadConfsBuffer = confsBuffer;
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


	public static void setHotBoundPct( String param, double value ) {
		if( value < 0.0 || value > 1.0 )
			throw new RuntimeException("ERROR: value of " + param + " must be in the range [0.0, 1.0]");

		hotBoundPct = value;
	}


	public static double getHotBoundPct() {
		return hotBoundPct;
	}


	public static void setHotTopRotsPct( String param, double value ) {
		if( value < 0.0 || value > 1.0 )
			throw new RuntimeException("ERROR: " + param + " must be in the range [0.0, 1.0]");

		hotTopRotsPct = value;
	}


	public static double getHotTopRotsPct() {
		return hotTopRotsPct;
	}


	public static void setHotNumRes( String param, int value ) {
		if( value < 3 || value > 4 ) 
			throw new RuntimeException("ERROR: allowed values of " + param + " are [3, 4]");

		hotNumRes = value;
	}


	public static int getHotNumRes() {
		return hotNumRes;
	}


	public static void setHotMethod( String param, String value ) {

		switch( value ) {

		case "none":
		case "error":
		case "manual":
			hotMethod = value;
			break;

		default:
			throw new RuntimeException("ERROR: allowed values of " + param + " are {none|error|manual}");
		}
	}


	public static String getHotMethod() {
		return hotMethod;
	}


	public static void setMaxKSconfs( long in ) {
		if( in < 1 ) in = 1;
		maxKSConfs = in;
	}


	public boolean maxKSConfsReached() {
		return useMaxKSConfs && getNumProcessed().longValue() >= maxKSConfs;
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

		case "traditional":
		case "ub":
		case "parallel0":
		case "parallel1":
		case "parallel2":
		case "parallelconf":
			pFuncCFGImpl = implementation;
			break;

		default:
			throw new RuntimeException("ERROR: specified value of implementation " + implementation.toLowerCase() + " is invalid");
		}
	}


	protected KSSearchProblem createReducedSP( boolean contSCFlex, int strand, 
			ArrayList<String> seq, ArrayList<Integer> absolutePos ) {

		KSSearchProblem reducedSP = panSP.getReducedSearchProblem(reducedSPName, 
				KSAbstract.list1D2ListOfLists(KSAllowedSeqs.getAAsFromSeq(seq)), 
				KSAllowedSeqs.getFlexResFromSeq(seq), 
				absolutePos);

		return reducedSP;
	}
	
	
	public void rePruneReducedSP(double pruningInterval) {	
		panSP = (KSSearchProblem) ObjectIO.deepCopy(panSP);
		if(panSP.emat == null) panSP.emat = (EnergyMatrix) ObjectIO.readObject(panSP.getMatrixFileName(panSP.getMatrixType()), false);
		cfp.setupPruning(panSP, pruningInterval, panSP.useEPIC, panSP.useTupExpForSearch).prune();
		reducedSP = createReducedSP(panSP.contSCFlex, strand, sequence, absolutePos);
	}

}