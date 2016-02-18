package edu.duke.cs.osprey.kstar;

import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.math.BigInteger;
import java.net.InetAddress;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.concurrent.ConcurrentHashMap;

import edu.duke.cs.osprey.confspace.AllowedSeqs;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.kstar.KSCalc.KSCalcType;
import edu.duke.cs.osprey.kstar.PFAbstract.EApproxReached;
import edu.duke.cs.osprey.pruning.PruningControl;
import edu.duke.cs.osprey.tools.ObjectIO;

/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 *
 */
public class KSImplementationLinear {

	HashMap<Integer, SearchProblem> strand2AllSearchProblem;
	HashMap<Integer, AllowedSeqs> strand2AllAllowedSequence;
	HashMap<Integer, PruningControl> strand2AllPruning;

	HashMap<ArrayList<String>, Integer> seq2SeqID = new HashMap<>();
	HashMap<Integer, ArrayList<String>> seqID2Seq = new HashMap<>();
	Map<Integer, KSCalc> calculations = new ConcurrentHashMap<Integer, KSCalc>();
	HashSet<Integer> observedSeqs = new HashSet<>();
	ArrayList< ArrayList<String> > prunedSeqs = new ArrayList<>();
	HashMap<Integer, HashMap<ArrayList<String>, PFAbstract>> strand2PrecomputedPFs =
			new HashMap<>();

	static String origRunName;
	String modifiedRunName;
	static boolean parallelExpansion = false;
	static String path = null;
	ConfigFileParser cfp;
	HashMap<KSCalcType, BigInteger> numMinimizedConfs = new HashMap<>();
	BigInteger observedMaxSearchSize = BigInteger.ZERO;
	KSCalc wtSeq = null;
	boolean wtObserved = false;

	long start;


	public KSImplementationLinear(ArrayList<ArrayList<String>> mutations,
			HashMap<Integer, SearchProblem> strand2AllSearchProblems,
			HashMap<Integer, AllowedSeqs> strand2AllAllowedSequences,
			HashMap<Integer, PruningControl> strand2AllPruning,
			ConfigFileParser cfp) {

		try {

			this.strand2AllSearchProblem = strand2AllSearchProblems;
			this.strand2AllAllowedSequence = strand2AllAllowedSequences;
			this.strand2AllPruning = strand2AllPruning;

			this.cfp = cfp;

			numMinimizedConfs.put(KSCalcType.FULL_KSTAR, BigInteger.ZERO);
			numMinimizedConfs.put(KSCalcType.PARTIAL_SUBSEQUENCE, BigInteger.ZERO);

			origRunName = cfp.getParams().getValue("runName");

			modifiedRunName = origRunName
					+ "." + InetAddress.getLocalHost().getHostName()
					+ "." + cfp.getParams().getValue("KStarMethod", "linear")
					+ "." + cfp.getParams().getValue("pFuncMethod", "1npcpmcache")
					+ ".txt";

			String expansionType = cfp.getParams().getValue("KStarExpansion", "serial");
			parallelExpansion = expansionType.equalsIgnoreCase("serial") ? false : true;

			deleteExistingOutputFile();

			path = cfp.getParams().getValue("ematdir", "emat.linear") + File.separator + "subSeq";
			ObjectIO.makeDir(path, false);

			mapSeqAndSeqID(Strand.COMPLEX);
			//printSeq2SeqID();

		} catch (Exception e) {
			System.out.println(e.getMessage());
			e.printStackTrace();
			System.exit(1);
		}
	}


	protected static String getSearchProblemName( int strand, ArrayList<String> seq ) {
		// terrible use of static here, but i wanted the code to appear once
		return path + File.separator + origRunName +"."+ Strand.getStrandString(strand) + "." + KSCalc.getSequenceSignature(seq);
	}


	protected void mapSeqAndSeqID(int strand) {

		int numSeqs = strand2AllAllowedSequence.get(strand).getStrandSeqList().size();

		for( int seq = 0; seq < numSeqs; ++seq ) {
			seq2SeqID.put(strand2AllAllowedSequence.get(strand).getStrandSeqList().get(seq), seq);
			seqID2Seq.put(seq, strand2AllAllowedSequence.get(strand).getStrandSeqList().get(seq));
		}
	}


	protected KSCalc createKSCalc(ArrayList<String> seq, int id) {

		System.out.print("\nInitializing fully defined sequence ");
		KSCalc.print(seq, System.out);
		System.out.println();

		KSCalc mutSeq = new KSCalc( id, strand2AllAllowedSequence, 
				strand2AllSearchProblem, strand2AllPruning, cfp );

		if(!mutSeq.isValid()) {

			System.out.print("\nThe sequence "); 
			mutSeq.printSequence(Strand.COMPLEX, System.out);
			System.out.println(" has zero RCs at one or more states and will be pruned");

			seq2SeqID.remove(seq);
			seqID2Seq.remove(id);

			mutSeq = null;
		}

		return mutSeq;
	}


	public void printSequences() {

		Set<Integer> sortedKeySet = new TreeSet<>(seqID2Seq.keySet());
		for( int seqID : sortedKeySet ) {

			System.out.print("Sequence " + seqID + ": ");
			KSCalc.print(seqID2Seq.get(seqID), System.out);
			System.out.println();
		}
	}


	protected void insertCalculation( int seqID, KSCalc calc ) {
		calculations.put(seqID, calc);
	}


	protected KSCalc removeCalculation( int seqID ) {
		KSCalc val = calculations.remove(seqID);

		if(val == null) System.out.println("WARNING: key " + seqID + " has no assigned value");

		return val;
	}


	protected TreeSet<Integer> getSequenceIDs() {
		return new TreeSet<>(seqID2Seq.keySet());
	}


	public int numSeqs() {
		return seqID2Seq.size();
	}


	public BigInteger getUnobservedMaxSearchSize() {

		System.out.println("\nCalculating total search space size");

		BigInteger ans = BigInteger.ZERO;
		ArrayList<KSCalc> tmpCalcs = new ArrayList<>();
		ArrayList<Integer> tmpIDs = new ArrayList<>();

		int capacity = Runtime.getRuntime().availableProcessors()-1;
		ArrayList<Integer> indexes = new ArrayList<>();

		for(int i = 0; i < capacity; ++i) {
			indexes.add(i);
			tmpCalcs.add(null);
		}
		indexes.trimToSize();
		tmpCalcs.trimToSize();

		ArrayList<Integer> unobservedSeqIDs =  new ArrayList<>();

		// compile list of unovserved ids
		for( int seqID : getSequenceIDs() ) {
			if( !observedSeqs.contains(seqID) )
				unobservedSeqIDs.add(seqID);
		}

		int total = unobservedSeqIDs.size();
		int count = 0;

		if( total == 0 ) return BigInteger.ZERO;

		while( true ) {

			// fill the unobserved calcs to capacity
			while( unobservedSeqIDs.size() > 0 ) {
				tmpIDs.add(unobservedSeqIDs.remove(0));
				if( tmpIDs.size() == capacity || unobservedSeqIDs.size() == 0 ) 
					break;
			}

			if( tmpIDs.size() == 0 ) break;

			while( indexes.size() > tmpIDs.size() ) {
				indexes.remove(indexes.size()-1);
				tmpCalcs.remove(tmpCalcs.size()-1);
			}

			count += tmpIDs.size();

			System.out.println("\nCalculating search space size for unobserved sequence # " + count + " / " + total);

			if( parallelExpansion ) {
				indexes.parallelStream().forEach( i -> {
					int seqID = tmpIDs.get(i);
					tmpCalcs.set( i, createKSCalc(seqID2Seq.get(seqID), seqID) );
				});
			}

			else {
				for( Integer i : indexes ) {
					int seqID = tmpIDs.get(i);
					tmpCalcs.set( i, createKSCalc(seqID2Seq.get(seqID), seqID) );
				}
			}

			for( KSCalc mutSeq : tmpCalcs ) {
				if( mutSeq != null )
					ans = ans.add(mutSeq.getNumInitialUnPrunedConfs());
			}

			tmpIDs.clear();
		}

		return ans;
	}


	protected void updateNumMinimizedConfs( KSCalcType type, BigInteger val ) {

		if( val.compareTo(BigInteger.ZERO) == 0 ) return;

		if( type != KSCalcType.FULL_KSTAR )
			type = KSCalcType.PARTIAL_SUBSEQUENCE;

		BigInteger ans = numMinimizedConfs.get(type);
		numMinimizedConfs.put(type, ans.add(val));
	}


	protected BigInteger getTotalMinimizedConfs() {
		BigInteger ans = BigInteger.ZERO;
		for( BigInteger value : numMinimizedConfs.values() ) {
			ans = ans.add(value);
		}

		return ans;
	}


	protected boolean isWT( KSCalc mutSeq ) {
		return mutSeq.complex.getSequence().equals(cfp.getWTSequence());
	}


	protected PFAbstract getPrecomputedPF( int strand, ArrayList<String> seq ) {

		HashMap<ArrayList<String>, PFAbstract> seq2PF = strand2PrecomputedPFs.get(strand);

		if(seq2PF == null)
			return null;

		PFAbstract pf = seq2PF.get(seq);

		if(pf == null)
			return null;

		System.out.print("Getting pre-computed " + Strand.getStrandString(strand) + " partition function: ");
		KSCalc.print(pf.getSequence(), System.out);
		System.out.println();

		return pf;
	}


	protected void setPrecomputedPF( int strand, PFAbstract pf ) {

		int occurrences = Collections.frequency(strand2AllAllowedSequence.get(strand).getStrandSeqList(), pf.getSequence());
		if( occurrences <= 1 )
			return;

		HashMap<ArrayList<String>, PFAbstract> seq2PF = null;

		if( !strand2PrecomputedPFs.containsKey(strand) )
			strand2PrecomputedPFs.put(strand, new HashMap<ArrayList<String>, PFAbstract>());

		seq2PF = strand2PrecomputedPFs.get(strand);

		seq2PF.put(pf.getSequence(), pf);
	}


	public void runWTSeq() {
		// prevent wt from entering subsequence tree
		ArrayList<String> wt = cfp.getWTSequence();
		int wtSeqID = seq2SeqID.remove(wt);
		seqID2Seq.remove(wtSeqID);

		// create and run wt
		System.out.println("\nCreating wild type sequence object");

		double oldPFInterval = PFAbstract.getInterval();

		PFAbstract.setInterval("max");

		wtSeq = createKSCalc(wt, wtSeqID);
		if( !wtSeq.isValid() ) {
			throw new RuntimeException("\nERROR: The wild type sequence "
					+ "has zero RCs at one or more positions and will be pruned");
		}

		wtSeq.expand(wtSeq);

		if( wtSeq.getEpsilon() == EApproxReached.NOT_POSSIBLE ) {

			System.out.print("K* aborting sequence ");
			wtSeq.printSequence(Strand.COMPLEX, System.out);
			System.out.println(" : cannot reach epsilon");

			throw new RuntimeException("ERROR: WT sequence cannot reach epsilon. "
					+ "Relax pruning criteria and try again.");
		}

		setPrecomputedPF(Strand.LIGAND, wtSeq.getPartitionFunction(Strand.LIGAND));
		setPrecomputedPF(Strand.PROTEIN, wtSeq.getPartitionFunction(Strand.PROTEIN));

		writeCalcSummary(wtSeq);

		PFAbstract.setInterval(String.valueOf(oldPFInterval));
	}


	public void run() {

		System.out.println("\nCalculating K* for the following " + numSeqs() + " sequence(s)");

		printSequences();

		KSCalc mutSeq = null;
		int total = numSeqs();

		start = System.currentTimeMillis();

		System.out.println("\nCalculating K* for sequence # " + 0 + " / " + (total-1) 
				+ " at " + (System.currentTimeMillis()-start)/1000 + " seconds"
				+ " after evaluating " + numMinimizedConfs.get(KSCalcType.FULL_KSTAR) + " conformations");
		System.out.println();

		runWTSeq();

		updateSearchSpaceStatistics(wtSeq);

		// run all K* calcuations
		double oldpFuncInterval = PFAbstract.getInterval();
		PFAbstract.setInterval("max");

		for( int seqID : getSequenceIDs() ) {

			System.out.println("\nCalculating K* for sequence # " + seqID + " / " + (total-1) 
					+ " at " + (System.currentTimeMillis()-start)/1000 + " seconds"
					+ " after evaluating " + numMinimizedConfs.get(KSCalcType.FULL_KSTAR) + " conformations");
			System.out.println();

			if( (mutSeq = createKSCalc(seqID2Seq.get(seqID), seqID)) == null )  {

				continue;
			}

			// get pf if it already exists
			PFAbstract ligandPF = getPrecomputedPF(Strand.LIGAND, mutSeq.getSequence(Strand.LIGAND));
			PFAbstract proteinPF = getPrecomputedPF(Strand.PROTEIN, mutSeq.getSequence(Strand.PROTEIN));
			if( ligandPF != null )
				mutSeq.setPartitionFunction(Strand.LIGAND, ligandPF);
			if( proteinPF != null)
				mutSeq.setPartitionFunction(Strand.PROTEIN, proteinPF);

			if( mutSeq.getEpsilon() == EApproxReached.FALSE ) {
				
				mutSeq.run(wtSeq);

				// store protein and ligand pf so it can be re-used later
				setPrecomputedPF(Strand.LIGAND, mutSeq.getPartitionFunction(Strand.LIGAND));
				setPrecomputedPF(Strand.PROTEIN, mutSeq.getPartitionFunction(Strand.PROTEIN));

				updateSearchSpaceStatistics(mutSeq);
			}

			if( mutSeq.getEpsilon() != EApproxReached.TRUE ) {

				System.out.print("K* aborting sequence ");
				mutSeq.printSequence(Strand.COMPLEX, System.out);

				switch(mutSeq.getEpsilon()) {
				case NOT_POSSIBLE:
					System.out.println(" : protein or ligand cannot reach epsilon");
					break;

				case NOT_STABLE:
					System.out.println(" : unbound protein or ligand is not stable");
					break;

				default:
					break;
				}

				continue;
			}

			writeCalcSummary(mutSeq);
		}

		PFAbstract.setInterval(String.valueOf(oldpFuncInterval));

		if( prunedSeqs.size() == total ) {
			throw new RuntimeException("\nERROR: All K* sequences were pruned. "
					+ "Relax pruning criteria and try again.");
		}

		long end = System.currentTimeMillis();

		BigInteger totalSearchSize = observedMaxSearchSize.add(getUnobservedMaxSearchSize());

		// total+1 for wt
		System.out.println("\nK* calculations computed: " + observedSeqs.size());
		System.out.println("\nK* conformations minimized: " + numMinimizedConfs.get(KSCalcType.FULL_KSTAR));
		System.out.println("\nTotal # of conformations in search space: " + totalSearchSize);
		System.out.println("\nK* running time: " + (end-start)/1000 + " seconds");
	}


	protected void updateSearchSpaceStatistics( KSCalc calc ) {

		EApproxReached eAppx = calc.getEpsilon();

		if( eAppx == EApproxReached.NOT_POSSIBLE ) {
			// don't count anything
			int id = calc.getSeqID();
			ArrayList<String> seq = calc.getSequence(Strand.COMPLEX);

			observedSeqs.remove(id);
			prunedSeqs.add(seq);

			seqID2Seq.remove(id);
			seq2SeqID.remove(seq);
		}

		else {
			observedSeqs.add(calc.getSeqID());
			updateNumMinimizedConfs(calc.getType(), calc.getNumMinimizedConfsDuringInterval());
			observedMaxSearchSize = observedMaxSearchSize.add(calc.getNumInitialUnPrunedConfs());

			if( !wtObserved && isWT(calc) ) {
				updateNumMinimizedConfs(calc.getType(), calc.ligand.getNumMinimizedConfsDuringInterval());
				calc.ligand.resetNumMinimizedConfsDuringInterval();
				observedMaxSearchSize = observedMaxSearchSize.add(calc.ligand.getNumInitialUnPrunedConfs());
				wtObserved = true;
			}
		}
	}


	public void printSummary() {

		System.out.println("\nSymmary");
		for( int seqID : getSequenceIDs() ) {
			removeCalculation(seqID).printSummary( System.out );
		}
	}


	protected void deleteExistingOutputFile() {
		try {

			File f = new File(modifiedRunName);
			if( f.exists() ) f.delete();

		} catch (Exception e) {
			System.out.println(e.getMessage());
			e.printStackTrace();
			System.exit(1);
		}
	}


	protected void writeCalcSummary(KSCalc calc) {
		try( PrintStream out = new PrintStream(new FileOutputStream(modifiedRunName, true)) ) {

			calc.printSummary(out);

		} catch (Exception e) {
			System.out.println(e.getMessage());
			e.printStackTrace();
			System.exit(1);
		}
	}

}
