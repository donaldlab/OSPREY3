package edu.duke.cs.osprey.kstar.impl;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ForkJoinPool;
import java.util.stream.IntStream;

import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.kstar.AllowedSeqs;
import edu.duke.cs.osprey.kstar.KSAbstract;
import edu.duke.cs.osprey.kstar.KSCalc;
import edu.duke.cs.osprey.kstar.Strand;
import edu.duke.cs.osprey.kstar.pfunc.PFAbstract;
import edu.duke.cs.osprey.kstar.pfunc.PFAbstract.EApproxReached;
import edu.duke.cs.osprey.parallelism.ThreadParallelism;

/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 *
 */
public class KSImplLinear extends KSAbstract {

	public KSImplLinear( ConfigFileParser cfp ) {
		super( cfp );
	}

	public void init( HashMap<Integer, AllowedSeqs> strand2AllowedSeqs ) {

		this.strand2AllowedSeqs = strand2AllowedSeqs;

		printSequences();

		createEmatDir();

		ArrayList<Boolean> contSCFlexVals = new ArrayList<Boolean>(Arrays.asList(true, false));
		createEmats(contSCFlexVals);
	}


	protected void prepareAllSingleSeqSPs( ArrayList<Boolean> contSCFlexVals ) {

		try {

			for( boolean contSCFlex : contSCFlexVals ) {

				ForkJoinPool forkJoinPool = new ForkJoinPool(ThreadParallelism.getNumThreads());
				forkJoinPool.submit(() ->

				IntStream.range(0, strand2AllowedSeqs.get(Strand.COMPLEX).getNumSeqs()).parallel().forEach(i -> {

					//for( int i = 0; i < strand2AllowedSeqs.get(Strand.COMPLEX).getNumSeqs(); ++i ) {

					System.out.println("\nCreating search problem for sequence " + 
							i + "/" + (strand2AllowedSeqs.get(Strand.COMPLEX).getNumSeqs()-1) + "\n");

					ConcurrentHashMap<Integer, SearchProblem> map = createSPsForSeq(i, contSCFlex);

					// put partition function in list, so we can parallelize energy matrix computation
					for(int strand : map.keySet()) {	
						SearchProblem sp = map.get(strand);
						name2SP.put(sp.name, sp);
					}

					//}
				})).get();
			}
			
			loadAndPruneMatrices(); 
			prunedSingleSeqs = true;

		} catch (Exception e) {
			System.out.println(e.getMessage());
			e.printStackTrace();
			System.exit(1);
		}
	}


	protected ConcurrentHashMap<Integer, SearchProblem> createSPsForSeq(int i, boolean contSCFlex) {
		// used to precompute energy matrices
		ConcurrentHashMap<Integer, SearchProblem> ans = new ConcurrentHashMap<>();

		ArrayList<Integer> strands = new ArrayList<>();
		strands.add(Strand.COMPLEX); strands.add(Strand.PROTEIN); strands.add(Strand.LIGAND);

		strands.parallelStream().forEach((strand) -> {

			AllowedSeqs strandSeqs = strand2AllowedSeqs.get(strand);

			ArrayList<String> seq = strandSeqs.getStrandSeqAtPos(i);
			
			ArrayList<Integer> flexResIndexes = strandSeqs.getFlexResIndexesFromSeq(seq);

			String spName = getSearchProblemName(contSCFlex, strand, seq);
			
			SearchProblem seqSP = null;			
			if( (seqSP = name2SP.get(spName)) == null ) {
				seqSP = createSingleSeqSPFast( contSCFlex, strand, seq, flexResIndexes );
			}

			// synchronized
			ans.put(strand, seqSP);
		});

		return ans;
	}


	@Override
	public String getKSMethod() {
		return "linear";
	}

	@Override
	public void run() {

		// each value corresponds to the desired flexibility of the 
		// pl, p, and l conformation spaces, respectively
		ArrayList<ArrayList<String>> strandSeqs = null;	
		ArrayList<Boolean> contSCFlexVals = new ArrayList<Boolean>(Arrays.asList(true, true, true));
		ArrayList<String> pfImplVals = new ArrayList<String>(Arrays.asList(PFAbstract.getImpl(), PFAbstract.getImpl(), PFAbstract.getImpl()));

		long begin = System.currentTimeMillis();

		if(strand2AllowedSeqs == null)
			throw new RuntimeException("ERROR: call init() method on this object before invoking run()");

		int numSeqs = strand2AllowedSeqs.get(Strand.COMPLEX).getNumSeqs();
		for( int i = 0; i < numSeqs; ++i ) {

			// wt is seq 0, mutants are others
			System.out.println("\nComputing K* for sequence " + i + "/" + 
					(numSeqs-1) + ": " + 
					arrayList1D2String(strand2AllowedSeqs.get(Strand.COMPLEX).getStrandSeqAtPos(i), " ") + "\n");

			// get sequences
			strandSeqs = getStrandStringsAtPos(i);

			// create partition functions
			ConcurrentHashMap<Integer, PFAbstract> pfs = createPFsForSeq(strandSeqs, contSCFlexVals, pfImplVals);

			// create K* calculation for sequence
			KSCalc seq = new KSCalc(i, pfs);

			// store wtSeq
			if(i == 0) wtKSCalc = seq;

			// compute partition functions
			seq.run(wtKSCalc);

			if(wtKSCalc.getEpsilonStatus() != EApproxReached.TRUE)
				throw new RuntimeException("ERROR: could not compute the wild-type sequence to an epsilon value of "
						+ PFAbstract.targetEpsilon + ". Change the value of epsilon." );
			
			// compute K* scores and print output if all 
			// partition functions are computed to epsilon accuracy
			if( seq.getEpsilonStatus() == EApproxReached.TRUE ) {
				seq.printSummary( getOputputFileName() );
			}
		}

		// print statistics
		System.out.println("\nK* calculations computed: " + numSeqs);
		System.out.println("K* conformations minimized: " + countMinimizedConfs());
		System.out.println("Total # of conformations in search space: " + countTotNumConfs());
		System.out.println("K* running time: " + (System.currentTimeMillis()-begin)/1000 + " seconds\n");

		// peace the fuck out ^_^
	}
}
