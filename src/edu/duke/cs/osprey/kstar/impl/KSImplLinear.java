package edu.duke.cs.osprey.kstar.impl;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
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
import edu.duke.cs.osprey.tools.ObjectIO;

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

		if(doCheckpoint)
			createCheckPointDir();

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

		if(strand2AllowedSeqs == null)
			throw new RuntimeException("ERROR: call init() method on this object before invoking run()");

		long begin = System.currentTimeMillis();

		if(doCheckpoint)
			runRR();

		else
			runFCFS();

		// print statistics
		System.out.println("\nK* calculations computed: " + strand2AllowedSeqs.get(Strand.COMPLEX).getNumSeqs());
		System.out.println("K* conformations minimized: " + countMinimizedConfs());
		System.out.println("Total # of conformations in search space: " + countTotNumConfs());
		System.out.println("K* running time: " + (System.currentTimeMillis()-begin)/1000 + " seconds\n");
	}


	protected void runFCFS() {

		// each value corresponds to the desired flexibility of the 
		// pl, p, and l conformation spaces, respectively
		ArrayList<ArrayList<String>> strandSeqs = null;	
		ArrayList<Boolean> contSCFlexVals = new ArrayList<Boolean>(Arrays.asList(true, true, true));
		ArrayList<String> pfImplVals = new ArrayList<String>(Arrays.asList(PFAbstract.getImpl(), PFAbstract.getImpl(), PFAbstract.getImpl()));

		wtKSCalc = computeWTCalc();

		int numSeqs = strand2AllowedSeqs.get(Strand.COMPLEX).getNumSeqs();
		for( int i = 1; i < numSeqs; ++i ) {

			// wt is seq 0, mutants are others
			System.out.println("\nComputing K* for sequence " + i + "/" + 
					(numSeqs-1) + ": " + 
					arrayList1D2String(strand2AllowedSeqs.get(Strand.COMPLEX).getStrandSeqAtPos(i), " ") + "\n");

			// get sequences
			strandSeqs = getStrandStringsAtPos(i);

			// create partition functions
			ConcurrentHashMap<Integer, PFAbstract> pfs = createPFsForSeq(strandSeqs, contSCFlexVals, pfImplVals);

			// create K* calculation for sequence
			KSCalc calc = new KSCalc(i, pfs);

			// compute partition functions
			calc.run(wtKSCalc);

			// compute K* scores and print output if all 
			// partition functions are computed to epsilon accuracy
			if( calc.getEpsilonStatus() == EApproxReached.TRUE ) {
				calc.printSummary( getOputputFilePath(), false );
			}
		}
	}


	protected void runRR() {
		
		// each value corresponds to the desired flexibility of the 
		// pl, p, and l conformation spaces, respectively
		ArrayList<ArrayList<String>> strandSeqs = null;	
		ArrayList<Boolean> contSCFlexVals = new ArrayList<Boolean>(Arrays.asList(true, true, true));
		ArrayList<String> pfImplVals = new ArrayList<String>(Arrays.asList(PFAbstract.getImpl(), PFAbstract.getImpl(), PFAbstract.getImpl()));

		// get all sequences		
		@SuppressWarnings("unchecked")
		HashSet<ArrayList<String>> seqSet = new HashSet<ArrayList<String>>((ArrayList<ArrayList<String>>) 
				ObjectIO.deepCopy(strand2AllowedSeqs.get(Strand.COMPLEX).getStrandSeqList()));

		// remove completed seqs from the set of calculations we must compute
		seqSet.removeAll(getSeqsFromFile(getOputputFilePath()));

		// run wt
		wtKSCalc = computeWTCalc();

		// remove completed sequences
		seqSet.remove(wtKSCalc.getPF(Strand.COMPLEX).getSequence());		

		int numSeqs = strand2AllowedSeqs.get(Strand.COMPLEX).getNumSeqs();

		do {

			for( int i = 1; i < numSeqs; ++i ) {

				// wt is seq 0, mutants are others
				ArrayList<String> seq = strand2AllowedSeqs.get(Strand.COMPLEX).getStrandSeqAtPos(i);

				if( !seqSet.contains(seq) ) continue;

				System.out.println("\nResuming K* for sequence " + i + ": " + arrayList1D2String(seq, " ") + "\n");

				// get sequences
				strandSeqs = getStrandStringsAtPos(i);

				// create partition functions
				ConcurrentHashMap<Integer, PFAbstract> pfs = createPFsForSeq(strandSeqs, contSCFlexVals, pfImplVals);

				// create K* calculation for sequence
				KSCalc calc = new KSCalc(i, pfs);

				// compute partition functions
				calc.run(wtKSCalc, KSAbstract.checkpointInterval);

				// serialize P and L; should only happen once
				for( int strand : Arrays.asList(Strand.LIGAND, Strand.PROTEIN) ) {
					PFAbstract pf = calc.getPF(strand);
					if(!pf.checkPointExists()) calc.serializePF(strand);
				}

				// remove entry from checkpoint file
				PFAbstract pf = calc.getPF(Strand.COMPLEX);
				calc.deleteSeqFromFile( pf.getSequence(), getCheckPointFilePath() );

				if( calc.getEpsilonStatus() != EApproxReached.FALSE ) {
					// we are not going to checkpoint this, so clean up
					calc.deleteCheckPointFile(Strand.COMPLEX);
					seqSet.remove(seq);

					if( calc.getEpsilonStatus() == EApproxReached.TRUE ) {
						calc.printSummary( getOputputFilePath(), false );
					}
				}

				else {
					// remove partition funtion from memory, write checkpoint
					name2PF.remove(pf.getSearchProblem().name);
					calc.serializePF(Strand.COMPLEX);
					calc.printSummary( getCheckPointFilePath(), false );
				}
			}

		} while(seqSet.size() > 0);

		// delete checkpoint dir and checkpoint file
		ObjectIO.delete(getCheckPointDir());
	}


}
