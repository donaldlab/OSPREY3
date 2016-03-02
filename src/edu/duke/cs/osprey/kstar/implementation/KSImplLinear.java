package edu.duke.cs.osprey.kstar.implementation;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.concurrent.ForkJoinPool;
import java.util.stream.IntStream;

import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.kstar.AllowedSeqs;
import edu.duke.cs.osprey.kstar.KSAbstract;
import edu.duke.cs.osprey.kstar.KSCalc;
import edu.duke.cs.osprey.kstar.Strand;
import edu.duke.cs.osprey.kstar.pfunction.PFAbstract;
import edu.duke.cs.osprey.kstar.pfunction.PFAbstract.EApproxReached;
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

		// get search problems per strand
		strand2AllSearchProblem.put(Strand.COMPLEX, cfp.getSearchProblem(Strand.COMPLEX, strand2AllowedSeqs.get(Strand.COMPLEX)));
		strand2AllSearchProblem.put(Strand.PROTEIN, cfp.getSearchProblem(Strand.PROTEIN, strand2AllowedSeqs.get(Strand.PROTEIN)));
		strand2AllSearchProblem.put(Strand.LIGAND, cfp.getSearchProblem(Strand.LIGAND, strand2AllowedSeqs.get(Strand.LIGAND)));

		printSequences();
		
		createEmatDir();
		
		createEnergyMatrices(true);
		// createEnergyMatrices(false);
	}
	
	
	public void createEnergyMatrices( boolean contSCFlex ) {

		System.out.println("\nCreating all energy matrices\n");

		long begin = System.currentTimeMillis();

		try {

			ForkJoinPool forkJoinPool = new ForkJoinPool(ThreadParallelism.getNumThreads());
			forkJoinPool.submit(() ->

			IntStream.range(0, strand2AllowedSeqs.get(Strand.COMPLEX).getNumSeqs()).parallel().forEach(i -> {

				System.out.println("\nCreating search problem for sequence " + 
						i + "/" + strand2AllowedSeqs.get(Strand.COMPLEX).getNumSeqs() + "\n");

				HashMap<Integer, SearchProblem> map = createSearchProblemsForSeq(i, contSCFlex);

				// put partition function in list, so we can parallelize energy matrix computation
				for(int strand : map.keySet()) {
					addSPToTmpList(strand, map.get(strand));
				}

			})).get();

			// create last of the energy matrices
			loadEnergyMatrices();

			// empty energy matrix list
			allSPNames.clear();

			System.out.println("\nFinished creating all energy matrices");
			System.out.println("Running time: " + (System.currentTimeMillis()-begin)/1000 + " seconds\n");

		} catch (Exception e) {
			System.out.println(e.getMessage());
			e.printStackTrace();
			System.exit(1);
		}
	}
	
	
	protected HashMap<Integer, SearchProblem> createSearchProblemsForSeq(int i, boolean contSCFlex) {
		// used to precompute energy matrices
		HashMap<Integer, SearchProblem> ans = new HashMap<Integer, SearchProblem>();

		ArrayList<Integer> strands = new ArrayList<>();
		strands.add(Strand.COMPLEX);
		strands.add(Strand.PROTEIN);
		strands.add(Strand.LIGAND);

		strands.parallelStream().forEach((strand) -> {

			AllowedSeqs strandSeqs = strand2AllowedSeqs.get(strand);

			ArrayList<String> seq = strandSeqs.getStrandSeqAtPos(i);

			String spName = getSearchProblemName(contSCFlex, strand, seq);

			if( createSP(spName) ) {

				SearchProblem seqSearchProblem = createSingleSequenceSearchProblem( contSCFlex, strand, seq );
				
				// synchronized
				addSPToLocalMap(strand, seqSearchProblem, ans);
			}		
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
		boolean[] contSCFlexVals = { true, true, true };
		
		long begin = System.currentTimeMillis();

		if(strand2AllowedSeqs == null)
			throw new RuntimeException("ERROR: call init() method on this object before invoking run()");

		int numSeqs = strand2AllowedSeqs.get(Strand.COMPLEX).getNumSeqs();
		for( int i = 0; i < numSeqs; ++i ) {

			// wt is seq 0, mutants are others
			System.out.println("\nComputing K* for sequence " + i + "/" + 
					(numSeqs-1) + ": " + 
					arrayList1D2String(strand2AllowedSeqs.get(Strand.COMPLEX).getStrandSeqAtPos(i), " ") + "\n");
			
			// create partition functions
			HashMap<Integer, PFAbstract> pfs = createPartitionFunctionsForSeq(i, contSCFlexVals);
			
			// create K* calculation for sequence
			KSCalc seq = new KSCalc(i, pfs);
			
			// store wtSeq
			if(i == 0) wtKSCalc = seq;
			
			// compute partition functions
			seq.run(wtKSCalc);
			
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
