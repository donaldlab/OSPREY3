package edu.duke.cs.osprey.kstar;

import java.util.HashMap;

import edu.duke.cs.osprey.confspace.AllowedSeqs;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.kstar.PFAbstract.EApproxReached;

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

		createEnergyMatrices();
	}

	@Override
	public String getKSMethod() {
		return "linear";
	}

	@Override
	public void run() {

		long begin = System.currentTimeMillis();

		if(strand2AllowedSeqs == null)
			throw new RuntimeException("ERROR: call init() method on this object before invoking run()");

		int numSeqs = strand2AllowedSeqs.get(Strand.COMPLEX).getNumSeqs();
		for( int i = 0; i < numSeqs; ++i ) {

			// wt is seq 0, mutants are others
			System.out.println("\nComputing K* for sequence " + i + "/" + 
					(numSeqs-1) + ": " + 
					arrayList1D2String(strand2AllowedSeqs.get(Strand.COMPLEX).getStrandSeq(i), " "));
			
			// create partition functions
			HashMap<Integer, PFAbstract> pfs = createPartitionFunctionsForSeq(i);
			
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

		long end = System.currentTimeMillis();

		// print statistics
		System.out.println("\nK* calculations computed: " + numSeqs);
		System.out.println("K* conformations minimized: " + countMinimizedConfs());
		System.out.println("Total # of conformations in search space: " + countTotNumConfs());
		System.out.println("K* running time: " + (end-begin)/1000 + " seconds\n");

		// peace the fuck out ^_^
	}

}
