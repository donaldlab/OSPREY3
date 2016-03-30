package edu.duke.cs.osprey.kstar.impl;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.ConcurrentHashMap;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.kstar.AllowedSeqs;
import edu.duke.cs.osprey.kstar.KSAbstract;
import edu.duke.cs.osprey.kstar.KSCalc;
import edu.duke.cs.osprey.kstar.Strand;
import edu.duke.cs.osprey.kstar.pfunc.PFAbstract;
import edu.duke.cs.osprey.kstar.pfunc.PFAbstract.EApproxReached;
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

		if(doCheckpoint) createCheckPointDir();

		ArrayList<Boolean> contSCFlexVals = new ArrayList<Boolean>(Arrays.asList(true));
		createEmats(contSCFlexVals);
	}


	protected void preLoadPFs( ArrayList<Boolean> contSCFlexVals ) {

		try {
			
			for( boolean contSCFlex : contSCFlexVals ) {
				
				for( int strand : Arrays.asList( Strand.COMPLEX, Strand.PROTEIN, Strand.LIGAND ) ) {
					
					HashSet<ArrayList<String>> seqs = new HashSet<>(strand2AllowedSeqs.get(strand).getStrandSeqList());
					
					seqs.parallelStream().forEach(seq -> {
						
						PFAbstract pf = createPF4Seq(contSCFlex, strand, seq, PFAbstract.getCFGImpl());
						
						// put partition function in list, so we can parallelize energy matrix computation
						name2PF.put(pf.getSearchProblemName(), pf);
					});
				}
			}

			loadAndPruneMatricesFromPFMap(); 
			prunedSingleSeqs = true;

		} catch (Exception e) {
			System.out.println(e.getMessage());
			e.printStackTrace();
			System.exit(1);
		}
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
		ArrayList<String> pfImplVals = new ArrayList<String>(Arrays.asList(PFAbstract.getCFGImpl(), 
				PFAbstract.getCFGImpl(), PFAbstract.getCFGImpl()));

		wtKSCalc = computeWTCalc();

		int numSeqs = strand2AllowedSeqs.get(Strand.COMPLEX).getNumSeqs();
		for( int i = 1; i < numSeqs; ++i ) {

			// wt is seq 0, mutants are others
			System.out.println("\nComputing K* for sequence " + i + "/" + 
					(numSeqs-1) + ": " + 
					list1D2String(strand2AllowedSeqs.get(Strand.COMPLEX).getStrandSeqAtPos(i), " ") + "\n");

			// get sequences
			strandSeqs = getStrandStringsAtPos(i);

			// create partition functions
			ConcurrentHashMap<Integer, PFAbstract> pfs = createPFs4Seqs(strandSeqs, contSCFlexVals, pfImplVals);

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
		ArrayList<String> pfImplVals = new ArrayList<String>(Arrays.asList(PFAbstract.getCFGImpl(), 
				PFAbstract.getCFGImpl(), PFAbstract.getCFGImpl()));

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

				System.out.println("\nResuming K* for sequence " + i + ": " + list1D2String(seq, " ") + "\n");

				// get sequences
				strandSeqs = getStrandStringsAtPos(i);

				// create partition functions
				ConcurrentHashMap<Integer, PFAbstract> pfs = createPFs4Seqs(strandSeqs, contSCFlexVals, pfImplVals);

				// create K* calculation for sequence
				KSCalc calc = new KSCalc(i, pfs);

				// compute partition functions
				calc.run(wtKSCalc);

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
					name2PF.remove(pf.getSearchProblemName());
					calc.serializePF(Strand.COMPLEX);
					calc.printSummary( getCheckPointFilePath(), false );
				}
			}

		} while(seqSet.size() > 0);

		// delete checkpoint dir and checkpoint file
		ObjectIO.delete(getCheckPointDir());
	}


}
