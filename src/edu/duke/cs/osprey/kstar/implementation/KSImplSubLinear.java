package edu.duke.cs.osprey.kstar.implementation;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.concurrent.ForkJoinPool;

import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.kstar.AllowedSeqs;
import edu.duke.cs.osprey.kstar.KSAbstract;
import edu.duke.cs.osprey.kstar.Strand;
import edu.duke.cs.osprey.parallelism.ThreadParallelism;

public class KSImplSubLinear extends KSAbstract {

	public KSImplSubLinear(ConfigFileParser cfp) {
		super(cfp);
		// TODO Auto-generated constructor stub
	}

	@Override
	public void init( HashMap<Integer, AllowedSeqs> strand2AllowedSeqs ) {

		this.strand2AllowedSeqs = strand2AllowedSeqs;

		// get search problems per strand
		strand2AllSearchProblem.put(Strand.COMPLEX, cfp.getSearchProblem(Strand.COMPLEX, strand2AllowedSeqs.get(Strand.COMPLEX)));
		strand2AllSearchProblem.put(Strand.PROTEIN, cfp.getSearchProblem(Strand.PROTEIN, strand2AllowedSeqs.get(Strand.PROTEIN)));
		strand2AllSearchProblem.put(Strand.LIGAND, cfp.getSearchProblem(Strand.LIGAND, strand2AllowedSeqs.get(Strand.LIGAND)));

		createEmatDir();
		
		createEnergyMatrices(true);
		createEnergyMatrices(false);
	}

	@Override
	public void run() {
		// TODO Auto-generated method stub

	}

	@Override
	public String getKSMethod() {
		return "sublinear";
	}

	@Override
	public void createEnergyMatrices( boolean contSCFlex ) {

		System.out.println("\nCreating all energy matrices\n");

		long begin = System.currentTimeMillis();

		try {

			int[] strands = { Strand.COMPLEX, Strand.PROTEIN, Strand.LIGAND };
			AllowedSeqs pSeqs = strand2AllowedSeqs.get(Strand.PROTEIN);
			AllowedSeqs lSeqs = strand2AllowedSeqs.get(Strand.LIGAND);

			for( int strand : strands ) {

				// create global search problem for strand
				String allSeqSPName = getSearchProblemName( contSCFlex, strand );

				if( createSP(allSeqSPName) ) {

					SearchProblem allSeqSearchProblem = createAllSequenceSearchProblem(contSCFlex, strand);

					addSPToTmpList(strand, allSeqSearchProblem);
				}

				AllowedSeqs seqs = strand2AllowedSeqs.get(strand);

				// ignore depth 0
				for( int depth = 1; depth <= seqs.getStrandSubSeqsMaxDepth(); ++depth ) {

					ArrayList<ArrayList<String>> subSeqsAtDepth = strand == Strand.COMPLEX ? 
							seqs.getStrandSubSeqsAtDepth( depth, pSeqs, lSeqs ) : seqs.getStrandSubSeqsAtDepth( depth );

							ForkJoinPool forkJoinPool = new ForkJoinPool(ThreadParallelism.getNumThreads());
							forkJoinPool.submit(() -> subSeqsAtDepth.parallelStream().forEach( subSeq -> {

								//for( ArrayList<String> subSeq : subSeqsAtDepth ) {

								String singleSeqSPName = getSearchProblemName( contSCFlex, strand, subSeq );

								if( createSP(singleSeqSPName) ) {

									SearchProblem seqSearchProblem = createSingleSequenceSearchProblem( contSCFlex, strand, subSeq );

									addSPToTmpList(strand, seqSearchProblem);
								}

								//}

							})).get();
				}
			}

			// create last of the energy matrices
			loadEnergyMatrices();

			// empty energy matrix list
			allSPNames.clear();

			System.out.println("\nFinished creating all energy matrices");
			System.out.println("Running time: " + ((System.currentTimeMillis()-begin)/1000) + " seconds\n");

		} catch (Exception e) {
			System.out.println(e.getMessage());
			e.printStackTrace();
			System.exit(1);
		}
	}

}
