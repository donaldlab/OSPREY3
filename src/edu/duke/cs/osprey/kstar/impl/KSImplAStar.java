package edu.duke.cs.osprey.kstar.impl;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.kstar.AllowedSeqs;
import edu.duke.cs.osprey.kstar.KSAbstract;
import edu.duke.cs.osprey.kstar.KUStarNode;
import edu.duke.cs.osprey.kstar.KUStarTree;
import edu.duke.cs.osprey.kstar.Strand;

public class KSImplAStar extends KSAbstract {

	public KSImplAStar(ConfigFileParser cfp) {
		super(cfp);
	}

	@Override
	public void init( HashMap<Integer, AllowedSeqs> strand2AllowedSeqs ) {

		this.strand2AllowedSeqs = strand2AllowedSeqs;

		createEmatDir();

		ArrayList<Boolean> contSCFlexVals = new ArrayList<Boolean>(Arrays.asList(true, false));
		createEmats(contSCFlexVals);
	}

	
	@Override
	protected void prepareAllSingleSeqSPs( ArrayList<Boolean> contSCFlexVals ) {

		try {

			ArrayList<Integer>strands = new ArrayList<Integer>(Arrays.asList(Strand.COMPLEX, Strand.PROTEIN, Strand.LIGAND));
			AllowedSeqs pSeqs = strand2AllowedSeqs.get(Strand.PROTEIN);
			AllowedSeqs lSeqs = strand2AllowedSeqs.get(Strand.LIGAND);

			for( boolean contSCFlex : contSCFlexVals ) {

				for( int strand : strands ) {

					AllowedSeqs seqs = strand2AllowedSeqs.get(strand);

					// ignore depth 0
					for( int depth = 1; depth <= seqs.getStrandSubSeqsMaxDepth(); ++depth ) {

						HashSet<ArrayList<String>> subSeqsAtDepth = strand == 
								Strand.COMPLEX ? seqs.getStrandSubSeqsAtDepth( depth, pSeqs, lSeqs ) : 
									seqs.getStrandSubSeqsAtDepth( depth );

								subSeqsAtDepth.parallelStream().forEach( seq -> {

									ArrayList<Integer> flexResIndexes = seqs.getFlexResIndexesFromSeq(seq);

									String spName = getSearchProblemName( contSCFlex, strand, seq );

									SearchProblem seqSP = null;			
									if( (seqSP = name2SP.get(spName)) == null ) {
										seqSP = createSingleSeqSP( contSCFlex, strand, seq, flexResIndexes, true );
									}

									name2SP.put(spName, seqSP);
								});
					}
				}
			}

			// create last of the energy matrices
			loadAndPruneMatrices();

		} catch (Exception e) {
			System.out.println(e.getMessage());
			e.printStackTrace();
			System.exit(1);
		}
	}
	
	
	@Override
	public String getKSMethod() {
		return "astar";
	}
	

	@Override
	public void run() {
		
		long begin = System.currentTimeMillis();
		
		// compute wt sequence for reference
		wtKSCalc = computeWTCalc();

		// initialize KUStar tree
		KUStarTree tree = new KUStarTree(this, strand2AllowedSeqs, wtKSCalc);

		// add root node
		tree.add( new KUStarNode(null, null, true) );

		int numSeqs = cfp.getParams().getInt("KStarNumSeqs", 5);
		int completed = 0;

		for( KUStarNode best = tree.poll(); best != null && completed < numSeqs; 
				best = tree.poll() ) {
			
			if( best.isFullyProcessed() ) {
				// run full k*: p and l completely, then pl as a stream
				completed++;
				continue;
			}

			ArrayList<KUStarNode> children = best.expand();
			tree.add(children);
		}

		System.out.println("completed: " + completed + " numExpanded: " + KUStarNode.getNumExpanded() 
			+ " numSubSeqs: " + strand2AllowedSeqs.get(Strand.COMPLEX).getNumSubSeqs()
			+ " numSeqs: " + strand2AllowedSeqs.get(Strand.COMPLEX).getNumSeqs());
		
		System.out.println("K* running time: " + (System.currentTimeMillis()-begin)/1000 + " seconds\n");
	}
}
