package edu.duke.cs.osprey.kstar.impl;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.kstar.AllowedSeqs;
import edu.duke.cs.osprey.kstar.KSAbstract;
import edu.duke.cs.osprey.kstar.KAStarNode;
import edu.duke.cs.osprey.kstar.KAStarTree;
import edu.duke.cs.osprey.kstar.Strand;
import edu.duke.cs.osprey.kstar.pfunc.PFAbstract;
import edu.duke.cs.osprey.kstar.pfunc.impl.PFnew00;

public class KSImplKAStar extends KSAbstract {

	public static boolean useTightBounds = true;

	public KSImplKAStar(ConfigFileParser cfp) {
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
	protected void preLoadPFs( ArrayList<Boolean> contSCFlexVals ) {

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

									String pfImpl = null;
									boolean flex = false;

									// here we use contscflex == true to mean upper bound and
									// == false means lower bound 

									if(contSCFlex) {
										// do upper bound k* calc
										if(strand == Strand.COMPLEX) { flex = false; pfImpl = PFAbstract.getCFGImpl(); }
										else { flex = true; pfImpl = new PFnew00().getImpl(); }
									}

									else {
										// do lower bound k* calc
										if(strand == Strand.COMPLEX) { flex = true; pfImpl = new PFnew00().getImpl(); }
										else { flex = false; pfImpl = PFAbstract.getCFGImpl(); }
									}

									PFAbstract pf = createPF4Seq(flex, strand, seq, pfImpl);

									name2PF.put(pf.getSearchProblemName(), pf);
								});
					}
				}
			}

			// create last of the energy matrices
			loadAndPruneMatricesFromPFMap();

		} catch (Exception e) {
			System.out.println(e.getMessage());
			e.printStackTrace();
			System.exit(1);
		}
	}


	@Override
	public String getKSMethod() {
		return "kastar";
	}


	@Override
	public void run() {

		setStartTime(System.currentTimeMillis());

		int numOutput = 0;

		if(useTightBounds)
			numOutput = runTB();

		else
			numOutput = runLB();

		System.out.println("\nseqsCompleted: " + KAStarNode.getNumLeavesCompleted() 
			+ " seqsCreated: " + KAStarNode.getNumLeavesCreated()
			+ " seqsPossible: " + strand2AllowedSeqs.get(Strand.COMPLEX).getNumSeqs()
			+ " seqsInOutput: " + numOutput);
		
		System.out.println("K* running time: " + (System.currentTimeMillis()-getStartTime())/1000 + " seconds\n");

		abortPFs();
	}


	private int runTB() {
		
		int completed = 0;
		
		// compute wt sequence for reference
		wtKSCalc = computeWTCalc();
		
		if( strand2AllowedSeqs.get(Strand.COMPLEX).getNumSeqs() <= 1 )
			return completed; // wt is sequence[0]
		
		// initialize KUStar tree
		KAStarTree tree = new KAStarTree(this, strand2AllowedSeqs, wtKSCalc);

		// add root node
		tree.add( new KAStarNode(null, null, true) );

		int target = cfp.getParams().getInt("KStarNumSeqs", 5);

		for( KAStarNode best = tree.poll(); best != null && completed < target; 
				best = tree.poll() ) {

			if( best.isFullyProcessed() ) {

				best.checkConsistency(best);

				best.lb.printSummary( getOputputFilePath(), getStartTime(), getNumSeqsCompleted(0) );
				completed++;
				continue;
			}

			ArrayList<KAStarNode> children = best.expand(useTightBounds);
			tree.add(children);
		}

		return completed;
	}


	private int runLB() {
		// run until the lower bound of the next completed sequence is greater 
		// than the upper bound of any previously completed sequence

		int completed = 0;
		
		// compute wt sequence for reference
		wtKSCalc = computeWTCalc();
		
		if( strand2AllowedSeqs.get(Strand.COMPLEX).getNumSeqs() <= 1 )
			return completed; // wt is sequence[0]

		// initialize KUStar tree
		KAStarTree tree = new KAStarTree(this, strand2AllowedSeqs, wtKSCalc);

		// add root node
		tree.add( new KAStarNode(null, null, true) );

		double gUB = Double.NEGATIVE_INFINITY;

		for( KAStarNode best = tree.poll(); best != null; best = tree.poll() ) {

			if( best.isFullyProcessed() ) {

				best.lb.printSummary( getOputputFilePath(), getStartTime(), getNumSeqsCompleted(0) );

				double bestUB = best.getUBScore();

				if(completed++ == 0) gUB = bestUB;

				else if(best.getLBScore() > gUB) break;

				else if(bestUB > gUB) gUB = bestUB;

				continue;
			}

			ArrayList<KAStarNode> children = best.expand(useTightBounds);
			tree.add(children);
		}

		return completed;
	}
}
