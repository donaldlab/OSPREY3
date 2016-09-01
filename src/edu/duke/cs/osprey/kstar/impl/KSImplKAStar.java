package edu.duke.cs.osprey.kstar.impl;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

import edu.duke.cs.osprey.kstar.KSAllowedSeqs;
import edu.duke.cs.osprey.kstar.KSAbstract;
import edu.duke.cs.osprey.kstar.KSCalc;
import edu.duke.cs.osprey.kstar.KSConfigFileParser;
import edu.duke.cs.osprey.kstar.KAStarNode;
import edu.duke.cs.osprey.kstar.KAStarTree;
import edu.duke.cs.osprey.kstar.KSTermini;

public class KSImplKAStar extends KSAbstract {

	public static boolean useTightBounds = true;
	public static String nodeExpansionMethod = "parallel1";
	protected BigInteger numMinimizedConfs = BigInteger.ZERO;
	protected HashSet<KSCalc> leafNodeKSCalcs = null; // nodes that get energy minimized

	public KSImplKAStar(KSConfigFileParser cfp) {
		super(cfp);
		leafNodeKSCalcs = new HashSet<>();
	}

	@Override
	public void init( HashMap<Integer, KSAllowedSeqs> strand2AllowedSeqs ) {

		this.strand2AllowedSeqs = strand2AllowedSeqs;

		printSequences();

		createOutputDir();
		
		createEmatDir();

		ArrayList<Boolean> contSCFlexVals = new ArrayList<Boolean>(Arrays.asList(true, false));
		createEmats(contSCFlexVals);
	}


	@Override
	public String getKSMethod() {
		return "kastar";
	}


	public void addLeafNode(KSCalc leaf) {
		leafNodeKSCalcs.add(leaf);
	}

	
	public boolean removeLeafNode(KSCalc leaf) {
		// also updates minimized confs
		for(int strand : Arrays.asList(KSTermini.LIGAND, KSTermini.PROTEIN, KSTermini.COMPLEX))
			numMinimizedConfs = numMinimizedConfs.add(leaf.getPF(strand).getNumMinimized4Output());
		
		return leafNodeKSCalcs.remove(leaf);
	}
	

	protected BigInteger countMinimizedConfs() {
		for(KSCalc leaf : leafNodeKSCalcs) {
			for(int strand : Arrays.asList(KSTermini.LIGAND, KSTermini.PROTEIN, KSTermini.COMPLEX))
				numMinimizedConfs = numMinimizedConfs.add(leaf.getPF(strand).getNumMinimized4Output());
		}

		return numMinimizedConfs;
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
		+ " seqsPossible: " + strand2AllowedSeqs.get(KSTermini.COMPLEX).getNumSeqs()
		+ " seqsInOutput: " + numOutput);

		System.out.println("K* leaf conformations processed: " + countMinimizedConfs());
		System.out.println("K* running time: " + (System.currentTimeMillis()-getStartTime())/1000 + " seconds\n");

		abortPFs();
	}


	private int runTB() {

		int completed = 0;

		// compute wt sequence for reference
		wtKSCalc = computeWTCalc();
		//setBestCalc(wtKSCalc);

		if( strand2AllowedSeqs.get(KSTermini.COMPLEX).getNumSeqs() <= 1 )
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
				
				//if(passesInterMutationPruning(best.lb))
				best.lb.printSummary( getOputputFilePath(), getStartTime(), getNumSeqsCreated(0), getNumSeqsCompleted(0) );
				
				completed++;
				
				if(!removeLeafNode(best.lb)) 
					throw new RuntimeException("ERROR: could not remove leaf node");
				
				continue;
			}

			ArrayList<KAStarNode> children = best.expand();
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

		if( strand2AllowedSeqs.get(KSTermini.COMPLEX).getNumSeqs() <= 1 )
			return completed; // wt is sequence[0]

		// initialize KUStar tree
		KAStarTree tree = new KAStarTree(this, strand2AllowedSeqs, wtKSCalc);

		// add root node
		tree.add( new KAStarNode(null, null, true) );

		double gUB = Double.NEGATIVE_INFINITY;

		for( KAStarNode best = tree.poll(); best != null; best = tree.poll() ) {

			if( best.isFullyProcessed() ) {

				double uB = best.getUBScore();

				if( completed++ == 0 ) gUB = uB;

				if( best.getLBScore() > gUB && gUB > Double.NEGATIVE_INFINITY ) 
					break;
				else
					best.lb.printSummary( getOputputFilePath(), getStartTime(), getNumSeqsCreated(0), getNumSeqsCompleted(0) );

				if( uB > gUB ) gUB = uB;

				continue;
			}

			ArrayList<KAStarNode> children = best.expand();
			tree.add(children);
		}

		return completed;
	}
}
