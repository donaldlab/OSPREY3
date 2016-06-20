package edu.duke.cs.osprey.kstar;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.concurrent.ConcurrentHashMap;
import java.util.stream.IntStream;

import edu.duke.cs.osprey.kstar.pfunc.PFAbstract;
import edu.duke.cs.osprey.kstar.pfunc.PFAbstract.EApproxReached;
import edu.duke.cs.osprey.kstar.pfunc.impl.PFTrad;
import edu.duke.cs.osprey.kstar.pfunc.impl.PFnew00;

/*
 * TODO
 * 1) only do upper bound stability check for depth 1?
 * 2) am i re-computing unbound pfs for tight bound leaf?
 */

public class KAStarNode {

	private static KSCalc wt;
	private static HashMap<Integer, AllowedSeqs> strand2AllowedSeqs;
	private static KSAbstract ksObj;
	private static int numCreated = 0;
	private static int numExpanded = 0;
	private static int numPruned = 0;
	public static int numLeavesCreated = 0;
	public static int numLeavesCompleted = 0;

	public KSCalc ub;
	public KSCalc lb;
	protected double plbScore; //parent lower bound
	protected double lbScore;
	protected double ubScore;
	protected boolean scoreNeedsRefinement;


	public static void init( KSAbstract ksObj, 
			HashMap<Integer, AllowedSeqs> strand2AllowedSeqs, KSCalc wt ) {

		KAStarNode.ksObj = ksObj;
		KAStarNode.strand2AllowedSeqs = strand2AllowedSeqs;
		KAStarNode.wt = wt;
	}


	public KAStarNode( KSCalc lb, KSCalc ub, boolean scoreNeedsRefinement ) {
		this.lb = lb;
		this.ub = ub;
		plbScore = Double.NEGATIVE_INFINITY;
		lbScore = Double.NEGATIVE_INFINITY;
		ubScore = Double.POSITIVE_INFINITY;
		this.scoreNeedsRefinement = scoreNeedsRefinement;
	}


	public static int getNumLeavesCreated() {
		return numLeavesCreated;
	}
	
	
	public static int getNumLeavesCompleted() {
		return numLeavesCompleted;
	}


	public static int getNumCreated() {
		return numCreated;
	}


	public static int getNumExpanded() {
		return numExpanded;
	}


	public static int getNumPruned() {
		return numPruned;
	}


	// only expand if scoreNeedsRefinement
	public ArrayList<KAStarNode> expand( boolean useTightBounds ) {

		ArrayList<KAStarNode> children = null;

		if( isFullyDefined() ) {

			if( !useTightBounds ) {
				// we will not use tight bounds, so the current node is re-designated a leaf node
				scoreNeedsRefinement = false;
			}

			if( !scoreNeedsRefinement() ) {

				if( !isFullyProcessed() ) {
					// we will use tight bounds; compute tight bound		
					computeScore(this);
				}

				children = new ArrayList<>();
				children.add(this);
				return children;
			}
		}

		numExpanded++;

		// using complex, p, l convention
		ArrayList<Integer> strands = new ArrayList<>(Arrays.asList(Termini.LIGAND, Termini.PROTEIN, Termini.COMPLEX));

		ArrayList<Integer> nextDepths = new ArrayList<>(); 
		for( int i = 0; i < strands.size(); ++i ) nextDepths.add(0);

		ArrayList<ArrayList<String>> seqs = new ArrayList<>(); 
		for( int i = 0; i < strands.size(); ++i ) seqs.add(new ArrayList<String>());

		AllowedSeqs pSeqs = strand2AllowedSeqs.get(Termini.PROTEIN);
		AllowedSeqs lSeqs = strand2AllowedSeqs.get(Termini.LIGAND);

		// get next depths for each strand
		for( int strand : strands ) {

			if( depth() != 0 ) {
				ArrayList<String> seq = lb.getPF(strand).getSequence();
				seqs.set(strand, seq);

				if( strand != Termini.COMPLEX )
					nextDepths.set(strand, Math.min(seq.size()+1, strand2AllowedSeqs.get(strand).getStrandSubSeqsMaxDepth()));
			}

			else {
				AllowedSeqs strandSeqs = strand2AllowedSeqs.get(strand);
				// go to the first non-zero depth
				int depth = 1;

				// set next depth
				if( strand != Termini.COMPLEX ) {
					// prepare subsequences
					strandSeqs.getStrandSubSeqsAtDepth( depth );
					nextDepths.set(strand, depth);
				}

				else
					// prepare subsequences
					strandSeqs.getStrandSubSeqsAtDepth( depth, pSeqs, lSeqs );
			}

			if( strand == Termini.COMPLEX )
				nextDepths.set(strand, nextDepths.get(Termini.PROTEIN) + nextDepths.get(Termini.LIGAND));
		}

		// get all sequences at next depth
		HashSet<ArrayList<String>> nextPSeqs = new HashSet<>(strand2AllowedSeqs.get(Termini.PROTEIN).getStrandSubSeqsAtDepth(nextDepths.get(Termini.PROTEIN)));
		HashSet<ArrayList<String>> nextLSeqs = new HashSet<>(strand2AllowedSeqs.get(Termini.LIGAND).getStrandSubSeqsAtDepth(nextDepths.get(Termini.LIGAND)));
		HashSet<ArrayList<String>> nextPLSeqs = new HashSet<>(strand2AllowedSeqs.get(Termini.COMPLEX).getStrandSubSeqsAtDepth(nextDepths.get(Termini.COMPLEX)));

		// in general, we want to add a single residue to the current protein and ligand
		// so our children are valid combinations of nextProtein x nextLigand.
		// however, from the root node, all reachable nodes are valid.

		if( depth() != 0 ) {
			// reduce next L to those subsequences reachable from this
			ArrayList<String> currentLSeq = seqs.get(Termini.LIGAND);
			for( Iterator<ArrayList<String>> iterator = nextLSeqs.iterator(); iterator.hasNext(); ) {
				ArrayList<String> subSeq = iterator.next();	    
				if( !subSeq.subList(0, currentLSeq.size()).equals(currentLSeq) ) iterator.remove();
			}

			// reduce next P to those subsequences reachable from this
			ArrayList<String> currentPSeq = seqs.get(Termini.PROTEIN);
			for( Iterator<ArrayList<String>> iterator = nextPSeqs.iterator(); iterator.hasNext(); ) {
				ArrayList<String> subSeq = iterator.next();
				if( !subSeq.subList(0, currentPSeq.size()).equals(currentPSeq) ) iterator.remove();
			}
		}

		children = getChildren( nextPLSeqs, nextPSeqs, nextLSeqs, useTightBounds );

		if(children.size() > 0) {

			//computeScoresSimple(children);
			computeScoresParallel(children);

			// remove children whose upper bound epsilon values are not false
			for( Iterator<KAStarNode> iterator = children.iterator(); iterator.hasNext(); ) {

				KAStarNode child = iterator.next();

				if( child.scoreNeedsRefinement() && child.ub.getEpsilonStatus() != EApproxReached.FALSE ) {

					// remove sequences that contain no valid conformations
					for( int strand : Arrays.asList(Termini.LIGAND, Termini.PROTEIN) ) {
						PFAbstract pf = child.ub.getPF(strand);

						if( pf.getEpsilonStatus() == EApproxReached.NOT_POSSIBLE && pf.getNumUnPruned().compareTo(BigInteger.ZERO) != 0 ) {
							System.out.println("ERROR: " + KSAbstract.list1D2String( pf.getSequence(), " " )  + " " + pf.getFlexibility() + " is wrong!!!");
						}

						if( pf.getEpsilonStatus() == EApproxReached.NOT_STABLE || 
								pf.getNumUnPruned().compareTo(BigInteger.ZERO) == 0 ) {

							ArrayList<String> seq = pf.getSequence();

							for( int strand2 : Arrays.asList(strand, Termini.COMPLEX) ) {
								pruneSequences(seq, strand2, nextDepths.get(strand2));
							}
						}
					}

					// epsilon is not going to be true, since we do not compute complex
					// epsilon cannot be not possible or not stable
					iterator.remove();
				}

				else if( !child.scoreNeedsRefinement() && !child.lb.canContinue() ) {
					iterator.remove();
				}
			}
		}

		return children;
	}


	private void pruneSequences( ArrayList<String> seq, int strand, int depth ) {
		HashSet<ArrayList<String>> set = null;
		int maxDepth = strand2AllowedSeqs.get(strand).getStrandSubSeqsMaxDepth();

		for( ; depth <= maxDepth; ++depth) {
			set = strand2AllowedSeqs.get(strand).getStrandSubSeqsAtDepth(depth);
			int oldSize = set.size();

			AllowedSeqs.deleteFromSet(seq, set);
			numPruned += (oldSize - set.size());
		}
	}


	private void computeScore(KAStarNode child) {

		if( child.scoreNeedsRefinement() ) {

			PFAbstract.suppressOutput = true;

			// compute protein and ligand upper bound pfs; prune if not stable
			// for stable pfs, compute lower bound
			for( int strand : Arrays.asList(Termini.LIGAND, Termini.PROTEIN) ) {
				if( child.ub.getEpsilonStatus() == EApproxReached.FALSE ) {
					child.ub.runPF(child.ub.getPF(strand), wt.getPF(strand), true, true);
				}
			}
			
			if(child.ub.getEpsilonStatus() != EApproxReached.FALSE) {
				PFAbstract.suppressOutput = false;
				return;
			}

			child.lb.run(wt, true, false);
			child.lbScore = -1.0 * child.lb.getKStarScoreLog10(true);

			checkConsistency(child);

			PFAbstract.suppressOutput = false;
		}

		else {
			KSAbstract.doCheckPoint = true;

			// we process leaf nodes as streams (in the complex case, at least)
			child.lb.run(wt, false, true);
			
			if( !child.lb.canContinue() ) return; // epsilon is not possible or not stable
			
			if( child.lb.getEpsilonStatus() == EApproxReached.TRUE ) {
				numLeavesCompleted = ksObj.getNumSeqsCompleted(1);
			}
			
			child.lbScore = -1.0 * child.lb.getKStarScoreLog10(true);

			KSAbstract.doCheckPoint = false;
		}
	}
	

	private class CalcParams {
		int strand; boolean complete; boolean stabilityCheck;
		public CalcParams( int strand, boolean complete, boolean stabilityCheck ) {
			this.strand = strand; this.complete = complete; this.stabilityCheck = stabilityCheck;
		}
	}

	
	protected void computeScoresSimple(ArrayList<KAStarNode> children) {

		if(children.size() == 0) return;

		boolean parallel = false;

		if(isUnique(children) && children.size() > 1)
			parallel = true;

		if(parallel) {
			children.parallelStream().forEach( child -> {

				if(!child.isFullyProcessed()) 
					computeScore(child);
			});
		}

		else {
			for(KAStarNode child : children) {

				if(!child.isFullyProcessed()) 
					computeScore(child);
			}
		}
	}
	

	protected void computeScoresParallel( ArrayList<KAStarNode> children ) {

		if( children.size() == 0 ) return;

		if( children.get(0).scoreNeedsRefinement() ) {
			PFAbstract.suppressOutput = true;

			// get upper bound protein and ligands
			ArrayList<KSCalc> calcs = new ArrayList<>();
			ArrayList<CalcParams> params = new ArrayList<>();

			for( int strand : Arrays.asList(Termini.LIGAND, Termini.PROTEIN) ) {
				ArrayList<KSCalc> ans = getCalcs4Strand( children, false, strand );
				for(KSCalc calc : ans) {
					calcs.add(calc);
					params.add(new CalcParams(strand, true, true));
				}
			}
			// run all pfs
			IntStream.range(0, calcs.size()).parallel().forEach(i -> {
				KSCalc calc = calcs.get(i);
				CalcParams param = params.get(i);
				calc.runPF(calc.getPF(param.strand), wt.getPF(param.strand), param.complete, param.stabilityCheck);
			});

			// clear calcs, params
			calcs.clear();
			params.clear();

			ArrayList<KAStarNode> children2 = new ArrayList<>();

			// select viable children
			// retain viable children for pruning
			for( KAStarNode child : children ) {
				if(child.ub.getEpsilonStatus() == EApproxReached.FALSE)
					children2.add(child);
			}

			// for viable children, run all pfs
			for( int strand : Arrays.asList(Termini.LIGAND, Termini.PROTEIN, Termini.COMPLEX) ) {
				ArrayList<KSCalc> ans = getCalcs4Strand( children2, true, strand );
				for(KSCalc calc : ans) {
					calcs.add(calc);
					params.add(new CalcParams(strand, true, false));
				}
			}
			IntStream.range(0, calcs.size()).parallel().forEach(i -> {
				KSCalc calc = calcs.get(i);
				CalcParams param = params.get(i);
				calc.runPF(calc.getPF(param.strand), wt.getPF(param.strand), param.complete, param.stabilityCheck);
			});

			// set scores for surviving children
			for( KAStarNode child : children2 ) {
				child.lbScore = -1.0 * child.lb.getKStarScoreLog10(true);

				checkConsistency(child);
			}

			PFAbstract.suppressOutput = false;
		}

		else {
			for( KAStarNode child : children ) 
				computeScore(child);
		}

	}


	private ArrayList<KSCalc> getCalcs4Strand( ArrayList<KAStarNode> children, boolean lb, int strand ) {

		// get minimal set of kscalcs that covers all desired pfs 
		HashMap<ArrayList<String>, KSCalc> seq2KSCalc = new HashMap<>();

		for( KAStarNode child : children ) {
			KSCalc calc = lb ? child.lb : child.ub;
			PFAbstract pf = calc.getPF(strand);
			seq2KSCalc.put(pf.getSequence(), calc);
		}

		ArrayList<KSCalc> calcs = new ArrayList<>(seq2KSCalc.values());
		return calcs;
	}


	private boolean isUnique(ArrayList<KAStarNode> children) {		
		ArrayList<ArrayList<String>> listPL = new ArrayList<>();
		ArrayList<ArrayList<String>> listP = new ArrayList<>();
		ArrayList<ArrayList<String>> listL = new ArrayList<>();

		for(int i = 0; i < children.size(); ++i) {
			KAStarNode child = children.get(i);
			listPL.add(child.lb.getPF(Termini.COMPLEX).getSequence());
			listP.add(child.lb.getPF(Termini.PROTEIN).getSequence());
			listL.add(child.lb.getPF(Termini.LIGAND).getSequence());
		}

		HashSet<ArrayList<String>> setPL = new HashSet<>(listPL);
		HashSet<ArrayList<String>> setP = new HashSet<>(listP);
		HashSet<ArrayList<String>> setL = new HashSet<>(listL);

		if(listPL.size() != setPL.size() || listP.size() != setP.size() || listL.size() != setL.size())
			return false;

		return true;
	}


	private ArrayList<KAStarNode> getChildren( HashSet<ArrayList<String>> nextPLSeqs,
			HashSet<ArrayList<String>> pSeqs, HashSet<ArrayList<String>> lSeqs, boolean useTightBounds ) {

		ArrayList<ArrayList<String>> strandSeqs = new ArrayList<>(Arrays.asList(null, null, null));
		// lower bound pf values
		ArrayList<Boolean> lbContSCFlexVals = new ArrayList<>(Arrays.asList(false, false, true));
		ArrayList<String> lbPFImplVals = new ArrayList<>(Arrays.asList(new PFTrad().getImpl(), new PFTrad().getImpl(), new PFnew00().getImpl()));

		// upper bound pf values
		ArrayList<Boolean> ubContSCFlexVals = new ArrayList<>(Arrays.asList(true, true, false));
		ArrayList<String> ubPFImplVals = new ArrayList<>(Arrays.asList(new PFnew00().getImpl(), new PFnew00().getImpl(), new PFTrad().getImpl()));

		// minimized pf values
		ArrayList<Boolean> tightContSCFlexVals = new ArrayList<>(Arrays.asList(true, true, true));
		ArrayList<String> tightPFImplVals = new ArrayList<>(Arrays.asList(PFAbstract.getCFGImpl(), PFAbstract.getCFGImpl(), PFAbstract.getCFGImpl()));

		ArrayList<KAStarNode> ans = new ArrayList<>();

		for( ArrayList<String> pSeq : pSeqs ) {

			for( ArrayList<String> lSeq : lSeqs ) {

				ArrayList<String> putativeNextPLSeq = new ArrayList<String>();

				putativeNextPLSeq.addAll(pSeq);
				putativeNextPLSeq.addAll(lSeq);
				putativeNextPLSeq.trimToSize();

				if( !nextPLSeqs.contains(putativeNextPLSeq) ) continue;

				numCreated++;

				// create partition functions for next sequences
				strandSeqs.set(Termini.COMPLEX, putativeNextPLSeq);
				strandSeqs.set(Termini.PROTEIN, pSeq);
				strandSeqs.set(Termini.LIGAND, lSeq);

				if( !isFullyDefined() ) {
					// create partition functions
					ConcurrentHashMap<Integer, PFAbstract> lbPFs = ksObj.createPFs4Seqs(strandSeqs, lbContSCFlexVals, lbPFImplVals);
					ConcurrentHashMap<Integer, PFAbstract> ubPFs = ksObj.createPFs4Seqs(strandSeqs, ubContSCFlexVals, ubPFImplVals);

					// create KUStar node with lower and upper bounds
					ans.add( new KAStarNode( new KSCalc(numCreated, lbPFs), new KSCalc(numCreated, ubPFs), childScoreNeedsRefinement() ) );
					ans.get(ans.size()-1).plbScore = this.lbScore;
				}

				else if( useTightBounds ) {
					// create a leaf node; we don't need upper or lower bounds ;we only need the minimized partition 
					// function our search problems exist, so we need only delete the lb, ub pfs from the table
					
					ConcurrentHashMap<Integer, PFAbstract> tightPFs = ksObj.createPFs4Seqs(strandSeqs, tightContSCFlexVals, tightPFImplVals);

					// assign sequence number from allowedSequences obj
					AllowedSeqs complexSeqs = ksObj.strand2AllowedSeqs.get(Termini.COMPLEX);
					int seqID = complexSeqs.getPosOfSeq(tightPFs.get(Termini.COMPLEX).getSequence());
					
					// create new KUStar node with tight score
					ans.add( new KAStarNode( new KSCalc(seqID, tightPFs), null, childScoreNeedsRefinement() ) );
					ans.get(ans.size()-1).plbScore = this.lbScore;

					// processed nodes are those whose confs will be minimized
					numLeavesCreated++;
				}

				else
					throw new RuntimeException("ERROR: cannot expand a fully assigned node");
			}
		}

		return ans;
	}


	private int depth() {
		if( lb == null ) return 0;
		// number of assigned residues is depth
		return lb.getPF(Termini.COMPLEX).getSequence().size();
	}


	public boolean scoreNeedsRefinement() {
		return scoreNeedsRefinement;
	}


	public double getParentLBScore() {
		return plbScore;
	}


	public double getLBScore() {
		return lbScore;
	}


	public double getUBScore() {

		if(ub == null) return ubScore;

		PFAbstract.suppressOutput = true;
		//ub.run(wt, false, true);
		ub.runPF(ub.getPF(Termini.COMPLEX), null, true, false);
		PFAbstract.suppressOutput = false;

		return ubScore = -1.0 * ub.getKStarScoreLog10(true);
	}


	private boolean isFullyDefined() {
		int maxDepth = wt.getPF(Termini.COMPLEX).getSequence().size();

		if( depth() < maxDepth )
			return false;

		return true;
	}


	public void checkConsistency(KAStarNode node) {
		// score must be > lb
		double nodeLB = node.getLBScore(), parentLB = node.getParentLBScore();
		
		if(nodeLB == Double.NEGATIVE_INFINITY)
			return;

		if( parentLB > nodeLB ) 
			throw new RuntimeException("ERROR: parentLB: " + parentLB + " must be <= nodeLB: " + nodeLB);
	}
	
	
	private boolean childScoreNeedsRefinement() {
		if( !isFullyDefined() ) return true;

		else if(scoreNeedsRefinement()) return false;

		else throw new RuntimeException("ERROR: cannot expand a leaf node");
	}


	public boolean isFullyProcessed() {
		return !scoreNeedsRefinement() && lb.getEpsilonStatus() == EApproxReached.TRUE;
	}


	//Comparator anonymous class implementation
	public static Comparator<KAStarNode> KUStarNodeComparator = new Comparator<KAStarNode>() {

		@Override
		public int compare( KAStarNode lhs, KAStarNode rhs ) {
			return lhs.lbScore >= rhs.lbScore ? 1 : -1;
		}
	};
}
