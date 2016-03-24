package edu.duke.cs.osprey.kstar;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.concurrent.ConcurrentHashMap;

import edu.duke.cs.osprey.kstar.pfunc.PFAbstract;
import edu.duke.cs.osprey.kstar.pfunc.PFAbstract.EApproxReached;

/*
 * TODO
 * 1) at leaf, run standard k*
 * 2) change checkpoint to a per-object variable? or enable/disable it for fully
 * defined sequences
 */

public class KUStarNode {

	private static KSCalc wt;
	private static HashMap<Integer, AllowedSeqs> strand2AllowedSeqs;
	private static KSAbstract ksObj;
	private static int numCreated = 0;
	private static int numExpanded = 0;

	public KSCalc ub;
	public KSCalc lb;
	protected double lbScore;
	protected double ubScore;
	protected boolean scoreNeedsRefinement;


	public static void init( KSAbstract ksObj, 
			HashMap<Integer, AllowedSeqs> strand2AllowedSeqs, KSCalc wt ) {

		KUStarNode.ksObj = ksObj;
		KUStarNode.strand2AllowedSeqs = strand2AllowedSeqs;
		KUStarNode.wt = wt;
	}


	public KUStarNode( KSCalc lb, KSCalc ub, boolean scoreNeedsRefinement ) {
		this.lb = lb;
		this.ub = ub;
		lbScore = Double.MIN_VALUE;
		ubScore = Double.MAX_VALUE;
		this.scoreNeedsRefinement = scoreNeedsRefinement;
	}


	public static int getNumCreated() {
		return numCreated;
	}


	public static int getNumExpanded() {
		return numExpanded;
	}


	// only expand if scoreNeedsRefinement
	public ArrayList<KUStarNode> expand(boolean useTightBounds) {

		ArrayList<KUStarNode> children = null;

		if(!scoreNeedsRefinement()) {
			// there is a single sequence for which we are running k* with energy minimization			
			computeScore(this);

			children = new ArrayList<>();
			children.add(this);

			return children;
		}

		numExpanded++;

		// using complex, p, l convention
		ArrayList<Integer> strands = new ArrayList<>(Arrays.asList(Strand.LIGAND, Strand.PROTEIN, Strand.COMPLEX));

		ArrayList<Integer> nextDepths = new ArrayList<>(); 
		for( int i = 0; i < strands.size(); ++i ) nextDepths.add(0);

		ArrayList<ArrayList<String>> seqs = new ArrayList<>(); 
		for( int i = 0; i < strands.size(); ++i ) seqs.add(new ArrayList<String>());

		AllowedSeqs pSeqs = strand2AllowedSeqs.get(Strand.PROTEIN);
		AllowedSeqs lSeqs = strand2AllowedSeqs.get(Strand.LIGAND);

		// get next depths for each strand
		for( int strand : strands ) {

			if( currentDepth() != 0 ) {
				ArrayList<String> seq = lb.getPF(strand).getSequence();
				seqs.set(strand, seq);

				if( strand != Strand.COMPLEX )
					nextDepths.set(strand, Math.min(seq.size()+1, strand2AllowedSeqs.get(strand).getStrandSubSeqsMaxDepth()));
			}

			else {
				AllowedSeqs strandSeqs = strand2AllowedSeqs.get(strand);
				// go to the first non-zero depth
				int depth = 1;

				// set next depth
				if( strand != Strand.COMPLEX ) {
					// prepare subsequences
					strandSeqs.getStrandSubSeqsAtDepth( depth );
					nextDepths.set(strand, depth);
				}

				else
					// prepare subsequences
					strandSeqs.getStrandSubSeqsAtDepth( depth, pSeqs, lSeqs );
			}

			if( strand == Strand.COMPLEX )
				nextDepths.set(strand, nextDepths.get(Strand.PROTEIN) + nextDepths.get(Strand.LIGAND));
		}

		// get all sequences at next depth
		HashSet<ArrayList<String>> nextPSeqs = new HashSet<>(strand2AllowedSeqs.get(Strand.PROTEIN).getStrandSubSeqsAtDepth(nextDepths.get(Strand.PROTEIN)));
		HashSet<ArrayList<String>> nextLSeqs = new HashSet<>(strand2AllowedSeqs.get(Strand.LIGAND).getStrandSubSeqsAtDepth(nextDepths.get(Strand.LIGAND)));
		HashSet<ArrayList<String>> nextPLSeqs = new HashSet<>(strand2AllowedSeqs.get(Strand.COMPLEX).getStrandSubSeqsAtDepth(nextDepths.get(Strand.COMPLEX)));

		// in general, we want to add a single residue to the current protein and ligand
		// so our children are valid combinations of nextProtein x nextLigand.
		// however, from the root node, all reachable nodes are valid.

		if( currentDepth() != 0 ) {
			// reduce next L to those subsequences reachable from this
			ArrayList<String> currentLSeq = seqs.get(Strand.LIGAND);
			for( Iterator<ArrayList<String>> iterator = nextLSeqs.iterator(); iterator.hasNext(); ) {
				ArrayList<String> subSeq = iterator.next();	    
				if( !subSeq.subList(0, currentLSeq.size()).equals(currentLSeq) ) iterator.remove();
			}

			// reduce next P to those subsequences reachable from this
			ArrayList<String> currentPSeq = seqs.get(Strand.PROTEIN);
			for( Iterator<ArrayList<String>> iterator = nextPSeqs.iterator(); iterator.hasNext(); ) {
				ArrayList<String> subSeq = iterator.next();
				if( !subSeq.subList(0, currentPSeq.size()).equals(currentPSeq) ) iterator.remove();
			}
		}

		children = getPutativeChildren( nextPLSeqs, nextPSeqs, nextLSeqs, useTightBounds );

		if(children.size() > 0) {

			computeScores(children);

			// remove children whose upper bound epsilon values are not false
			for( Iterator<KUStarNode> iterator = children.iterator(); iterator.hasNext(); ) {

				KUStarNode child = iterator.next();

				if( child.scoreNeedsRefinement() && child.ub.getEpsilonStatus() != EApproxReached.FALSE ) {
					// epsilon is not going to be true, since we do not compute complex
					// epsilon cannot be not possible or not stable
					iterator.remove();
				}
			}
		}

		return children;
	}


	private void computeScores(ArrayList<KUStarNode> children) {

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
			for(KUStarNode child : children) {
				
				if(!child.isFullyProcessed()) 
					computeScore(child);
			}
		}
	}


	private void computeScore(KUStarNode child) {

		if(child.scoreNeedsRefinement()) {

			PFAbstract.suppressOutput = true;

			// compute protein and ligand upper bound pfs; prune if not stable
			// for stable pfs, compute lower bound
			child.ub.runPF(child.ub.getPF(Strand.LIGAND), wt.getPF(Strand.LIGAND), true);
			child.ub.runPF(child.ub.getPF(Strand.PROTEIN), wt.getPF(Strand.PROTEIN), true);

			if(child.ub.getEpsilonStatus() == EApproxReached.FALSE) {
				child.lb.run(wt);
				child.lbScore = -1.0 * child.lb.getKStarScore();
			}

			PFAbstract.suppressOutput = false;
		}

		else {
			KSAbstract.doCheckpoint = true;

			// we process leaf nodes as streams (in the complex case, at least)
			child.lb.run(wt);
			child.lbScore = -1.0 * child.lb.getKStarScore();

			KSAbstract.doCheckpoint = false;
		}
	}


	private boolean isUnique(ArrayList<KUStarNode> children) {		
		ArrayList<ArrayList<String>> listPL = new ArrayList<>();
		ArrayList<ArrayList<String>> listP = new ArrayList<>();
		ArrayList<ArrayList<String>> listL = new ArrayList<>();

		for(int i = 0; i < children.size(); ++i) {
			KUStarNode child = children.get(i);
			listPL.add(child.lb.getPF(Strand.COMPLEX).getSequence());
			listP.add(child.lb.getPF(Strand.PROTEIN).getSequence());
			listL.add(child.lb.getPF(Strand.LIGAND).getSequence());
		}

		HashSet<ArrayList<String>> setPL = new HashSet<>(listPL);
		HashSet<ArrayList<String>> setP = new HashSet<>(listP);
		HashSet<ArrayList<String>> setL = new HashSet<>(listL);

		if(listPL.size() != setPL.size() || listP.size() != setP.size() || listL.size() != setL.size())
			return false;

		return true;
	}


	private ArrayList<KUStarNode> getPutativeChildren( HashSet<ArrayList<String>> nextPLSeqs,
			HashSet<ArrayList<String>> pSeqs, HashSet<ArrayList<String>> lSeqs,
			boolean useTightBounds ) {

		ArrayList<ArrayList<String>> strandSeqs = new ArrayList<>(Arrays.asList(null, null, null));
		// lower bound pf values
		ArrayList<Boolean> lbContSCFlexVals = new ArrayList<>(Arrays.asList(false, false, true));
		ArrayList<String> lbPFImplVals = new ArrayList<>(Arrays.asList("trad", "trad", "new00"));

		// upper bound pf values
		ArrayList<Boolean> ubContSCFlexVals = new ArrayList<>(Arrays.asList(true, true, false));
		ArrayList<String> ubPFImplVals = new ArrayList<>(Arrays.asList("new00", "new00", "trad"));

		// minimized pf values
		ArrayList<Boolean> tightContSCFlexVals = new ArrayList<>(Arrays.asList(true, true, true));
		ArrayList<String> tightPFImplVals = new ArrayList<>(Arrays.asList(PFAbstract.getCFGImpl(), PFAbstract.getCFGImpl(), PFAbstract.getCFGImpl()));

		ArrayList<KUStarNode> ans = new ArrayList<>();

		for( ArrayList<String> pSeq : pSeqs ) {

			for( ArrayList<String> lSeq : lSeqs ) {

				ArrayList<String> putativeNextPLSeq = new ArrayList<String>();

				putativeNextPLSeq.addAll(pSeq);
				putativeNextPLSeq.addAll(lSeq);
				putativeNextPLSeq.trimToSize();

				if( !nextPLSeqs.contains(putativeNextPLSeq) ) continue;

				// create partition functions for next sequences
				strandSeqs.set(Strand.COMPLEX, putativeNextPLSeq);
				strandSeqs.set(Strand.PROTEIN, pSeq);
				strandSeqs.set(Strand.LIGAND, lSeq);

				if( childScoreNeedsRefinement() ) {
					// create partition functions
					ConcurrentHashMap<Integer, PFAbstract> lbPFs = ksObj.createPFs4Seq(strandSeqs, lbContSCFlexVals, lbPFImplVals);
					ConcurrentHashMap<Integer, PFAbstract> ubPFs = ksObj.createPFs4Seq(strandSeqs, ubContSCFlexVals, ubPFImplVals);

					numCreated++;

					// create KUStar node with lower and upper bounds
					ans.add( new KUStarNode( new KSCalc(numCreated, lbPFs), new KSCalc(numCreated, ubPFs), childScoreNeedsRefinement() ) );
				}

				else if( useTightBounds ) {

					// create a leaf node; we don't need upper or lower bounds ;we only need the minimized partition 
					// function our search problems exist, so we need only delete the lb, ub pfs from the table
					ArrayList<Integer> strands = new ArrayList<Integer>(Arrays.asList(Strand.LIGAND, 
							Strand.PROTEIN, Strand.COMPLEX));

					for(int strand : strands) {

						boolean contSCFlex = tightContSCFlexVals.get(strand);
						ArrayList<String> seq = strandSeqs.get(strand);

						String spName = ksObj.getSearchProblemName(contSCFlex, strand, seq);
						ksObj.removeFromMap(spName, false, true);
					}

					ConcurrentHashMap<Integer, PFAbstract> tightPFs = ksObj.createPFs4Seq(strandSeqs, tightContSCFlexVals, tightPFImplVals);

					numCreated++;

					// create new KUStar node with tight score
					ans.add( new KUStarNode( new KSCalc(numCreated, tightPFs), null, childScoreNeedsRefinement() ) );
				}

				else {
					// we will not use tight bounds, so the current node is re-designated a leaf node
					scoreNeedsRefinement = false;
					ans.add(this);
				}
			}
		}

		return ans;
	}


	private int currentDepth() {
		if( lb == null ) return 0;
		// number of assigned residues is depth
		return lb.getPF(Strand.COMPLEX).getSequence().size();
	}


	public boolean scoreNeedsRefinement() {
		return scoreNeedsRefinement;
	}
	
	
	public double getLBScore() {
		return lbScore;
	}
	
	
	public double getUBScore() {
		
		if(ub == null) return ubScore;
		
		ub.runPF(ub.getPF(Strand.COMPLEX), null, true);
		return ubScore = -1.0 * ub.getKStarScore();
	}


	private boolean childScoreNeedsRefinement() {

		int maxDepth = wt.getPF(Strand.COMPLEX).getSequence().size();

		int depth = currentDepth();

		if( depth < maxDepth )
			return true;

		return false;
	}


	public boolean isFullyProcessed() {
		return !scoreNeedsRefinement() && lb.getEpsilonStatus() == EApproxReached.TRUE;
	}


	//Comparator anonymous class implementation
	public static Comparator<KUStarNode> KUStarNodeComparator = new Comparator<KUStarNode>() {

		@Override
		public int compare( KUStarNode lhs, KUStarNode rhs ) {
			return lhs.lbScore >= rhs.lbScore ? 1 : -1;
		}
	};
}
