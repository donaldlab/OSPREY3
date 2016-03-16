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
import edu.duke.cs.osprey.tools.ExpFunction;

/*
 * TODO
 * 1) set score
 * 2) when is this node fully assigned?
 * 3) when does this score need refinement?
 */

public class KUStarNode {

	private static KSCalc wt;
	private static HashMap<Integer, AllowedSeqs> strand2AllowedSeqs;
	private static KSAbstract ksObj;
	private static int numCreated = 0;
	private static int numExpanded = 0;

	private KSCalc calc;
	protected double score;
	protected boolean scoreNeedsRefinement;


	public static void init( KSAbstract ksObj, 
			HashMap<Integer, AllowedSeqs> strand2AllowedSeqs, KSCalc wt ) {

		KUStarNode.ksObj = ksObj;
		KUStarNode.strand2AllowedSeqs = strand2AllowedSeqs;
		KUStarNode.wt = wt;
	}


	public KUStarNode( KSCalc calc, boolean scoreNeedsRefinement ) {
		this.calc = calc;
		score = Double.MAX_VALUE;
		this.scoreNeedsRefinement = scoreNeedsRefinement;
	}


	public static int getNumCreated() {
		return numCreated;
	}


	public static int getNumExpanded() {
		return numExpanded;
	}


	protected void setScore() {

		calc.run(wt);

		if( calc.getEpsilonStatus() == EApproxReached.TRUE ) {

			PFAbstract pl = calc.getPF(Strand.COMPLEX);
			PFAbstract p = calc.getPF(Strand.PROTEIN);
			PFAbstract l = calc.getPF(Strand.LIGAND);

			ExpFunction e = new ExpFunction();

			score = e.log10(p.getQStar()) + e.log10(l.getQStar()) - e.log10(pl.getQStar());
		}

		else
			score = Double.MAX_VALUE;
	}


	// only expand if scoreNeedsRefinement
	public ArrayList<KUStarNode> expand() {

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
				ArrayList<String> seq = calc.getPF(strand).getSequence();
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

		ArrayList<KUStarNode> successors = null;

		// in general, we want to add a single residue to either the protein xor the ligand
		// so our successor is nextProtein+thisLigand, thisProtein+nextLigand
		// however, to expand the root node, we add two residues to get nextProtein+nextLigand

		if( currentDepth() == 0 ) {
			successors = getPutativeSuccessors( nextPLSeqs, nextPSeqs, nextLSeqs );
		}

		else {
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

			successors = getPutativeSuccessors( nextPLSeqs, nextPSeqs, nextLSeqs );
		}

		// parallel for?
		// compute all p and l UPPER bound partition functions
		// 		if not stable, then do not add?

		return successors;
	}


	private ArrayList<KUStarNode> getPutativeSuccessors( HashSet<ArrayList<String>> nextPLSeqs,
			HashSet<ArrayList<String>> pSeqs,
			HashSet<ArrayList<String>> lSeqs ) {

		ArrayList<ArrayList<String>> strandSeqs = new ArrayList<>(Arrays.asList(null, null, null));
		ArrayList<Boolean> contSCFlexVals = new ArrayList<>(Arrays.asList(false, false, true));
		ArrayList<String> pfImplVals = new ArrayList<>(Arrays.asList("1nnocache", "1nnocache", "1nubnm"));

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

				// create partition functions
				ConcurrentHashMap<Integer, PFAbstract> pfs = ksObj.createPFsForSeq(strandSeqs, contSCFlexVals, pfImplVals);

				// create KUStar node
				ans.add( new KUStarNode( new KSCalc(++numCreated, pfs), childScoreNeedsRefinement() ) );
			}
		}

		return ans;
	}


	protected int currentDepth() {
		if( calc == null ) return 0;
		// number of assigned residues is depth
		return calc.getPF(Strand.COMPLEX).getSequence().size();
	}


	public boolean scoreNeedsRefinement() {
		return scoreNeedsRefinement;
	}


	private boolean childScoreNeedsRefinement() {

		int maxDepth = wt.getPF(Strand.COMPLEX).getSequence().size();

		int depth = currentDepth();

		if( depth != maxDepth )
			return true;

		return false;
	}


	//Comparator anonymous class implementation
	public static Comparator<KUStarNode> KUStarNodeComparator = new Comparator<KUStarNode>() {

		@Override
		public int compare( KUStarNode lhs, KUStarNode rhs ) {
			return lhs.score >= rhs.score ? 1 : -1;
		}
	};
}
