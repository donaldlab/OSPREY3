package edu.duke.cs.osprey.kstar;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
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
	static int numExpanded = 0;

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

		// using complex, p, l convention
		ArrayList<Integer> strands = new ArrayList<Integer>(Arrays.asList(Strand.COMPLEX, Strand.PROTEIN, Strand.LIGAND));

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
				nextDepths.set(strand, Math.min(seq.size()+1, strand2AllowedSeqs.get(strand).getStrandSubSeqsMaxDepth()));
			}

			else {
				AllowedSeqs strandSeqs = strand2AllowedSeqs.get(strand);
				// go to the first non-zero depth
				for(int depth = 0; depth <= strandSeqs.getStrandSubSeqsMaxDepth(); ++depth) {

					HashSet<ArrayList<String>> subSeqsAtDepth = strand == 
							Strand.COMPLEX ? strandSeqs.getStrandSubSeqsAtDepth( depth, pSeqs, lSeqs ) : 
								strandSeqs.getStrandSubSeqsAtDepth( depth );

							if( subSeqsAtDepth.size() != 0 ) {
								nextDepths.set(strand, depth);
								break;
							}
				}
			}
		}

		// get all sequences at next depth
		HashSet<ArrayList<String>> nextPSeqs = strand2AllowedSeqs.get(Strand.PROTEIN).getStrandSubSeqsAtDepth(nextDepths.get(Strand.PROTEIN));
		HashSet<ArrayList<String>> nextLSeqs = strand2AllowedSeqs.get(Strand.LIGAND).getStrandSubSeqsAtDepth(nextDepths.get(Strand.LIGAND));
		HashSet<ArrayList<String>> nextPLSeqs = new HashSet<ArrayList<String>>(strand2AllowedSeqs.get(Strand.COMPLEX).getStrandSubSeqsAtDepth(nextDepths.get(Strand.COMPLEX)));

		ArrayList<KUStarNode> successors;

		// in general, we want to add a single residue to either the protein xor the ligand
		// so our successor is nextProtein+thisLigand, thisProtein+nextLigand
		// however, to expand the root node, we add two residues to get nextProtein+nextLigand

		if( currentDepth() == 0 ) {
			successors = getPutativeSuccessors( nextPLSeqs, nextPSeqs, nextLSeqs );
		}

		else {
			HashSet<ArrayList<String>> currentPSeq = new HashSet<>();
			currentPSeq.add(seqs.get(Strand.PROTEIN));
			successors = getPutativeSuccessors( nextPLSeqs, currentPSeq, nextLSeqs );

			HashSet<ArrayList<String>> currentLSeq = new HashSet<>();
			currentLSeq.add(seqs.get(Strand.LIGAND));
			successors.addAll( getPutativeSuccessors( nextPLSeqs, nextPSeqs, currentLSeq ) );
		}

		// compute all p and l UPPER bound partition functions
		// 		if not stable, then do not add

		return successors;
	}


	private ArrayList<KUStarNode> getPutativeSuccessors( HashSet<ArrayList<String>> nextPLSeqs,
			HashSet<ArrayList<String>> pSeqs,
			HashSet<ArrayList<String>> lSeqs ) {

		ArrayList<ArrayList<String>> strandSeqs = new ArrayList<ArrayList<String>>(Arrays.asList(null, null, null));
		ArrayList<Boolean> contSCFlexVals = new ArrayList<Boolean>(Arrays.asList(true, false, false));
		ArrayList<String> pfImplVals = new ArrayList<String>(Arrays.asList("1nubnm", "1nnocache", "1nnocache"));

		ArrayList<KUStarNode> ans = new ArrayList<>();

		for( ArrayList<String> pSeq : pSeqs ) {

			for( ArrayList<String> lSeq : lSeqs ) {

				ArrayList<String> putativeNextPLSeq = new ArrayList<String>();

				putativeNextPLSeq.addAll(pSeq);
				putativeNextPLSeq.addAll(lSeq);

				if( !nextPLSeqs.contains(putativeNextPLSeq) ) continue;

				// create partition functions for next sequences
				strandSeqs.set(Strand.COMPLEX, putativeNextPLSeq);
				strandSeqs.set(Strand.PROTEIN, pSeq);
				strandSeqs.set(Strand.LIGAND, lSeq);

				// create partition functions
				ConcurrentHashMap<Integer, PFAbstract> pfs = ksObj.createPFsForSeq(strandSeqs, contSCFlexVals, pfImplVals);

				// create KUStar node
				ans.add( new KUStarNode( new KSCalc(++numExpanded, pfs), childScoreNeedsRefinement() ) );
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
