package edu.duke.cs.osprey.kstar;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.concurrent.ConcurrentHashMap;

import edu.duke.cs.osprey.kstar.pfunc.PFAbstract;
import edu.duke.cs.osprey.kstar.pfunc.PFAbstract.EApproxReached;
import edu.duke.cs.osprey.kstar.pfunc.impl.PFnew00;

/*
 * TODO
 * 1) only do upper bound stability check for depth 1
 */

public class KAStarNode {

	private static KSCalc wt;
	private static HashMap<Integer, AllowedSeqs> strand2AllowedSeqs;
	private static KSAbstract ksObj;
	private static int numCreated = 0;
	private static int numExpanded = 0;
	private static int numPruned = 0;
	public static int numProcessed = 0;

	public KSCalc ub;
	public KSCalc lb;
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
		lbScore = Double.MIN_VALUE;
		ubScore = Double.MAX_VALUE;
		this.scoreNeedsRefinement = scoreNeedsRefinement;
	}


	public static int getNumProcessed() {
		return numProcessed;
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
	public ArrayList<KAStarNode> expand(boolean useTightBounds) {

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
		ArrayList<Integer> strands = new ArrayList<>(Arrays.asList(Strand.LIGAND, Strand.PROTEIN, Strand.COMPLEX));

		ArrayList<Integer> nextDepths = new ArrayList<>(); 
		for( int i = 0; i < strands.size(); ++i ) nextDepths.add(0);

		ArrayList<ArrayList<String>> seqs = new ArrayList<>(); 
		for( int i = 0; i < strands.size(); ++i ) seqs.add(new ArrayList<String>());

		AllowedSeqs pSeqs = strand2AllowedSeqs.get(Strand.PROTEIN);
		AllowedSeqs lSeqs = strand2AllowedSeqs.get(Strand.LIGAND);

		// get next depths for each strand
		for( int strand : strands ) {

			if( depth() != 0 ) {
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

		if( depth() != 0 ) {
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

		children = getChildren( nextPLSeqs, nextPSeqs, nextLSeqs, useTightBounds );

		if(children.size() > 0) {

			computeScores(children);

			// remove children whose upper bound epsilon values are not false
			for( Iterator<KAStarNode> iterator = children.iterator(); iterator.hasNext(); ) {

				KAStarNode child = iterator.next();

				if( child.scoreNeedsRefinement() && child.ub.getEpsilonStatus() != EApproxReached.FALSE ) {

					// remove sequences that contain no valid conformations
					for( int strand : Arrays.asList(Strand.LIGAND, Strand.PROTEIN) ) {
						PFAbstract pf = child.ub.getPF(strand);

						if( pf.getEpsilonStatus() == EApproxReached.NOT_POSSIBLE && pf.getNumUnPruned().compareTo(BigInteger.ZERO) != 0 ) {
							System.out.println("ERROR: " + KSAbstract.list1D2String( pf.getSequence(), " " ) + " is wrong!!!");
						}
						
						if( pf.getEpsilonStatus() == EApproxReached.NOT_STABLE || 
								pf.getNumUnPruned().compareTo(BigInteger.ZERO) == 0 ) {

							ArrayList<String> seq = pf.getSequence();

							for( int strand2 : Arrays.asList(strand, Strand.COMPLEX) ) {
								pruneSequences(seq, strand2, nextDepths.get(strand2));
							}
						}

					}

					// epsilon is not going to be true, since we do not compute complex
					// epsilon cannot be not possible or not stable
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


	private void computeScores(ArrayList<KAStarNode> children) {

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


	private void computeScore(KAStarNode child) {

		if( child.scoreNeedsRefinement() ) {

			PFAbstract.suppressOutput = true;

			// compute protein and ligand upper bound pfs; prune if not stable
			// for stable pfs, compute lower bound
			for( int strand : Arrays.asList(Strand.LIGAND, Strand.PROTEIN) ) {
				if( child.ub.getEpsilonStatus() == EApproxReached.FALSE ) {
					child.ub.runPF(child.ub.getPF(strand), wt.getPF(strand), true, true);
				}
			}

			if(child.ub.getEpsilonStatus() != EApproxReached.FALSE)
				return;

			child.lb.run(wt, true, false);
			child.lbScore = -1.0 * child.lb.getKStarScoreLog10();

			PFAbstract.suppressOutput = false;
		}

		else {
			KSAbstract.doCheckPoint = true;

			// we process leaf nodes as streams (in the complex case, at least)
			child.lb.run(wt, false, true);
			child.lbScore = -1.0 * child.lb.getKStarScoreLog10();

			KSAbstract.doCheckPoint = false;
		}
	}


	private boolean isUnique(ArrayList<KAStarNode> children) {		
		ArrayList<ArrayList<String>> listPL = new ArrayList<>();
		ArrayList<ArrayList<String>> listP = new ArrayList<>();
		ArrayList<ArrayList<String>> listL = new ArrayList<>();

		for(int i = 0; i < children.size(); ++i) {
			KAStarNode child = children.get(i);
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


	private ArrayList<KAStarNode> getChildren( HashSet<ArrayList<String>> nextPLSeqs,
			HashSet<ArrayList<String>> pSeqs, HashSet<ArrayList<String>> lSeqs, boolean useTightBounds ) {

		ArrayList<ArrayList<String>> strandSeqs = new ArrayList<>(Arrays.asList(null, null, null));
		// lower bound pf values
		ArrayList<Boolean> lbContSCFlexVals = new ArrayList<>(Arrays.asList(false, false, true));
		ArrayList<String> lbPFImplVals = new ArrayList<>(Arrays.asList(PFAbstract.getCFGImpl(), PFAbstract.getCFGImpl(), PFnew00.getImpl()));

		// upper bound pf values
		ArrayList<Boolean> ubContSCFlexVals = new ArrayList<>(Arrays.asList(true, true, false));
		ArrayList<String> ubPFImplVals = new ArrayList<>(Arrays.asList(PFnew00.getImpl(), PFnew00.getImpl(), PFAbstract.getCFGImpl()));

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
				strandSeqs.set(Strand.COMPLEX, putativeNextPLSeq);
				strandSeqs.set(Strand.PROTEIN, pSeq);
				strandSeqs.set(Strand.LIGAND, lSeq);

				if( !isFullyDefined() ) {
					// create partition functions
					ConcurrentHashMap<Integer, PFAbstract> lbPFs = ksObj.createPFs4Seqs(strandSeqs, lbContSCFlexVals, lbPFImplVals);
					ConcurrentHashMap<Integer, PFAbstract> ubPFs = ksObj.createPFs4Seqs(strandSeqs, ubContSCFlexVals, ubPFImplVals);

					// create KUStar node with lower and upper bounds
					ans.add( new KAStarNode( new KSCalc(numCreated, lbPFs), new KSCalc(numCreated, ubPFs), childScoreNeedsRefinement() ) );
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
						ksObj.removeFromMap(spName, true, true);
					}

					ConcurrentHashMap<Integer, PFAbstract> tightPFs = ksObj.createPFs4Seqs(strandSeqs, tightContSCFlexVals, tightPFImplVals);

					// create new KUStar node with tight score
					ans.add( new KAStarNode( new KSCalc(numCreated, tightPFs), null, childScoreNeedsRefinement() ) );

					// processed nodes are those whose confs will be minimized
					numProcessed++;
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

		PFAbstract.suppressOutput = true;
		//ub.run(wt, false, true);
		ub.runPF(ub.getPF(Strand.COMPLEX), null, true, false);
		PFAbstract.suppressOutput = false;

		return ubScore = -1.0 * ub.getKStarScoreLog10();
	}


	private boolean isFullyDefined() {
		int maxDepth = wt.getPF(Strand.COMPLEX).getSequence().size();

		if( depth() < maxDepth )
			return false;

		return true;
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
