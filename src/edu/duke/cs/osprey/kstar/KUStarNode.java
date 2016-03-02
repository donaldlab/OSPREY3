package edu.duke.cs.osprey.kstar;

import edu.duke.cs.osprey.kstar.pfunction.PFAbstract;
import edu.duke.cs.osprey.kstar.pfunction.PFAbstract.EApproxReached;
import edu.duke.cs.osprey.tools.ExpFunction;

/*
 * TODO
 * 1) set score
 * 2) when is this node fully assigned?
 * 3) when does this score need refinement?
 */

public class KUStarNode {

	protected static KSCalc wt;
	
	private KSCalc calc;
	protected double score;
	protected boolean scoreNeedsRefinement;
	
	
	protected static void init( KSCalc wtCalc ) {
		wt = wtCalc;
	}
	
	
	public KUStarNode( boolean scoreNeedsRefinement ) {
		
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
	protected void expand() {
	}
	
	
	protected int depth() {
		// number of assigned residues is depth
		return calc.getPF(Strand.COMPLEX).getSequence().size();
	}
	
	
	private boolean childScoreNeedsRefinement() {
		
		int maxDepth = wt.getPF(Strand.COMPLEX).getSequence().size();
		
		int depth = depth();
		
		if( depth != maxDepth )
			return true;
		
		return false;
	}
}
