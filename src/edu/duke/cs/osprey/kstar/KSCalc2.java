package edu.duke.cs.osprey.kstar;

import java.util.HashMap;

import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.kstar.PFAbstract.EApproxReached;

public class KSCalc2 {

	private HashMap<Integer, PFAbstract> strand2PF = new HashMap<Integer, PFAbstract>();

	public KSCalc2( PFAbstract pl, PFAbstract p, PFAbstract l) {

		strand2PF.put(Strand.COMPLEX, pl);
		strand2PF.put(Strand.PROTEIN, p);
		strand2PF.put(Strand.LIGAND, l);
	}

	protected PFAbstract getPF(int strand) {
		return strand2PF.get(strand);
	}

	protected boolean pIsStable(PFAbstract wtP) {
		PFAbstract p = strand2PF.get(Strand.PROTEIN);

		if(p.getUpperBound().compareTo( wtP.getQStar().multiply(PFAbstract.getStabilityThreshold()) ) >= 0)
			return true;

		return false;
	}
	
	protected EApproxReached getEpsilonStatus() {

		PFAbstract pl = strand2PF.get(Strand.COMPLEX);
		PFAbstract p = strand2PF.get(Strand.PROTEIN);
		PFAbstract l = strand2PF.get(Strand.LIGAND);

		if(pl.eAppx == EApproxReached.NOT_POSSIBLE 
				|| p.eAppx == EApproxReached.NOT_POSSIBLE 
				|| l.eAppx == EApproxReached.NOT_POSSIBLE) 
			return EApproxReached.NOT_POSSIBLE;
		
		if(p.eAppx == EApproxReached.NOT_STABLE)
			return EApproxReached.NOT_STABLE;

		else if(pl.eAppx == EApproxReached.TRUE
				&& p.eAppx == EApproxReached.TRUE
				&& l.eAppx == EApproxReached.TRUE) 
			return EApproxReached.TRUE;

		return EApproxReached.FALSE;
	}
}
