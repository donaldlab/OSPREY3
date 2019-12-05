package edu.duke.cs.osprey.energy.compiled;


import edu.duke.cs.osprey.confspace.compiled.ConfSpace;
import edu.duke.cs.osprey.confspace.compiled.PosInter;
import edu.duke.cs.osprey.confspace.compiled.PosInterDist;
import edu.duke.cs.osprey.ematrix.SimpleReferenceEnergies;

import java.util.List;


/**
 * Position Interaction Generator
 * Generates position interactions for conformation energy calculators.
 */
public class PosInterGen {

	public final PosInterDist dist;
	public final SimpleReferenceEnergies eref;

	public PosInterGen(PosInterDist dist, SimpleReferenceEnergies eref) {
		this.dist = dist;
		this.eref = eref;
	}

	public List<PosInter> single(ConfSpace confSpace, int posi1, int confi1) {
		return dist.single(confSpace, eref, posi1, confi1);
	}

	public List<PosInter> pair(ConfSpace confSpace, int posi1, int confi1, int posi2, int confi2) {
		return dist.pair(confSpace, eref, posi1, confi1, posi2, confi2);
	}

	public List<PosInter> all(ConfSpace confSpace, int[] conf) {
		return PosInterDist.all(confSpace, eref, conf);
	}

	public List<PosInter> dynamic(ConfSpace confSpace, int[] conf) {
		return PosInterDist.dynamic(confSpace, eref, conf);
	}
}
