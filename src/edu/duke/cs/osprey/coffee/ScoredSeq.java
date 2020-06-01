package edu.duke.cs.osprey.coffee;

import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.tools.MathTools;


public class ScoredSeq {

	public final Sequence sequence;
	public final MathTools.DoubleBounds lmfe;
	public final MathTools.DoubleBounds[] stateFreeEnergies;

	public ScoredSeq(Sequence sequence, MathTools.DoubleBounds lmfe, MathTools.DoubleBounds[] stateFreeEnergies) {
		this.sequence = sequence;
		this.lmfe = lmfe;
		this.stateFreeEnergies = stateFreeEnergies;
	}
}
