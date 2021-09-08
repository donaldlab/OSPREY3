package edu.duke.cs.osprey.coffee.seqdb;

import edu.duke.cs.osprey.astar.conf.ConfAStarNode;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.tools.MathTools.DoubleBounds;

import java.util.Arrays;


public class SeqFreeEnergies {

	public final Sequence seq;
	public final DoubleBounds[] freeEnergies;

	public SeqFreeEnergies(Sequence seq, DoubleBounds[] freeEnergies) {
		this.seq = seq;
		this.freeEnergies = freeEnergies;
	}

	@Override
	public boolean equals(Object other) {
		return other instanceof SeqFreeEnergies && equals((SeqFreeEnergies)other);
	}

	public boolean equals(SeqFreeEnergies other) {
		return this.seq.equals(other.seq)
			&& Arrays.equals(this.freeEnergies, other.freeEnergies);
	}
}
