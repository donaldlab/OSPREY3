package edu.duke.cs.osprey.kstar;

import edu.duke.cs.osprey.confspace.Sequence;

public class ScoredSequence {

    public final Sequence sequence;
    public final KStarScore score;

    public ScoredSequence(Sequence sequence, KStarScore score) {
        this.sequence = sequence;
        this.score = score;
    }

    @Override
    public String toString() {
        return "sequence: " + sequence + "   K*(log10): " + score;
    }

    public String toString(Sequence wildtype) {
        return "sequence: " + sequence.toString(Sequence.Renderer.AssignmentMutations) + "   K*(log10): " + score;
    }
}
