package edu.duke.cs.osprey.kstar;

import edu.duke.cs.osprey.confspace.Sequence;

public record ScoredSequence(Sequence sequence, KStarScore score) {

    @Override
    public String toString() {
        return "sequence: " + sequence + "   K*(log10): " + score;
    }
}
