package edu.duke.cs.osprey.ematrix;

import edu.duke.cs.osprey.confspace.*;
import edu.duke.cs.osprey.tools.MathTools;

import java.util.List;

/**
 * Score corrector for quantities that are independent of assignments outside the tuple
 */
public class IndependentScoreCorrector extends ScoreCorrector<TupE> {
    public IndependentScoreCorrector(List<SimpleConfSpace.Position> positions, MathTools.Optimizer opt) {
        super(positions, opt);
    }

    @Override
    protected TupleTrieImplementations.TupETrie makeTrie(List<SimpleConfSpace.Position> positions) {
        return new TupleTrieImplementations.TupETrie(positions);
    }

    @Override
    protected double correctionSize(TupE correction) {
        // By default this value represents CorrectedScore - UnCorrectedScore
        return correction.E;
    }

}
