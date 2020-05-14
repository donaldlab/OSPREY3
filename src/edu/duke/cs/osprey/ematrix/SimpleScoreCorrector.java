package edu.duke.cs.osprey.ematrix;

import edu.duke.cs.osprey.confspace.*;
import edu.duke.cs.osprey.tools.MathTools;

import java.util.List;

/**
 * Score corrector for quantities that are independent of assignments outside the tuple
 */
public class SimpleScoreCorrector extends ScoreCorrector<TupE> {
    public SimpleScoreCorrector(List<SimpleConfSpace.Position> positions, MathTools.Optimizer opt) {
        super(positions, opt);
    }

    @Override
    protected double correctionSize(TupE correction) {
        // By default this value represents CorrectedScore - UnCorrectedScore
        return correction.E;
    }

}
