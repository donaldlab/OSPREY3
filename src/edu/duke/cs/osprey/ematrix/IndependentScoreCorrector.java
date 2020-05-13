package edu.duke.cs.osprey.ematrix;

import edu.duke.cs.osprey.confspace.*;
import edu.duke.cs.osprey.tools.MathTools;

import java.util.List;

/**
 * Score corrector for quantities that are independent of assignments outside the tuple
 */
public class IndependentScoreCorrector extends ScoreCorrector<TupE> {
    public IndependentScoreCorrector(SimpleConfSpace confSpace, MathTools.Optimizer opt) {
        super(confSpace, opt);
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

    public double getCorrection(int[] assignments){
        return getCorrection(new RCTuple(assignments));
    }

    public double getCorrection(RCTuple query){
        List<TupE> cover = getBestCorrectionsFor(query);
        double sum = 0.0;
        for (TupE correction : cover){
            //TODO: again, if we store these correction sizes per access this should be speedier
            sum += correctionSize(correction);
        }
        // We never want to loosen bounds
        if(opt.isBetter(sum, 0.0))
            return sum;
        else
            return 0.0;
    }
}
