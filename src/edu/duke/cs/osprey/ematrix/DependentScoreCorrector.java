package edu.duke.cs.osprey.ematrix;

import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.astar.conf.scoring.AStarScorer;
import edu.duke.cs.osprey.astar.conf.scoring.ScorerFactory;
import edu.duke.cs.osprey.confspace.*;
import edu.duke.cs.osprey.tools.MathTools;

import java.util.Arrays;
import java.util.List;

/**
 * Score corrector for quantities that are dependent on assignments outside the tuple, and
 * therefore need to score on a new energy matrix
 */
public class DependentScoreCorrector extends ScoreCorrector<MappableTupE> {
    protected EnergyMatrix mappedEmat;
    protected AStarScorer corrGScorer;
    protected AStarScorer corrHScorer;

    // Needed to calculate new scores
    protected ConfIndex index;  // this will not be a shared object, so no need to be polite
    protected RCs rcs;
    protected double unCorrectedScore;

    // Needed to be polite
    protected int[] oldAssignments;

    public DependentScoreCorrector(SimpleConfSpace confSpace, MathTools.Optimizer opt) {
        super(confSpace, opt);
        this.mappedEmat = null;
        this.corrGScorer = null;
        this.corrHScorer = null;

        this.index = null;
        this.rcs = null;
        this.unCorrectedScore = Double.NaN;

        this.oldAssignments = new int[confSpace.getNumPos()];
        Arrays.fill(this.oldAssignments, -1);
    }

    public DependentScoreCorrector(SimpleConfSpace confSpace, MathTools.Optimizer opt, EnergyMatrix mappedEmat, ScorerFactory gScorerFactory, ScorerFactory hScorerFactory){
        super(confSpace, opt);
        this.mappedEmat = mappedEmat;
        this.corrGScorer = gScorerFactory.make(mappedEmat);
        this.corrHScorer = hScorerFactory.make(mappedEmat);

        this.index = null;
        this.rcs = null;
        this.unCorrectedScore = Double.NaN;

        this.oldAssignments = new int[confSpace.getNumPos()];
        Arrays.fill(this.oldAssignments, -1);
    }

    @Override
    protected TupETrie<MappableTupE> makeTrie(List<SimpleConfSpace.Position> positions) {
        return new MappableTupETrie(positions);
    }

    @Override
    protected double correctionSize(MappableTupE correction) {
        correction.mappedTup.pasteToIndex(this.index);
        double correctScore = this.corrGScorer.calc(this.index, this.rcs) + this.corrHScorer.calc(this.index, this.rcs);
        return correctScore - this.unCorrectedScore;
    }

    public double getCorrection(RCTuple query, double unCorrectedScore, ConfIndex index, RCs rcs){
        // Set info for scorer, needed for the correctionSize method
        setRequiredInformation(unCorrectedScore, index, rcs);

        // go ahead and get the corrections
        List<MappableTupE> cover = getBestCorrectionsFor(query);
        double sum = 0.0;
        for (MappableTupE correction : cover){
            //TODO: again, if we store these correction sizes per access this should be speedier
            sum += correctionSize(correction);
        }

        // Unset info to get a clean slate for next correction
        unsetRequiredInformation();
        // We never want to loosen bounds
        if(opt.isBetter(sum, 0.0))
            return sum;
        else
            return 0.0;
    }

    public double getCorrection(int[] assignments, double unCorrectedScore, ConfIndex index, RCs rcs){
        return getCorrection(new RCTuple(assignments), unCorrectedScore, index, rcs);
    }

    protected void setRequiredInformation(double unCorrectedScore, ConfIndex index, RCs rcs){
        this.unCorrectedScore = unCorrectedScore;
        setConfIndex(index);
        this.rcs = rcs;

    }

    protected void unsetRequiredInformation(){
        this.unCorrectedScore = Double.NaN;
        this.index = null;
        this.rcs = null;
    }

    protected void setConfIndex(ConfIndex index){
        this.index = new ConfIndex(index);
    }

    public void setMappedEmat(EnergyMatrix emat){
        this.mappedEmat = emat;
    }

    public void setGScorerFactory(ScorerFactory fact){
        this.corrGScorer = fact.make(mappedEmat);
    }

    public void setHScorerFactory(ScorerFactory fact){
        this.corrHScorer = fact.make(mappedEmat);
    }
}
