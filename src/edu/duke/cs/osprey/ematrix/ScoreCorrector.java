package edu.duke.cs.osprey.ematrix;

import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.confspace.*;
import edu.duke.cs.osprey.tools.MathTools;

import java.io.File;
import java.util.*;

/**
 * Class to store and retrieve corrections to energy scores
 *
 * @param <T>   extends TupE, the class that we use to map RCTuples to score corrections
 */
public abstract class ScoreCorrector <T extends RCTupleContainer> implements Correctable<T>{
    protected static boolean debug = false;

    protected TupleTrie<T> corrections;
    protected MathTools.Optimizer opt; // are we correcting up or down?
    protected int numPos;

    //Abstract methods that will depend on the type of corrector we are making
    protected abstract double correctionSize(T correction);

    public ScoreCorrector(List<SimpleConfSpace.Position> positions, MathTools.Optimizer opt) {
        this.corrections = new TupleTrie<>(positions);
        this.numPos = positions.size();
        this.opt = opt;
    }

    public double getCorrection(int[] assignments){
        return getCorrection(new RCTuple(assignments));
    }

    public double getCorrection(ConfIndex index){
        return getCorrection(index.makeConf());
    }

    public double getCorrection(RCTuple query){
        List<T> cover = getBestCorrectionsFor(query);
        double sum = 0.0;
        for (T correction : cover){
            //TODO: again, if we store these correction sizes per access this should be speedier
            sum += correctionSize(correction);
        }
        // We never want to loosen bounds
        if(opt.isBetter(sum, 0.0))
            return sum;
        else
            return 0.0;
    }

    public List<T> getAllCorrections(){
        return corrections.getAllEntries();
    }

    public int getNumCorrections(){
        return this.corrections.size();
    }

    protected List<T> getBestCorrectionsFor(RCTuple query){
        return largeWeightedSetPacking(corrections.getEntries(query));
    }

    /**
     * Returns a high-weight set packing of the list of corrections.
     *
     * NB: These corrections must by definition be independent of the unassigned conformation.
     *
     * Although ideally we would want this to be the *heaviest* possible set packing,
     * AFAIK that currently requires solving the maximum weighted set-packing problem,
     * which is in NP-hard.
     *
     * So, this uses a greedy strategy for now. This does not always return the optimal solution
     * in practice. Beware.
     *
     */
    protected List<T> largeWeightedSetPacking(List<T> confCorrections) {
        // Sort the corrections by correction size
        //TODO: This can involve computing hscores, so optimize this probably
        confCorrections.sort((a, b) -> opt.compare(correctionSize(a),correctionSize(b)));
        // Attempt 1: be greedy and start from the largest correction you
        // can get instead of trying to solve the NP-Complete problem.
        Set<Integer> usedPositions = new HashSet<>();
        List<T> usedCorrections = new ArrayList<>();
        for(T correction: confCorrections) {
            // If we've covered the set of residues already, break
            if (usedPositions.size() >= numPos) {
                break;
            }
            // If this correction overlaps one that we've used already,
            // then it's incompatible, so break and record that we have overlap
            Collection<Integer> positions = correction.tup.pos;
            boolean noIntersections = true;
            for(int position : positions) {
                if(usedPositions.contains(position)) {
                    noIntersections = false;
                    break;
                }
            }
            // if it doesn't overlap, then use it
            if(noIntersections) {
                usedPositions.addAll(correction.tup.pos);
                usedCorrections.add(correction);
            }
        }
        if(debug)
        {
            for(T correction: usedCorrections) {
                for(int pos: correction.tup.pos) {
                    if (usedPositions.contains(pos))
                        System.err.println("REUSING POSITION "+pos);
                }
            }
        }
        return usedCorrections;
    }

    public void insertCorrection(T correction){
        corrections.insert(correction);
    }

    public void insertAllCorrections(List<T> corrList){
        for (T correction : corrList){
            corrections.insert(correction);
        }
    }

    public boolean containsCorrectionFor(RCTuple tup){
        return corrections.contains(tup);
    }

    public void writeCorrectionsToFile(String filename){
        corrections.writeEntriesToFile(filename);
    }
    public void writeCorrectionsToFile(File f) {
        corrections.writeEntriesToFile(f);
    }

}
