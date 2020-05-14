package edu.duke.cs.osprey.ematrix;

import edu.duke.cs.osprey.confspace.*;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.sharkstar.SHARKStarNodeScorer;
import edu.duke.cs.osprey.tools.MathTools;

import java.util.List;

/**
 * Used to compute on-the-fly extensions to an energy matrix
 *
 */
public class MappingEnergyMatrix extends ProxyEnergyMatrix implements Correctable<TupMapping>{

    protected TupleTrie<TupMapping> tupTrie;
    protected EnergyMatrix mappedEmat;

    public MappingEnergyMatrix(SimpleConfSpace confSpace, EnergyMatrix target) {
        super(confSpace, target);
    }

    @Override
    public void insertCorrection(TupMapping tupMapping) {
        tupTrie.insert(tupMapping);
    }

    @Override
    public Double getOneBody(int pos, int rc) {
        return super.getOneBody(pos, rc);
    }

    @Override
    public double getEnergy(int pos1, int rc1, int pos2, int rc2) {
        return super.getEnergy(pos1, rc1, pos2, rc2);
    }

    @Override
    public Double getPairwise(int pos1, int rc1, int pos2, int rc2) {
        return super.getPairwise(pos1, rc1, pos2, rc2);
    }

    @Override
    public double getEnergy(int pos, int rc) {
        return super.getEnergy(pos, rc);
    }

    @Override
    public boolean containsCorrectionFor(RCTuple tup) {
        return tupTrie.contains(tup);
    }
}

