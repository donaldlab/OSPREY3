/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.tupexp;

import edu.duke.cs.osprey.confspace.HigherTupleFinder;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SuperRCTuple;
import edu.duke.cs.osprey.confspace.SearchProblem;
import java.util.ArrayList;

/**
 * Implementation of a tuple expander that expands conformational energies--
 * including continuous flexibility, non-pairwise terms, etc.
 *
 * @author mhall44
 */
public class ConfETupleExpander extends TupleExpander {

    SearchProblem sp;

    public ConfETupleExpander(SearchProblem sp) {
        super(sp.confSpace.numPos, sp.confSpace.getNumRCsAtPos(), sp.pruneMat.getPruningInterval());
        this.sp = sp;
    }

    @Override
    double scoreAssignmentList(int[] assignmentList) {
        if (sp.useEPIC) {//Faster if we can score by EPIC

            double E = sp.EPICMinimizedEnergy(assignmentList);

            if (E == Double.POSITIVE_INFINITY) {//this is going to be a problem if used as a true value
                RCTuple tup = new RCTuple(assignmentList);
                if (isPruned(tup)) {
                    throw new RuntimeException("ERROR: Scoring pruned conformation: " + tup.stringListing());
                } else {
                    throw new RuntimeException("ERROR: Infinite E for unpruned conf: " + tup.stringListing());
                }
            }

            return E;
        } else//Score by minimizing energy function directly
        {
            return sp.minimizedEnergy(assignmentList);
        }
    }

    @Override
    boolean isPruned(RCTuple tup) {
        return sp.pruneMat.isPruned(tup);
    }

    @Override
    void pruneTuple(RCTuple tup) {
        sp.pruneMat.markAsPruned(tup);
    }

    @Override
    ArrayList<RCTuple> higherOrderPrunedTuples(RCTuple tup) {
        //list higher-order pruned tuples containing the pair tup

        if (tup.pos.size() != 2) {
            throw new RuntimeException("ERROR: higherOrderPrunedTuples is meant to take an RC pair as argument");
        }

        HigherTupleFinder<Boolean> htf
                = sp.pruneMat.getHigherOrderTerms(tup.pos.get(0), tup.RCs.get(0), tup.pos.get(1), tup.RCs.get(1));

        if (htf != null) {
            ArrayList<RCTuple> otherTups = htf.listInteractionsWithValue(true);
            //otherTups are recorded as what tup interacts with...add in tup to get the whole pruned tuple
            for (RCTuple otherTup : otherTups) {
                for (int i = 0; i < 2; i++) {
                    otherTup.pos.add(tup.pos.get(i));
                    otherTup.RCs.add(tup.RCs.get(i));
                }
            }

            return otherTups;
        } else//no interactions
        {
            return new ArrayList<>();
        }
    }

    ArrayList<SuperRCTuple> higherOrderPrunedTuples(SuperRCTuple tup) {
        //list higher-order pruned tuples containing the pair tup

        if (tup.pos.size() != 2) {
            throw new RuntimeException("ERROR: higherOrderPrunedTuples is meant to take an RC pair as argument");
        }

        HigherTupleFinder<Boolean> htf
                = sp.pruneMat.getHigherOrderTerms(tup.pos.get(0), tup.superRCs.get(0), tup.pos.get(1), tup.superRCs.get(1));

        if (htf != null) {
            ArrayList<SuperRCTuple> otherTups = htf.listInteractionsWithValueSuper(true);
            //otherTups are recorded as what tup interacts with...add in tup to get the whole pruned tuple
            for (SuperRCTuple otherTup : otherTups) {
                for (int i = 0; i < 2; i++) {
                    otherTup.pos.add(tup.pos.get(i));
                    otherTup.superRCs.add(tup.superRCs.get(i));
                }
            }
            return otherTups;
        } else//no interactions
        {
            return new ArrayList<>();
        }
    }
}
