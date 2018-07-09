/*
** This file is part of OSPREY 3.0
** 
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
** 
** OSPREY relies on grants for its development, and since visibility
** in the scientific literature is essential for our success, we
** ask that users of OSPREY cite our papers. See the CITING_OSPREY
** document in this distribution for more information.
** 
** Contact Info:
**    Bruce Donald
**    Duke University
**    Department of Computer Science
**    Levine Science Research Center (LSRC)
**    Durham
**    NC 27708-0129
**    USA
**    e-mail: www.cs.duke.edu/brd/
** 
** <signature of Bruce Donald>, Mar 1, 2018
** Bruce Donald, Professor of Computer Science
*/

package edu.duke.cs.osprey.tupexp;

import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.confspace.HigherTupleFinder;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.plug.LPChecks;
import edu.duke.cs.osprey.plug.PolytopeMatrix;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import java.util.ArrayList;
import org.apache.commons.math3.optim.linear.LinearConstraint;

/**
 *
 * @author mhall44
 */
public class BasicPruningTupleExpander extends TupleExpander {

    PruningMatrix pruneMat;
    PolytopeMatrix plugMat;//if null just don't use PLUG

    public BasicPruningTupleExpander(PruningMatrix pruneMat, PolytopeMatrix plugMat){
        super(pruneMat.getNumPos(), pruneMat.getNumConfAtPos(), pruneMat.getPruningInterval(), new LUTESettings());
        this.pruneMat = pruneMat;
        this.plugMat = plugMat;
        if(plugMat!=null)
            throw new RuntimeException("ERROR: Expected rejectInf to be on...");
    }

    @Override
    public double scoreAssignmentList(int[] assignmentList) {
        if(isPruned(new RCTuple(assignmentList)))
            return Double.POSITIVE_INFINITY;
        if(plugMat!=null){
            ArrayList<LinearConstraint> fullPolytope = plugMat.getFullPolytope(new RCTuple(assignmentList));
            if(fullPolytope==null)//clash
                return Double.POSITIVE_INFINITY;
            if( ! LPChecks.polytopeHasFeasiblePt(fullPolytope) )
                return Double.POSITIVE_INFINITY;
        }
        return 0;
    }

    @Override
    public boolean isPruned(RCTuple tup) {
        return pruneMat.isPruned(tup);
    }

    @Override
    public void pruneTuple(RCTuple tup) {
        pruneMat.markAsPruned(tup);
    }

    @Override
    public ArrayList<RCTuple> higherOrderPrunedTuples(RCTuple tup) {
        //Just like ConfETupleExpander
        //list higher-order pruned tuples containing the pair tup

        if(tup.pos.size() != 2)
            throw new RuntimeException("ERROR: higherOrderPrunedTuples is meant to take an RC pair as argument");

        HigherTupleFinder<Boolean> htf =
                pruneMat.getHigherOrderTerms(tup.pos.get(0), tup.RCs.get(0), tup.pos.get(1), tup.RCs.get(1));

        if(htf!=null){
            ArrayList<RCTuple> otherTups = htf.listInteractionsWithValue(true);
            //otherTups are recorded as what tup interacts with...add in tup to get the whole pruned tuple
            for(RCTuple otherTup: otherTups){
                for(int i=0; i<2; i++){
                    otherTup.pos.add(tup.pos.get(i));
                    otherTup.RCs.add(tup.RCs.get(i));
                }
            }

            return otherTups;
        }
        else//no interactions
            return new ArrayList<>();
    }

}
