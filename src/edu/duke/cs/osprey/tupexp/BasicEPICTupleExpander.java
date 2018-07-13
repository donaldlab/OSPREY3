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

import edu.duke.cs.osprey.confspace.HigherTupleFinder;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.epic.EPICMatrix;
import edu.duke.cs.osprey.ematrix.epic.NewEPICMatrix;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import java.util.ArrayList;

/**
 *
 * @author mhall44
 */
public class BasicEPICTupleExpander extends TupleExpander {
    
    NewEPICMatrix epicMat;
    PruningMatrix pruneMat;
    
    public BasicEPICTupleExpander(SimpleConfSpace confSpace, double pruningInterval,
            LUTESettings luteSettings, NewEPICMatrix epicMat, PruningMatrix pruneMat){
        super(confSpace.getNumPos(), confSpace.getNumResConfsByPos(), pruningInterval, luteSettings);
        this.epicMat = epicMat;
        this.pruneMat = pruneMat;
    }

    
    @Override
    double scoreAssignmentList(int[] assignmentList) {
        RCTuple tup = new RCTuple(assignmentList);
        double E = epicMat.minimizeEnergy(tup, true);
        
        if(E==Double.POSITIVE_INFINITY){//this is going to be a problem if used as a true value
            if(isPruned(tup))
                throw new RuntimeException("ERROR: Scoring pruned conformation: "+tup.stringListing());
            else
                throw new RuntimeException("ERROR: Infinite E for unpruned conf: "+tup.stringListing());
        }
        
        return E;
    }
          

    
    @Override
    boolean isPruned(RCTuple tup) {
        return pruneMat.isPruned(tup);
    }

    
    @Override
    void pruneTuple(RCTuple tup) {
        pruneMat.markAsPruned(tup);
    }

    
    @Override
    ArrayList<RCTuple> higherOrderPrunedTuples(RCTuple tup) {
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
