/*
** This file is part of OSPREY.
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation, either version 2 of the License, or
** (at your option) any later version.
** 
** OSPREY is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
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

