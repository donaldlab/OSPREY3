/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
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
