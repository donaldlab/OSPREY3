/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.confspace;

import java.util.ArrayList;

/**
 *
 * @author mhall44
 */
public class RCTuple {
    
    //a tuple of RCs
    public ArrayList<Integer> pos;//which flexible positions
    public ArrayList<Integer> RCs;//the RCs themselves (residue-specific numbering, as in the TupleMatrices)

    
    
    public RCTuple(ArrayList<Integer> pos, ArrayList<Integer> RCs) {
        this.pos = pos;
        this.RCs = RCs;
    }
    
    
    //one-RC tuple
    public RCTuple(int pos1, int RC1){
        pos = new ArrayList<>();
        RCs = new ArrayList<>();
        pos.add(pos1);
        RCs.add(RC1);
    }
    
    
    //a pair
    public RCTuple(int pos1, int RC1, int pos2, int RC2){
        pos = new ArrayList<>();
        RCs = new ArrayList<>();
        
        pos.add(pos1);
        pos.add(pos2);
        
        RCs.add(RC1);
        RCs.add(RC2);
    }
    
    //Sometimes we'll want to generate an RC tuple from a conformation, specified as RCs for all positions
    //in order.  
    //In this case, negative values are not (fully) defined, so the tuple contains all positions
    //with positive values in conf
    public RCTuple(int[] conf){
        
        pos = new ArrayList<>();
        RCs = new ArrayList<>();
        
        for(int posNum=0; posNum<conf.length; posNum++){
            if(conf[posNum]>=0){//RC fully defined
                pos.add(posNum);
                RCs.add(conf[posNum]);
            }
        }
    }
    
    
    public boolean isSameTuple(RCTuple tuple2){
        //do the two tuple objects specify the same tuple of RCs?
        if( (pos.size()!=RCs.size()) || (tuple2.pos.size()!=tuple2.RCs.size()) )
            throw new RuntimeException("ERROR: Ill-defined RC tuple");
        
        
        if(pos.size()!=tuple2.pos.size())
            return false;
        
        //tuples are well-defined and same size...check position by position
        for(int index=0; index<pos.size(); index++){
            if(pos.get(index)!=tuple2.pos.get(index))
                return false;
            if(RCs.get(index)!=tuple2.RCs.get(index))
                return false;
        }
        
        //if we get here they're the same
        return true;
    }
    
    
}
