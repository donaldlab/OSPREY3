/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.confspace;

import java.io.Serializable;
import java.util.ArrayList;
/**
 *
 * @author hmn5
 */
public class SuperRCTuple extends RCTuple{
    public ArrayList<Integer> superRCs;//sthe RCs themselves (residue-specific numbering, as in the TupleMatrices)
    
    public SuperRCTuple(){
        super();
        pos = new ArrayList<>();
        superRCs = new ArrayList<>();
    }
    
    public SuperRCTuple(ArrayList<Integer> pos, ArrayList<Integer> superRCs){
        this.pos = pos;
        this.superRCs = superRCs;
        this.RCs = superRCs;//we can treate superRCs as RCs in PruningMatrix so 
                            //we don't need to reqrite PruningMatrix
    }
    
    //one-RC
    public SuperRCTuple(int pos1, int superRC1){
        this.pos = new ArrayList<>();
        this.superRCs = new ArrayList<>();
        this.pos.add(pos1);
        this.superRCs.add(superRC1);
        this.RCs = this.superRCs;
    }
    
    //a pair
    public SuperRCTuple(int pos1, int superRC1, int pos2, int superRC2){
        this.pos = new ArrayList<>();
        this.superRCs = new ArrayList<>();
        
        this.pos.add(pos1);
        this.pos.add(pos2);
        
        this.superRCs.add(superRC1);
        this.superRCs.add(superRC2);
        
        this.RCs = this.superRCs;
    }
    //When we want to generate an RC tuple from a conformation, it will specify
    //super RCs at each position
    public SuperRCTuple(int[] conf){
        
        this.pos = new ArrayList<>();
        this.superRCs = new ArrayList<>();
        
        for(int posNum=0; posNum<conf.length; posNum++){
            if(conf[posNum]>=0){//RC fully defined
                this.pos.add(posNum);
                this.superRCs.add(conf[posNum]);
            }
        }        
    }

     public boolean isSameTuple(SuperRCTuple tuple2){
        //do the two tuple objects specify the same tuple of RCs?
        if( (pos.size()!= this.superRCs.size()) || (tuple2.pos.size()!=tuple2.superRCs.size()) )
            throw new RuntimeException("ERROR: Ill-defined RC tuple");
        
        
        if(pos.size()!=tuple2.pos.size())
            return false;
        
        //tuples are well-defined and same size...check position by position
        for(int index=0; index<pos.size(); index++){
            if(pos.get(index)!=tuple2.pos.get(index))
                return false;
            if(superRCs.get(index)!=tuple2.superRCs.get(index))
                return false;
        }
        
        //if we get here they're the same
        return true;
     }
     
   public String stringListing(){
        //Listing the RCs in a string
        String ans = "";
        
        for(int posNum=0; posNum<pos.size(); posNum++){
            ans = ans + "Pos " + pos.get(posNum) + " SuperRC " + superRCs.get(posNum) + " ";
        }
        
        return ans;
    }
    
}
