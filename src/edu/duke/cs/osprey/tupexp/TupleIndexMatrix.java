/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.tupexp;

import java.util.ArrayList;

import edu.duke.cs.osprey.confspace.HigherTupleFinder;
import edu.duke.cs.osprey.confspace.TupleMatrixGeneric;

/**
 *
 * This TupleMatrix gives the indices for all the tuples in a TupleExpander
 * It is intended for fast lookup of what tuples are in a given sample
 * (for use in fitting: in large designs, storing all tuples for each samples
 * is taking too much memory)
 * 
 * @author mhall44
 */
public class TupleIndexMatrix extends TupleMatrixGeneric<Integer> {
    
    //default is -1 (tuple absent from expansion)
    
    public TupleIndexMatrix(int numPos, int[] numRCsAtPos, double pruningInterval){
        super(numPos, numRCsAtPos, pruningInterval, -1);
    }
    
    
    ArrayList<Integer> calcSampleTuples(int[] sample){
        //Convert list of samples (RC indices at each position) to list of its tuples in this matrix
        ArrayList<Integer> ans = new ArrayList<>();
        int numPos = sample.length;
        
        for(int pos=0; pos<numPos; pos++){
            int RCNum = sample[pos];
            
            int oneBodyIndex = getOneBody(pos,RCNum);
            if(oneBodyIndex != -1)
                ans.add(oneBodyIndex);
            
            for(int pos2=0; pos2<pos; pos2++){
                int rc2 = sample[pos2];
                
               int pairwiseIndex = getPairwise(pos,RCNum,pos2,rc2);
               if(pairwiseIndex != -1)
                   ans.add(pairwiseIndex);
                
               //now look for higher-order interactions of (RCNum,rc2) with lower-numbered residues (<pos2)
                HigherTupleFinder<Integer> htf = getHigherOrderTerms(pos,RCNum,pos2,rc2);
                if(htf != null)
                    calcHigherOrderTuples(ans,sample,pos2,htf);
            }
        }
        
        return ans;
    }
    
    
    void calcHigherOrderTuples(ArrayList<Integer> sampleTuples, int sample[], int maxPos, HigherTupleFinder<Integer> htf){
        //Look for higher-order tuples that are interactions of htf with positions <maxPos in sample.
        //Add these tuple's indices to sampleTuples
        ArrayList<Integer> interactingPos = htf.getInteractingPos();
        
        for(int ipos : interactingPos){
            
            if(ipos<maxPos){
                int iTupIndex = htf.getInteraction(ipos, sample[ipos]);
                
                if( iTupIndex != -1)//tuple found
                    sampleTuples.add(iTupIndex);
                
                //see if we need to go to even higher order
                HigherTupleFinder htf2 = htf.getHigherInteractions(ipos, sample[ipos]);
                if(htf2!=null){
                    calcHigherOrderTuples(sampleTuples, sample, ipos, htf2);
                }
            }
        }
    }
    
}