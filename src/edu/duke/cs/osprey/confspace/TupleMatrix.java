/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.confspace;

import java.io.Serializable;
import java.util.ArrayList;

/**
 *
 * @author mhall44
 */
public class TupleMatrix<T> implements Serializable {
    //We will need "matrices" of quantities defined
    //for example, the energy matrix (T=Double) stores single, pairwise, and higher-order energies
    //and we'll also have pruning (T=Boolean) and EPIC (T=EPoly) matrices
    //we'll store things as ArrayLists to make it easier to merge residues, partition RCs, etc.
    //and also to facilitate generics
    
    //note: tuples are sets not ordered pairs, i.e. E(i_r,j_s) = E(j_s,i_r), and pruning (i_r,j_s) means pruning (j_s,i_r)
    
    public ArrayList<ArrayList<ArrayList<ArrayList<T>>>> pairwise;//pairwise energies.  Set up 4D to save space
    //indices: res1, res2, RC1, RC2 where res1>res2
    public ArrayList<ArrayList<T>> oneBody;//intra+shell
    
    public ArrayList<ArrayList<ArrayList<ArrayList<HigherTupleFinder<T>>>>> higherTerms;//look up higher terms by pair
    //same indices as pairwise
    //can be null if no interactions
    
    //maybe separate intra too?
    
    
    //The above all use RCs that are indexed by residues
    //the following lists indicate exactly what RCs those are
    //private because changing it may require re-indexing everything else
    //METHODS TO ADD RCs?
    
    //private ArrayList<ArrayList<RC>> RCList;
    
    /*
    //reverse lookup: for each residue, a lookup first by AA type, then by 
    //Return the residue-based RC number
    ArrayList<TreeMap<String,TreeMap<Integer,
    */
    
    
    double pruningInterval;//This matrix needs to hold entries for all RCs
    //that cannot be pruned with the specified pruning interval (Ew + Ival)
    //i.e. the matrix must describe all conformations within pruningInterval 
    //of the lowest pairwise lower bound
    
    T defaultHigherInteraction;//We only mark sparse higher interactions;
    //if unmarked we assume this value (e.g., 0 for energy, false for pruning)
    
    
    public TupleMatrix(ArrayList<ArrayList<ArrayList<ArrayList<HigherTupleFinder<T>>>>> higherTerms, 
    		ArrayList<ArrayList<ArrayList<ArrayList<T>>>> pairwise, 
    		ArrayList<ArrayList<T>> oneBody){
    	this.oneBody = oneBody;
    	this.pairwise = pairwise;
    	this.higherTerms = higherTerms;
    }
    
    
    public TupleMatrix(T defaultHigherInteraction){
        //no allocation (for overriding)
        this.defaultHigherInteraction = defaultHigherInteraction;
    }
    
    
    public TupleMatrix(ConfSpace cSpace, double pruningInterval, T defaultHigherInteraction){
        //allocate the matrix based on the provided conformational space
        init(cSpace.numPos, cSpace.getNumRCsAtPos(), pruningInterval, defaultHigherInteraction);
    }
    
    
    public TupleMatrix(int numPos, int[] numAllowedAtPos, double pruningInterval, T defaultHigherInteraction){
        //allocate the matrix based on the provided conformational space size
        //also specify what pruningInterval it's valid up to
        init(numPos, numAllowedAtPos, pruningInterval, defaultHigherInteraction);
    }
    
    
    private void init(int numPos, int[] numAllowedAtPos, double pruningInterval, T defaultHigherInteraction) {
        
        this.pruningInterval = pruningInterval;
        this.defaultHigherInteraction = defaultHigherInteraction;
        
        oneBody = new ArrayList<>();
        pairwise = new ArrayList<>();
        
        higherTerms = new ArrayList<>();//preallocate these too, but all null for now (no higher-order terms yet)
        
        for(int pos=0; pos<numPos; pos++){
                        
            int numRCs = numAllowedAtPos[pos];
            
            //preallocate oneBody for this position
            ArrayList<T> oneBodyAtPos = new ArrayList<>();
            for(int rc=0; rc<numRCs; rc++)//preallocate oneBody
                oneBodyAtPos.add(null);
            
            oneBodyAtPos.trimToSize();//we may need to save space so we'll trim everything to size
            oneBody.add(oneBodyAtPos);

            
            ArrayList<ArrayList<ArrayList<T>>> pairwiseAtPos = new ArrayList<>();
            //may want to leave some pairs of positions null if negligible interaction expected...
            //handle later though
            
            ArrayList<ArrayList<ArrayList<HigherTupleFinder<T>>>> higherOrderAtPos = new ArrayList<>();
            
            for(int pos2=0; pos2<pos; pos2++){
                
                int numRCs2 = numAllowedAtPos[pos2];

                ArrayList<ArrayList<T>> pairwiseAtPair = new ArrayList<>();
                ArrayList<ArrayList<HigherTupleFinder<T>>> higherOrderAtPair = new ArrayList<>();
                
                for(int rc=0; rc<numRCs; rc++){
                    ArrayList<T> pairwiseAtRC = new ArrayList<>();
                    ArrayList<HigherTupleFinder<T>> higherOrderAtRC = new ArrayList<>();
                    
                    for(int rc2=0; rc2<numRCs2; rc2++){
                        pairwiseAtRC.add(null);
                        higherOrderAtRC.add(null);
                    }
                    
                    pairwiseAtRC.trimToSize();
                    pairwiseAtPair.add(pairwiseAtRC);
                    
                    higherOrderAtRC.trimToSize();
                    higherOrderAtPair.add(higherOrderAtRC);
                }
                
                pairwiseAtPair.trimToSize();
                pairwiseAtPos.add(pairwiseAtPair);
                
                higherOrderAtPair.trimToSize();
                higherOrderAtPos.add(higherOrderAtPair);
            }
            
            pairwiseAtPos.trimToSize();
            pairwise.add(pairwiseAtPos);
            
            higherOrderAtPos.trimToSize();
            higherTerms.add(higherOrderAtPos);
        }
        
        oneBody.trimToSize();
        pairwise.trimToSize();
        higherTerms.trimToSize();
    }
    
    
    public T getPairwise(int res1, int index1, int res2, int index2){
        //working with residue-specific RC indices directly.  
        if(res1>res2)
            return pairwise.get(res1).get(res2).get(index1).get(index2);
        else
            return pairwise.get(res2).get(res1).get(index2).get(index1);
    }
    
    public T getOneBody(int res, int index){
        return oneBody.get(res).get(index);
    }
    
    
    public void setPairwise(int res1, int index1, int res2, int index2, T val){

        if(res1>res2)
            pairwise.get(res1).get(res2).get(index1).set(index2,val);
        else
            pairwise.get(res2).get(res1).get(index2).set(index1,val);
    }
    
    
    public void setOneBody(int res, int index, T val){
        oneBody.get(res).set(index,val);
    }
    
    
    public int numRCsAtPos(int pos){
        return oneBody.get(pos).size();
    }
    
    public int numPos(){
        return oneBody.size();
    }
    
    public void setTupleValue(RCTuple tup, T val){
        //assign the given value to the specified RC tuple
        int tupSize = tup.pos.size();
        
        if(tupSize==1)//just a one-body quantity
            setOneBody( tup.pos.get(0), tup.RCs.get(0), val);
        else if(tupSize==2)//two-body
            setPairwise( tup.pos.get(0), tup.RCs.get(0), tup.pos.get(1), tup.RCs.get(1), val );
        else if(tupSize>2){//higher-order
            setHigherOrder(tup,val);
        }
        else
            throw new UnsupportedOperationException( "ERROR: Not supporting tuple size " + tupSize );
    }
    
    
    public void setHigherOrder(RCTuple tup, T val){
        //set a higher-order term
        //we need all pairs contained in tup to know about it
        
        //loop over pairs
        for(int index1=0; index1<tup.pos.size(); index1++){
            for(int index2=0; index2<index1; index2++){
                
                int pos1 = tup.pos.get(index1);
                int rc1 = tup.RCs.get(index1);
                int pos2 = tup.pos.get(index2);
                int rc2 = tup.RCs.get(index2);
                
                //put tup into the HigherTupleFinder for this pair
                HigherTupleFinder<T> htf = getHigherOrderTerms(pos1,rc1,pos2,rc2);
                
                //create a HigherTupleFinder if there is none yet
                if(htf==null){
                    htf = new HigherTupleFinder(defaultHigherInteraction);
                    
                    if(pos1>pos2)
                        higherTerms.get(pos1).get(pos2).get(rc1).set(rc2,htf);
                    else
                        higherTerms.get(pos2).get(pos1).get(rc2).set(rc1,htf);
                }
                
                RCTuple subTup = tup.subtractMember(index1).subtractMember(index2);
                
                htf.setInteraction(subTup,val);
            }
        }
        
    }
    
    public double getPruningInterval() {
        return pruningInterval;
    }
    
    public void setPruningInterval(double pruningInterval) {
        this.pruningInterval = pruningInterval;
    }
    
    public HigherTupleFinder<T> getHigherOrderTerms(int res1, int index1, int res2, int index2){
        //working with residue-specific RC indices directly.  
        if(res1>res2)
            return higherTerms.get(res1).get(res2).get(index1).get(index2);
        else
            return higherTerms.get(res2).get(res1).get(index2).get(index1);
    }
    
    
    public TupleMatrix<T> singleSeqMatrix(ArrayList<String> seq, ConfSpace origConfSpace){
    	//cut down the matrix to only include RCs with the appropriate amino acid types for this sequence
    	//the original (all-sequences) conformational space is provided
    	int numPos = oneBody.size();
    	
    	ArrayList<ArrayList<T>> newOneBody = new ArrayList<>();
    	ArrayList<ArrayList<ArrayList<ArrayList<T>>>> newPairwise = new ArrayList<>();
    	ArrayList<ArrayList<ArrayList<ArrayList<HigherTupleFinder<T>>>>> newHigherTerms = new ArrayList<>();
    
    	
    	for(int pos=0; pos<numPos; pos++){
    		
    		String AAType = seq.get(pos);
    		ArrayList<T> newOneBodyAtPos = new ArrayList<>();
    		ArrayList<ArrayList<ArrayList<T>>> newPairwiseAtPos = new ArrayList<>();
    		ArrayList<ArrayList<ArrayList<HigherTupleFinder<T>>>> newHigherOrderAtPos = new ArrayList<>();
    		
    		for(int RCNum=0; RCNum<oneBody.get(pos).size(); RCNum++){
    			if(origConfSpace.posFlex.get(pos).RCs.get(RCNum).AAType.equalsIgnoreCase(AAType)){
    				newOneBodyAtPos.add(getOneBody(pos,RCNum));
    			}
    		}
    		
    		newOneBodyAtPos.trimToSize();
    		newOneBody.add(newOneBodyAtPos);
    		
    		for(int pos2=0; pos2<pos; pos2++){
    			String AAType2 = seq.get(pos2);
    			ArrayList<ArrayList<T>> newPairwiseAtPair = new ArrayList<>();
    			ArrayList<ArrayList<HigherTupleFinder<T>>> newHigherOrderAtPair = new ArrayList<>();
    			
    			for(int RCNum=0; RCNum<oneBody.get(pos).size(); RCNum++){
        			if(origConfSpace.posFlex.get(pos).RCs.get(RCNum).AAType.equalsIgnoreCase(AAType)){
        				
        				ArrayList<T> newPairwiseAtRC = new ArrayList<>();
        				ArrayList<HigherTupleFinder<T>> newHigherOrderAtRC = new ArrayList<>();
        				
        				for(int RCNum2=0; RCNum2<oneBody.get(pos2).size(); RCNum2++){
        					if(origConfSpace.posFlex.get(pos2).RCs.get(RCNum2).AAType.equalsIgnoreCase(AAType2)){
        						newPairwiseAtRC.add(getPairwise(pos,RCNum,pos2,RCNum2));
        						newHigherOrderAtRC.add(getHigherOrderTerms(pos,RCNum,pos2,RCNum2));
        					}
        				}
        				
        				newPairwiseAtRC.trimToSize();
        				newPairwiseAtPair.add(newPairwiseAtRC);
        				
        				newHigherOrderAtRC.trimToSize();
        				newHigherOrderAtPair.add(newHigherOrderAtRC);
        			}
        		}
    			
    			newPairwiseAtPair.trimToSize();
    			newPairwiseAtPos.add(newPairwiseAtPair);
    			
    			newHigherOrderAtPair.trimToSize();
    			newHigherOrderAtPos.add(newHigherOrderAtPair);
    		}
    		
    		newPairwise.add(newPairwiseAtPos);
    		newHigherTerms.add(newHigherOrderAtPos);
    	}
    	
    	newOneBody.trimToSize();
    	newPairwise.trimToSize();
    	newHigherTerms.trimToSize();

    	return new TupleMatrix<T>(newHigherTerms, newPairwise, newOneBody);
    }
    
    
    public TupleMatrix<T> singleSeqMatrix(ArrayList<String> seq, 
    		ArrayList<Integer> flexPos, ConfSpace origConfSpace){
    	//cut down the matrix to only include RCs with the appropriate amino acid types for this sequence
    	//the original (all-sequences) conformational space is provided
    	
    	ArrayList<ArrayList<T>> newOneBody = new ArrayList<>();
    	ArrayList<ArrayList<ArrayList<ArrayList<T>>>> newPairwise = new ArrayList<>();
    	ArrayList<ArrayList<ArrayList<ArrayList<HigherTupleFinder<T>>>>> newHigherTerms = new ArrayList<>();
    	
    	
    	for(int index = 0; index < seq.size(); ++index) {

    		int pos = flexPos.get(index);
    		
    		String AAType = seq.get(index);
    		ArrayList<T> newOneBodyAtPos = new ArrayList<>();
    		ArrayList<ArrayList<ArrayList<T>>> newPairwiseAtPos = new ArrayList<>();
    		ArrayList<ArrayList<ArrayList<HigherTupleFinder<T>>>> newHigherOrderAtPos = new ArrayList<>();
    		
    		for(int RCNum=0; RCNum<oneBody.get(pos).size(); RCNum++){
    			if(origConfSpace.posFlex.get(pos).RCs.get(RCNum).AAType.equalsIgnoreCase(AAType)){
    				newOneBodyAtPos.add(getOneBody(pos,RCNum));
    			}
    		}
    		
    		newOneBodyAtPos.trimToSize();
    		newOneBody.add(newOneBodyAtPos);
    		
    		for(int index2 = 0; index2 < index; ++index2) {
    			
    			int pos2 = flexPos.get(index2);
    			
    			String AAType2 = seq.get(index2);
    			ArrayList<ArrayList<T>> newPairwiseAtPair = new ArrayList<>();
    			ArrayList<ArrayList<HigherTupleFinder<T>>> newHigherOrderAtPair = new ArrayList<>();
    			
    			for(int RCNum=0; RCNum<oneBody.get(pos).size(); RCNum++){
        			if(origConfSpace.posFlex.get(pos).RCs.get(RCNum).AAType.equalsIgnoreCase(AAType)){
        				
        				ArrayList<T> newPairwiseAtRC = new ArrayList<>();
        				ArrayList<HigherTupleFinder<T>> newHigherOrderAtRC = new ArrayList<>();
        				
        				for(int RCNum2=0; RCNum2<oneBody.get(pos2).size(); RCNum2++){
        					if(origConfSpace.posFlex.get(pos2).RCs.get(RCNum2).AAType.equalsIgnoreCase(AAType2)){
        						newPairwiseAtRC.add(getPairwise(pos,RCNum,pos2,RCNum2));
        						newHigherOrderAtRC.add(getHigherOrderTerms(pos,RCNum,pos2,RCNum2));
        					}
        				}
        				
        				newPairwiseAtRC.trimToSize();
        				newPairwiseAtPair.add(newPairwiseAtRC);
        				
        				newHigherOrderAtRC.trimToSize();
        				newHigherOrderAtPair.add(newHigherOrderAtRC);
        			}
        		}
    			
    			newPairwiseAtPair.trimToSize();
    			newPairwiseAtPos.add(newPairwiseAtPair);
    			
    			newHigherOrderAtPair.trimToSize();
    			newHigherOrderAtPos.add(newHigherOrderAtPair);
    		}
    		
    		newPairwise.add(newPairwiseAtPos);
    		newHigherTerms.add(newHigherOrderAtPos);
    	}
    	
    	newOneBody.trimToSize();
    	newPairwise.trimToSize();
    	newHigherTerms.trimToSize();

    	return new TupleMatrix<T>(newHigherTerms, newPairwise, newOneBody);
    }
  
}
