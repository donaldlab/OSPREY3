/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.ematrix;

import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.confspace.HigherTupleFinder;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.TupleMatrix;
import edu.duke.cs.osprey.pruning.PruningMatrix;

import java.util.ArrayList;

/**
 *
 * @author mhall44
 */
public class EnergyMatrix extends TupleMatrix<Double> {
    
    private double constTerm = 0;
    
    //we may want to have reference energies associated with this matrix
    ReferenceEnergies eRefMat = null;
    
    
    public EnergyMatrix(TupleMatrix<Double> tupleMatrix) {
    	super(tupleMatrix.higherTerms, tupleMatrix.pairwise, tupleMatrix.oneBody);
    }
    
    public EnergyMatrix(ConfSpace cSpace, double pruningInterval){
        //For a normal, precomputed energy matrix we expect infinite pruning interval
        //(matrix valid for all RCs),
        //but EnergyMatrix objects from tup-exp may only be valid for a finite pruning interval
        super(cSpace, pruningInterval, 0.);//default higher interaction for energies is 0
    }
    
    public EnergyMatrix(int numPos, int[] numRCsAtPos, double pruningInterval){
        super(numPos, numRCsAtPos, pruningInterval, 0.);
    }
    
    public double confE(int conf[]){
        //value of energy represented in energy matrix, for specified conformation
        //expressed as residue-specific RC indices (as in the storage matrices)
        return getInternalEnergy(new RCTuple(conf)) + constTerm;
    }
    
    
    public double getInternalEnergy(RCTuple tup){
        //internal energy of a tuple of residues when they're in the specified RCs
        
        int numPosInTuple = tup.pos.size();
        double E = 0;
        
        for(int indexInTuple=0; indexInTuple<numPosInTuple; indexInTuple++){
            int posNum = tup.pos.get(indexInTuple);
            int RCNum = tup.RCs.get(indexInTuple);
            
            double intraE = getOneBody(posNum,RCNum);
            E += intraE;
            
            for(int index2=0; index2<indexInTuple; index2++){
                int pos2 = tup.pos.get(index2);
                int rc2 = tup.RCs.get(index2);
                
                double pairwiseE = getPairwise(posNum,RCNum,pos2,rc2);
                E += pairwiseE;
                
                HigherTupleFinder<Double> htf = getHigherOrderTerms(posNum,RCNum,pos2,rc2);
                if(htf != null)
                    E += internalEHigherOrder(tup,index2,htf);
            }
        }
        
        return E;
    }
    
    
    double internalEHigherOrder(RCTuple tup, int curIndex, HigherTupleFinder<Double> htf){
        //Computes the portion of the internal energy for tuple tup
        //that consists of interactions in htf (corresponds to some sub-tuple of tup)
        //with RCs whose indices in tup are < curIndex
        double E = 0;
        ArrayList<Integer> interactingPos = htf.getInteractingPos();
        
        for(int ipos : interactingPos){
            
            //see if ipos is in tup with index < curIndex
            int iposIndex = -1;
            for(int ind=0; ind<curIndex; ind++){
                if(tup.pos.get(ind)==ipos){
                    iposIndex = ind;
                    break;
                }
            }

            if(iposIndex > -1){//ipos interactions need to be counted
                int iposRC = tup.RCs.get(iposIndex);
                E += htf.getInteraction(ipos, iposRC);
                
                //see if need to go up to highers order again...
                HigherTupleFinder htf2 = htf.getHigherInteractions(ipos,iposRC);
                if(htf2!=null){
                    E += internalEHigherOrder(tup,iposIndex,htf2);
                }
            }
        }
        
        return E;
    }
        
    /*
    public double getPairwiseE(int res1, int AA1, int rot1, int res2, int AA2, int rot2){
    //lookup by residue number, amino-acid index (for the residue's rotamer library),
    //and RC index for that residue and AA type (rotamer index for AA type if rigid backbone,
    //otherwise defined in the ConfSearchSpace)
    //RETURN ERROR IF RES1 AND RES2 ARE NOT SINGLE RESIDUES
    return getPairwise(res1, index1, res2, index2);
    }
     */
    //intra+shell similar...
    public double getConstTerm() {
        return constTerm;
    }

    public void setConstTerm(double constTerm) {
        this.constTerm = constTerm;
    }
    
    
    
    public double[][] topPairwiseInteractions(){
        //return the top absolute values of the pairwise interactions
        //between all pairs of positions
        int numPos = oneBody.size();
        
        double strongestPairE[][] = new double[numPos][numPos];
        
        for(int pos=0; pos<numPos; pos++){
            for(int pos2=0; pos2<pos; pos2++){
                for(int rc=0; rc<numRCsAtPos(pos); rc++){
                    for(int rc2=0; rc2<numRCsAtPos(pos2); rc2++){
                        strongestPairE[pos][pos2] = Math.max( strongestPairE[pos][pos2], Math.abs(getPairwise(pos, rc, pos2, rc2)) );
                        strongestPairE[pos2][pos] = strongestPairE[pos][pos2];
                    }
                }
            }
        }
        
        return strongestPairE;
    }

    public ReferenceEnergies geteRefMat() {
        return eRefMat;
    }

    
    public EnergyMatrix singleSeqMatrix(ArrayList<String> seq, ConfSpace origConfSpace){
    	// return (EnergyMatrix) super.singleSeqMatrix(seq, origConfSpace);
    	TupleMatrix<Double> tmat = super.singleSeqMatrix(seq, origConfSpace);
    	return new EnergyMatrix(tmat);
    }
    
    
	public EnergyMatrix singleSeqMatrix(ArrayList<String> seq, ArrayList<Integer> flexPos, ConfSpace origConfSpace){
		//return (PruningMatrix) super.singleSeqMatrix(seq, origConfSpace);
		TupleMatrix<Double> tmat = super.singleSeqMatrix(seq, flexPos, origConfSpace);
		return new EnergyMatrix(tmat);
	}
   
}
