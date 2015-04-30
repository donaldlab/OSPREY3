/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.ematrix;

import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.TupleMatrix;

/**
 *
 * @author mhall44
 */
public class EnergyMatrix extends TupleMatrix<Double> {
    
    public EnergyMatrix(ConfSpace cSpace){
        super(cSpace);
    }
    
    
    public double confE(int conf[]){
        //value of energy represented in energy matrix, for specified conformation
        //expressed as residue-specific RC indices (as in the storage matrices)
        return getInternalEnergy(new RCTuple(conf));
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
                
                //higher-order energies will go here
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
    
}
