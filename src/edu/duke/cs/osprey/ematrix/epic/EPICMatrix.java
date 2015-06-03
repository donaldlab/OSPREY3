/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.ematrix.epic;

import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.TupleMatrix;
import edu.duke.cs.osprey.minimization.CCDMinimizer;
import edu.duke.cs.osprey.minimization.Minimizer;
import edu.duke.cs.osprey.minimization.MolecEObjFunction;
import java.util.ArrayList;

/**
 *
 * A matrix of 1-body, pairwise, etc. energies as polynomials in internal coordinates (EPIC)
 * Together with the EnergyMatrix, this allows precomputation of the complete energy of a continuously flexible system.  
 * 
 * @author mhall44
 */
public class EPICMatrix extends TupleMatrix<EPoly> {
    
    ConfSpace confSpace = null;//the conformational space for which the energy polynomials are defined
    //we need this to be able to get conformational energies from the polynomials

    
    public EPICMatrix(ConfSpace confSpace, double pruningInterval) {
        super(confSpace.numPos, confSpace.getNumRCsAtPos(),pruningInterval);
        this.confSpace = confSpace;
    }
    
    
    public double minContE(int[] conf){
        //Given a list of RCs for each position, compute the minimized continuous component of the energy
        //(i.e., conformationally minimized sum of EPIC term values, without minE)
        //negative RC numbers indicate undefined positions (exclude these positions: EPIC terms always positive, so this gives
        //a valid lower bound)
        //return infinity if the conformation is pruned
                
        RCTuple RCTup = new RCTuple(conf);
        
        EPICEnergyFunction efunc = internalEnergyFunction(RCTup);
        
        if(efunc==null)//pruned!
            return Double.POSITIVE_INFINITY;
        
        MolecEObjFunction objFcn = new MolecEObjFunction(efunc,confSpace,RCTup);
        
        Minimizer minim = new CCDMinimizer(objFcn,false);
        DoubleMatrix1D bestDOFVals = minim.minimize();
        double E = objFcn.getValue(bestDOFVals);
        
        return E;
    }
    
    
    public EPICEnergyFunction internalEnergyFunction(RCTuple tup){
        //Make an energy function representing the internal energy of the tuple
        //kind of an EPIC analog of EnergyMatrix.getInternalEnergy
        //return null if the tuple is pruned
        
        int numPosInTuple = tup.pos.size();
        
        ArrayList<EPoly> terms = new ArrayList<>();
        
        for(int indexInTuple=0; indexInTuple<numPosInTuple; indexInTuple++){
            int posNum = tup.pos.get(indexInTuple);
            int RCNum = tup.RCs.get(indexInTuple);
            
            EPoly intraE = getOneBody(posNum,RCNum);
            if(intraE==null)//pruned
                return null;
            
            terms.add( intraE );
            
            for(int index2=0; index2<indexInTuple; index2++){
                int pos2 = tup.pos.get(index2);
                int rc2 = tup.RCs.get(index2);
                
                EPoly pairwiseE = getPairwise(posNum,RCNum,pos2,rc2);
                if(pairwiseE==null)//pruned
                    return null;
                
                terms.add( pairwiseE );
                
                //higher-order energies will go here
            }
        }
        
        EPICEnergyFunction efunc = new EPICEnergyFunction(terms);
        return efunc;
    }
    
    
    
    
}
