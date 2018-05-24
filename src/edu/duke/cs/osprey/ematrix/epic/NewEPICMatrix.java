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

package edu.duke.cs.osprey.ematrix.epic;

import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.confspace.ParametricMolecule;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.TupleMatrixGeneric;
import edu.duke.cs.osprey.minimization.CCDMinimizer;
import edu.duke.cs.osprey.minimization.Minimizer;
import edu.duke.cs.osprey.minimization.MoleculeObjectiveFunction;
import java.util.ArrayList;

/**
 *
 * SimpleConfSpace/MoleculeObjectiveFunction version of EPICMatrix
 * also made to provide full EPIC energies instead of E-mtx lower bounds + EPIC cont part,
 * thus avoiding error if the minimizer found different minima when computing the E mtx and the EPIC mtx
 * 
 * @author mhall44
 */
public class NewEPICMatrix extends TupleMatrixGeneric<EPoly> {
    
    SimpleConfSpace confSpace = null;//the conformational space for which the energy polynomials are defined
    //we need this to be able to get conformational energies from the polynomials

    
    public NewEPICMatrix(SimpleConfSpace confSpace, double pruningInterval) {
        super(confSpace.getNumPos(), confSpace.getNumResConfsByPos(), pruningInterval, null);
        this.confSpace = confSpace;
    }
    
    
    
   public double minimizeEnergy(RCTuple RCTup, boolean includeMinE){
       //Get the minimized energy for the tuple (internal plus interactions of 
       //positions in the tuple with the shell), as approximated by EPIC
       //if not includeMinE then it's just the continuous portion
       //note this could introduce error if the minimizer found different minima
       //when making E and EPIC mtx's, so when possible (e.g., LUTE fitting) includeMinE is better 
       
       
            //137DEBUG!!!  seeing if min taking too long
            long startTime = System.currentTimeMillis();
       
       
       
        EPICEnergyFunction efunc = internalEnergyFunction(RCTup,includeMinE);
        
        //137DEBUG!!
        long gotEFcnTime = System.currentTimeMillis()-startTime;
        
        
        if(efunc==null)//pruned!
            return Double.POSITIVE_INFINITY;
        
        ParametricMolecule bpmol = confSpace.makeMolecule(RCTup);
        MoleculeObjectiveFunction objFcn = new MoleculeObjectiveFunction(bpmol, efunc);
        
        Minimizer minim = new CCDMinimizer(objFcn,false);
        
        
        //137DEBUG!!
        long preMinTime = System.currentTimeMillis()-startTime;
        
        
        
        DoubleMatrix1D bestDOFVals = minim.minimize().dofValues;
        double E = objFcn.getValue(bestDOFVals);
        
        //DEBUG!!!!!!!!
        /*String objFcnFile = "OBJFCN"+System.currentTimeMillis()+".dat";
        ObjectIO.writeObject(objFcn, objFcnFile);
        System.out.println("OUTPUTTED OBJ FCN TO "+objFcnFile+".  E: "+E+" bestDOFVals: "+bestDOFVals);
        */
        
        //137DEBUG!!!
        long sampTime = System.currentTimeMillis() - startTime;
        //taking over 10 s is going to be an issue
        if(sampTime > 10000){
            System.out.println();
            System.out.println("Minimization took over 10 s (ms shown): "+sampTime);
            System.out.println("Time to get E Fcn: "+gotEFcnTime+" to start min: "+preMinTime);
            System.out.println("Sample: "+RCTup.stringListing());
            System.out.println("Energy: "+E);
            CCDMinimizer ccdMin = (CCDMinimizer)minim;
            
            System.out.println("GVCounts.  Estmin: "+ccdMin.GVCountEstmin+" Edge: "+ccdMin.GVCountEdge
                    +" Bigger: "+ccdMin.GVCountBigger+" Smaller: "+ccdMin.GVCountSmaller);
            
            System.out.println();
        }
        
        efunc.unassignSharedMolec();
        
        return E;
    }
   
    public double minimizeEnergy(int[] conf){
        ////Given a list of RCs for each position, compute the minimized energy
        //negative RC numbers indicate undefined positions 
        //(exclude these positions: EPIC terms always positive, so this gives
        //a valid lower bound).  again, return infinity if conf pruned
        return minimizeEnergy(new RCTuple(conf), true);
    }
    
    
    public double minContE(int[] conf){
        return minimizeEnergy(new RCTuple(conf), false);
    }
        
    
    public EPICEnergyFunction internalEnergyFunction(RCTuple tup, boolean includeMinE){
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
                //though currently not supporting these for EPIC
                if( getHigherOrderTerms(posNum,RCNum,pos2,rc2) != null ){
                    throw new UnsupportedOperationException("ERROR: Higher-order EPIC terms "
                            + "not yet supported in internalEnergyFunction");
                }
            }
        }
        
        EPICEnergyFunction efunc = new EPICEnergyFunction(terms, includeMinE);
        return efunc;
    }

    public SimpleConfSpace getConfSpace() {
        return confSpace;
    }
    
    
    
}
