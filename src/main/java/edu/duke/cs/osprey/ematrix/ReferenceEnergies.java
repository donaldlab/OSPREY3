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

package edu.duke.cs.osprey.ematrix;

import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.confspace.RC;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.TreeMap;

/**
 *
 * This class computes a reference energy for each residue position & AA type
 * (which is the minimum intra-RC energy at that position & AA type)
 * and then corrects the energy matrix accordingly
 * 
 * @author mhall44
 */
public class ReferenceEnergies implements Serializable {
    
    
    ConfSpace cSpace;
    ArrayList<TreeMap<String,Double>> eRefMatrix = new ArrayList<>();//list (by position) of AA type --> Eref maps
    
    
    public ReferenceEnergies(ConfSpace cs){
        cSpace = cs;
        computeERef();
    }
    
    
    private void computeERef(){
        
        for(int pos=0; pos<cSpace.numPos; pos++){
            
            System.out.println("Starting intra energy calculations for position "+pos);
            
            //no shell residues need to be included...intra only
            TermECalculator intraECalc = new TermECalculator(cSpace, new ArrayList<>(), false,
                    true, null, null, false, pos);
            
            ArrayList<Double> intraEList = (ArrayList<Double>) intraECalc.doCalculation();//intraE for each RC at this pos
            
            //OK now build the Eref mapping
            TreeMap<String,Double> posERef = new TreeMap<String,Double>();
            
            ArrayList<RC> RCList = cSpace.posFlex.get(pos).RCs;
            
            for(int rcNum=0; rcNum<RCList.size(); rcNum++){
                
                String AAType = RCList.get(rcNum).AAType;
                double intraE = intraEList.get(rcNum);
                
                if(posERef.containsKey(AAType)){
                    posERef.put( AAType, Math.min(intraE,posERef.get(AAType)) );
                }
                else
                    posERef.put(AAType, intraE);
            }
            
            eRefMatrix.add(posERef);
        }
    }
    
    
    void correctEnergyMatrix(EnergyMatrix emat){
        //subtract the right eRef off of each intra+shell energy
        for(int pos=0; pos<cSpace.numPos; pos++){
            
            ArrayList<RC> RCList = cSpace.posFlex.get(pos).RCs;
            for(int rcNum=0; rcNum<RCList.size(); rcNum++){
                String AAType = RCList.get(rcNum).AAType;
                double eRef = eRefMatrix.get(pos).get(AAType);
                
                if(eRef == Double.POSITIVE_INFINITY){
                    //if all RCs for a residue type have infinite one-body energy 
                    //(i.e., are impossible),
                    //then they stay at infinity after eRef correction
                    eRef = 0;
                }
                
                double Euncorr = emat.getOneBody(pos, rcNum);
                double Ecorr = Euncorr - eRef;
                emat.setOneBody(pos, rcNum, Ecorr);
            }
        }
    }
    
    
    public double confERef(int[] conf){
        //Given RC's for each position, add up the total reference energy
        double totERef = 0;
        
        for(int pos=0; pos<cSpace.numPos; pos++){
        	totERef += posERef(pos, conf[pos]);
        }
        
        return totERef;
    }
    
    public double posERef(int pos, int rc) {
        return eRefMatrix.get(pos).get(cSpace.posFlex.get(pos).RCs.get(rc).AAType);
    }
}
