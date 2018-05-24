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

package edu.duke.cs.osprey.dof.deeper;

import edu.duke.cs.osprey.restypes.HardCodedResidueInfo;
import edu.duke.cs.osprey.structure.Atom;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.RigidBodyMotion;
import java.io.Serializable;
import java.util.HashMap;
import java.util.Set;

/**
 *
 * Records backbone coordinates for a residue,
 * so we can restore them when performing perturbations
 * Sidechain is moved as a rigid body (including CA)
 * 
 * @author mhall44
 */
public class ResBBState implements Serializable {
    
    HashMap<String,double[]> coords = new HashMap<>();
    //Residue res;
    
    double CBCoord[] = null;//do we know where to put the sidechain
    
    final static double HNProRatio = 0.692297285699751;
    //idealized ratio of N-H to Pro N-CD bond lengths
    
    
    
    public ResBBState(Residue res){
        //Record the state of res
                
        for(Atom at : res.atoms){
            boolean isBBAtom = false;
            for(String BBAtomName : HardCodedResidueInfo.possibleBBAtoms){
                if(at.name.equalsIgnoreCase(BBAtomName)){
                    isBBAtom = true;
                    break;
                }
            }
            
            if(isBBAtom)
                coords.put(at.name, at.getCoords());
            
            else if(at.name.equalsIgnoreCase("CB")){
                CBCoord = at.getCoords();
            }
            else if(res.fullName.startsWith("GLY")){
                if(at.name.equalsIgnoreCase("HA3"))//CB-like HA
                    CBCoord = at.getCoords();//since we're just using it to get a rotation
                    //about CA, it's OK if it's a little closer than usual to CA
            }
            
            
            if(res.fullName.startsWith("PRO")){
                if(at.name.equalsIgnoreCase("CD")){
                    double HCoord[] = at.getCoords();
                    //convert to H so can use if there's a mutation
                    rescaleBondLen(HCoord, res.getCoordsByAtomName("N"), HNProRatio);
                    coords.put("H", HCoord);
                }
            }
        }
    }
    
    
    public ResBBState(ResBBState state2){//deep copy
        if(state2.CBCoord!=null)
            CBCoord = state2.CBCoord.clone();
        
        for(String s : state2.coords.keySet())
            coords.put(s, state2.coords.get(s).clone());
    }
    
    
    public void putInState(Residue res){
        //Put res in the conformational state defined by the backbone coordinates recorded here
        
        Set<String> BBAtomNames = coords.keySet();
        
        //If there is a CB, move the sidechain and HA as a rigid body including CA
        //gly HA's can be placed exactly by sidechain idealization, so don't worry about them
        if(CBCoord!=null){
            double[] resCBCoord = res.getCoordsByAtomName("CB");
            
            if(resCBCoord!=null){
                
                double resCACoord[] = res.getCoordsByAtomName("CA");
                
                RigidBodyMotion sidechainMotion = new RigidBodyMotion(
                        new double[][] { resCACoord, resCBCoord, new double[3] },
                        new double[][] { coords.get("CA"), CBCoord, new double[3] } );
                //last one is arbitrary (6th DOF will be handled by gen chi1 adjustment)
                
                for(int atomIndex=0; atomIndex<res.atoms.size(); atomIndex++){
                    if( ! BBAtomNames.contains(res.atoms.get(atomIndex).name) ){
                        //sidechain atom (or HA)
                        sidechainMotion.transform(res.coords, atomIndex);
                    }
                }
            }
        }
        
        for(String atomName : BBAtomNames){
            int atomIndex = res.getAtomIndexByName(atomName);
            double atomCoords[] = coords.get(atomName);
            
            if(atomIndex==-1){
                if(res.fullName.startsWith("PRO") && atomName.equalsIgnoreCase("H")){//mutating to PRO...use H for CD
                    atomIndex = res.getAtomIndexByName("CD");
                    atomCoords = atomCoords.clone();//will modify to make CD coords
                    rescaleBondLen(atomCoords, coords.get("N"), 1./HNProRatio);
                }
                else
                    throw new RuntimeException("ERROR: Didn't find atom "+atomName+" in residue "+res.fullName);
            }
            
            System.arraycopy(atomCoords, 0, res.coords, 3*atomIndex, 3);
        }
        
    }
    
    
    
    void rescaleBondLen(double coord[], double refCoord[], double distRatio){
        //Rescale coord to or away from refCoord so the coord-refCoord distance
        //is multiplied by distRatio
        for(int dim=0; dim<3; dim++)
            coord[dim] = refCoord[dim] + distRatio*(coord[dim]-refCoord[dim]);
    }

}
