/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.confspace;

import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.Residue;
import java.io.Serializable;
import java.util.ArrayList;

/**
 *
 * @author mhall44
 */
public class PositionConfSpace implements Serializable {
    //This class defines the conformational space of a flexible residue
    //including allowed amino-acid types, and rotamers/RCs for each type 
    //subclass PositionConfSpace to make super-residues with super-RCs...
    
    
    public ArrayList<RC> RCs = new ArrayList<>();
    
    public Residue res;//The residue involved
    
    static double dihedFlexInterval = 9;// +/- 9 degree sidechain dihedral continuous flexibility...
    //later can allow this to vary across different dihedrals
    
    public PositionConfSpace(Residue res, ArrayList<DegreeOfFreedom> resDOFs, ArrayList<String> allowedAAs, boolean contSCFlex) {
        
        //We'll start with just one RC for each rotamer
        //But in general there are a lot of options for RCs...
        
        ResidueTemplateLibrary templateLib = EnvironmentVars.resTemplates;
        this.res = res;
        
        for(String AAType : allowedAAs){
            int numDihedrals = templateLib.numDihedralsForResType(AAType);
            int numRot = templateLib.numRotForResType(AAType);
            
            //resDOFs is all sidechain DOFs, for now
            ArrayList<DegreeOfFreedom> dofListForRC = new ArrayList<>();
            for(int dih=0; dih<numDihedrals; dih++)//get the first numDihedrals dihedrals
                dofListForRC.add(resDOFs.get(dih));
            
            if(numRot==0){//ALA or GLY: no rotamers or dihedrals, so create a single rigid RC
                RC newRC = new RC(AAType, -1, dofListForRC, new ArrayList<Double>(), new ArrayList<Double>(), RCs.size());
                RCs.add(newRC);
            }
                
            for(int rot=0; rot<numRot; rot++){
                //create RC
                ArrayList<Double> dofLB = new ArrayList<>();//lower bounds on each DOF for this RC
                ArrayList<Double> dofUB = new ArrayList<>();//upper bounds

                for(int dih=0; dih<numDihedrals; dih++){
                    double dihedralValForRot = templateLib.getDihedralForRotamer(AAType,rot,dih);
                    
                    if(contSCFlex){//allow continuous flexibility up to dihedFlexInterval in each direction
                        dofLB.add(dihedralValForRot-dihedFlexInterval);
                        dofUB.add(dihedralValForRot+dihedFlexInterval);
                    }
                    else {
                        dofLB.add(dihedralValForRot);
                        dofUB.add(dihedralValForRot);
                    }
                }
                
                RC newRC = new RC(AAType, rot, dofListForRC, dofLB, dofUB, RCs.size());
                RCs.add(newRC);
            }
        }
        
    }
    
    
    
}
