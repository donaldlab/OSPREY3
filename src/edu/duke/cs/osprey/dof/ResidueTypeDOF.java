/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.dof;

import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.dof.deeper.SidechainIdealizer;
import edu.duke.cs.osprey.restypes.HardCodedResidueInfo;
import edu.duke.cs.osprey.restypes.ResidueTemplate;
import edu.duke.cs.osprey.restypes.GenericResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.Atom;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.RigidBodyMotion;
import java.util.ArrayList;

/**
 *
 * @author mhall44
 */
public class ResidueTypeDOF extends DegreeOfFreedom {
    //This degree of freedom is the residue type (e.g., AA type) at a particular position
    //So applying values of it means mutating the residue
    
    Residue res;//what residue in the molecule we are talking about

    
    
    public ResidueTypeDOF(Residue res) {
        this.res = res;
    }
    
    
    public void mutateTo(String resType) {
        //paramVal is the index in the ResidueTemplateLibrary of the new parameter type
        //so it must be an integer...
        GenericResidueTemplateLibrary templateLib = EnvironmentVars.resTemplates;
        
        ResidueTemplate oldTemplate = res.template;
        ResidueTemplate newTemplate = templateLib.getTemplateForMutation(resType,res,true);
        
        //the residue's going to change some, so break its inter-residue bonds
        res.removeInterResBonds();
        res.intraResBondsMarked = false;//we'll need to redo these too
        
        res.template = newTemplate;
        
        res.fullName = newTemplate.name + res.fullName.substring(3);
        //res type name is first three characters of full name
        
        
        //coordinates will come from the template,
        //but we'll move them as a rigid body to match the backbone atoms
        int[][] mutAlignAtoms = HardCodedResidueInfo.findMutAlignmentAtoms(oldTemplate,newTemplate);
        //2x4 array of which atoms are used in the old and new residues to align
        //to do the mutation
        
        double oldMAACoords[][] = extractCoords(mutAlignAtoms[0],res.coords);
        //coordinates of the atoms in the old residue to align
        
        double newCoords[] = newTemplate.templateRes.coords.clone();
        double templateMAACords[][] = extractCoords(mutAlignAtoms[1],newCoords);
        
        //we now construct a rigid-body motion that will map the sidechain (or generally,
        //the non-"backbone" part, if not a standard amino acid) from the template to the residue's
        //curent frame of reference
        RigidBodyMotion templMotion = new RigidBodyMotion(templateMAACords,oldMAACoords);
        
        templMotion.transform(newCoords);
        
        //the backbone atoms will be kept exactly as before the mutation
        //if the sidechain attaches only to the first mutAlignAtom, this method keeps bond lengths
        //exactly as in the template for sidechain, and as in the old backbone otherwise
        ArrayList<String> BBAtomNames =  HardCodedResidueInfo.listBBAtomsForMut(newTemplate,oldTemplate);
        for(String BBAtomName : BBAtomNames){
            int BBAtomIndexOld = oldTemplate.templateRes.getAtomIndexByName(BBAtomName);
            int BBAtomIndexNew = newTemplate.templateRes.getAtomIndexByName(BBAtomName);
            
            //copy coordinates of the BB atom from old to new coordinates
            System.arraycopy(res.coords, 3*BBAtomIndexOld, newCoords, 3*BBAtomIndexNew, 3);
        }
        
        res.coords = newCoords;
        
        //finally, update atoms in res to match new template
        ArrayList<Atom> newAtoms = new ArrayList<>();
        for(Atom at : newTemplate.templateRes.atoms){
            Atom newAtom = at.copy();
            newAtom.res = res;
            newAtoms.add(newAtom);
        }
        res.atoms = newAtoms;    
        
        //reconnect all bonds
        res.markIntraResBondsByTemplate();
        HardCodedResidueInfo.reconnectInterResBonds(res);
        
        //special case if sidechain loops back in additional place to backbone...
        if(oldTemplate.name.equalsIgnoreCase("PRO") || newTemplate.name.equalsIgnoreCase("PRO")){
            SidechainIdealizer.idealizeSidechain(res);
            if(!newTemplate.name.equalsIgnoreCase("PRO")){//if mutating from Pro, no ring closure issues possible anymore
                if(res.pucker!=null){
                    if(res.pucker.puckerProblem != null){
                        res.pucker.puckerProblem.removeFromRes();
                        res.pucker.puckerProblem = null;
                    }
                }
            }
        }
    }
    
    
    static double[][] extractCoords(int[] index, double[] allCoords){
        //allCoords is the concatenated coords of a bunch of atoms
        //get the coords (3-D) of the atoms with the specified indices
        double[][] ans = new double[index.length][3];
        
        for(int i=0; i<index.length; i++){
            System.arraycopy(allCoords, 3*index[i], ans[i], 0, 3);
        }
        
        return ans;
    }
    
    
    public String getCurResType(){//current residue type
        return res.fullName.substring(0,3);//this is always the first three letters of the full name
    }
    

    @Override
    public void apply(double paramVal) {
        throw new IllegalArgumentException("ERROR: ResidueTypeDOF takes an residue type name"
                + " as argument; can't take "+paramVal);
    }
    
    
    
    
    
    
}
