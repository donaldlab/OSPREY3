/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.confspace;

import edu.duke.cs.osprey.bbfree.BBFreeBlock;
import edu.duke.cs.osprey.bbfree.BBFreeDOF;
import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.dof.MoveableStrand;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.dof.deeper.perts.Perturbation;
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
    
    public PositionConfSpace(Residue res, ArrayList<DegreeOfFreedom> resDOFs, ArrayList<String> allowedAAs, 
            boolean contSCFlex, ArrayList<DegreeOfFreedom> strandDOFs, 
            ArrayList<Perturbation> perts, ArrayList<ArrayList<double[]>> pertIntervals, 
            ArrayList<ArrayList<int[]>> pertStates, BBFreeBlock bfb ) {
        
        //We'll start with just one RC for each rotamer
        //But in general there are a lot of options for RCs...
        
        ResidueTemplateLibrary templateLib = EnvironmentVars.resTemplates;
        this.res = res;
        
        
        if(pertStates==null){//no DEEPer flexibility...
            pertStates = new ArrayList<>();
            pertStates.add(null);
        }
        
        
        for(String AAType : allowedAAs){
            int numDihedrals = templateLib.numDihedralsForResType(AAType);
            int numRot = templateLib.numRotForResType(AAType);
            
            //resDOFs is all sidechain DOFs, for now
            ArrayList<DegreeOfFreedom> dofListForRot = new ArrayList<>();
            for(int dih=0; dih<numDihedrals; dih++)//get the first numDihedrals dihedrals
                dofListForRot.add(resDOFs.get(dih));
                       
            
            for(ArrayList<int[]> pertState : pertStates){
                
                if(AAType.equalsIgnoreCase("PRO")){//special case: has ring pucker
                    //If PRO is set to be flexible we'll assume this includes pucker flexibility
                    //(the main flexibility of proline)
                    for( int puckerVal : new int[]{0,1} ){//create both puckers
                        createRC(0, AAType, -1, contSCFlex, dofListForRot, puckerVal,
                                strandDOFs, bfb, pertState, perts, pertIntervals);
                    }
                }
                
                else if(numRot==0){//ALA or GLY: no rotamers or dihedrals, so create a single rigid RC (or one for each pert state)
                    createRC(0, AAType, -1, contSCFlex, dofListForRot, -1, strandDOFs, bfb, 
                            pertState, perts, pertIntervals);
                }
                
                else {
                    for(int rot=0; rot<numRot; rot++){
                        createRC(numDihedrals, AAType, rot, contSCFlex, dofListForRot, -1, strandDOFs, 
                                bfb, pertState, perts, pertIntervals);
                    }
                }
            }
        }
        
    }
    
    
    private void createRC(int numDihedrals, String AAType, int rot, boolean contSCFlex, ArrayList<DegreeOfFreedom> dofListForRot,
            int proPucker, ArrayList<DegreeOfFreedom> strandDOFs, BBFreeBlock bfb,
            ArrayList<int[]> pertState, ArrayList<Perturbation> perts, ArrayList<ArrayList<double[]>> pertIntervals){
        
        //Create an RC with the specified rotamer and (if set) perturbation state
        
        ResidueTemplateLibrary templateLib = EnvironmentVars.resTemplates;
        //create RC
        ArrayList<Double> dofLB = new ArrayList<>();//lower bounds on each DOF for this RC
        ArrayList<Double> dofUB = new ArrayList<>();//upper bounds
        
        //we'll start with the sidechain dihedral DOFs
        
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
        
        ArrayList<DegreeOfFreedom> dofListForRC = new ArrayList<>();
        dofListForRC.addAll(dofListForRot);
        
        //Put in proline pucker if this is a proline
        if(AAType.equalsIgnoreCase("PRO")){
            dofListForRC.add(res.pucker);
            dofLB.add( (double) proPucker );
            dofUB.add( (double) proPucker );
        }
        
        //add strand degrees of freedom if any
        //these are DOFs whose strand includes this residue!
        for(DegreeOfFreedom strandDOF : strandDOFs){
            dofListForRC.add(strandDOF);
            double[] bounds = MoveableStrand.getStrandDOFBounds(strandDOF);
            dofLB.add(bounds[0]);
            dofUB.add(bounds[1]);
        }
        
        if(bfb!=null){//This res is in a free-backbone block
            ArrayList<BBFreeDOF> bbFreeDOFs = bfb.getDOFs();
            double freeDOFVoxel[][] = bfb.getFreeDOFVoxel();
            
            for(int dnum=0; dnum<bbFreeDOFs.size(); dnum++){
                dofListForRC.add(bbFreeDOFs.get(dnum));
                dofLB.add(freeDOFVoxel[0][dnum]);
                dofUB.add(freeDOFVoxel[1][dnum]);
            }
        }
            
        
        if(pertState != null) {
            //need to add DEEPer DOFs
            
            for(int[] singlePertState : pertState){
                int pertNum = singlePertState[0];//number (in perts) of the perturbation we're adding
                int pertStateNum = singlePertState[1];//state of this perturbation
                dofListForRC.add(perts.get(pertNum));
                
                double[] pertInterval = pertIntervals.get(pertNum).get(pertStateNum);//interval for this perturbation in this state
                dofLB.add(pertInterval[0]);
                dofUB.add(pertInterval[1]);
            }
        }
        
        RC newRC = new RC(AAType, rot, dofListForRC, dofLB, dofUB, RCs.size());
        
        RCs.add(newRC);
    }
    
    
    
}
