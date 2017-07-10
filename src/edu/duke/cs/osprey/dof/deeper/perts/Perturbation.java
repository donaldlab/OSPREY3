/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.dof.deeper.perts;

import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.dof.DOFBlock;
import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.dof.deeper.GenChi1Calc;
import edu.duke.cs.osprey.dof.deeper.ResBBState;
import edu.duke.cs.osprey.dof.deeper.SidechainIdealizer;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.RigidBodyMotion;
import java.util.ArrayList;

/**
 *
 * A perturbation, as in DEEPer
 * 
 * Note: This DEEPer package, unlike the rest of the software,
 * consists of perturbations specific to proteins
 * and so it will make explicit reference to protein atoms outside HardCodedResidueInfo
 * DEEPer perturbations can be used in the protein part of a mixed protein-nonprotein system
 * 
 * @author mhall44
 */
public abstract class Perturbation extends DegreeOfFreedom {
    
    PerturbationBlock block = null;//block to which this perturbation belongs.  
    //Keeps track of non-commutativity issues & sidechain idealization
    //Perturbation needs to be assigned to a block before being used.
    
    int indexInBlock = -1;//index of this perturbation in block
    
    
    double curParamVal = 0;
    
    
    //this field should be set by the constructor (the three above need not be
    //as long as a PerturbationBlock is being created)
    ArrayList<Residue> resDirectlyAffected;

    
    
    public Perturbation(ArrayList<Residue> resDirectlyAffected) {
        this.resDirectlyAffected = resDirectlyAffected;
    }
    
    
    @Override
    public void apply(double paramVal) {
        
        if(paramVal==curParamVal)
            return;//no change needed
        
        curParamVal = paramVal;
        
        ArrayList<Residue> dependentResidues = block.dependentResidues.get(indexInBlock);
        //residues that are moved by this perturbation or its successors
        
        //OK now we will revert all successors to the state they were at before this perturbation
        //and any successors were applied
        //The backbone is reverted exactly, and the sidechain follows as a rigid body
        ArrayList<Double> dependentGenChi1 = new ArrayList<>();
        for(Residue res : dependentResidues){
            dependentGenChi1.add( GenChi1Calc.getGenChi1(res) );//record gen chi1 so we can restore it later
            ResBBState prePertState = block.prePertBBStates.get(indexInBlock).get(res);
            prePertState.putInState(res);//revert BB atoms
        }
        
        //OK now we can actually apply the perturbation motion, and update pre-pert states for any successors
        doPerturbationMotion(paramVal);
        block.updateSuccessorPrePertStates(indexInBlock);
        
        //OK we now have to restore the correct parameter values for each of the successors,
        //and if they have successors the pre-pert states need to be updated too
        
        for(Perturbation successor : block.successors.get(indexInBlock)){
            successor.doPerturbationMotion(successor.curParamVal);
            block.updateSuccessorPrePertStates(successor.indexInBlock);
        }
        
        //OK now the backbone atoms are all correct!
        //We now place the sidechains (by idealization, and appropriate setting of gen chi1)
        //other aspects of sidechain geometry (chi2, etc.) will be correct because the sidechain
        //is treated as a rigid body (except for Pro, which has no other sidechain DOFs)
        
        for(int resNum=0; resNum<dependentResidues.size(); resNum++){
            Residue res = dependentResidues.get(resNum);
            SidechainIdealizer.idealizeSidechain(EnvironmentVars.resTemplates, res);
            GenChi1Calc.setGenChi1(res, dependentGenChi1.get(resNum));
        }
    }
    
    
    
    public abstract boolean doPerturbationMotion(double paramVal);
    //actually handle the perturbation motion, moving the backbone atoms to the desired position
    //If the perturbation was geometrically impossible, return false
    //and mark any residues with messed-up conformations as invalidBB
    //This motion should move the sidechain + CA as a rigid body
    //(this way, following up with a sidechain idealization rigid-body motion
    //defined by the new CA position, idealized CB position, and gen chi1
    //will be sure to get the sidechain in the right pose)
    
    
    
    void movePeptidePlane(RigidBodyMotion motion, int startingRes, boolean includeFinalSCH){
        //Apply motion to the peptide plane between 
        //resDirectlyAffected[startingRes] and resDirectlyAffected[startingRes+1]
        //Transform the sidechain, CA, and HA of the latter if indicated
        //Used in several perturbations
        

        //handle carbonyl of first residue
        Residue firstRes = resDirectlyAffected.get(startingRes);
        for(String atomName : new String[] {"C","O"})
            motion.transform( firstRes.coords, firstRes.getAtomIndexByName(atomName) );
        
        //handle amide and possibly CA of second residue
        Residue secondRes = resDirectlyAffected.get(startingRes+1);
        
        
        if(includeFinalSCH){
            //move everything except the carbonyl
            for(int atomIndex=0; atomIndex<secondRes.atoms.size(); atomIndex++){
                
                String atomName = secondRes.atoms.get(atomIndex).name;
                
                if( ! ( atomName.equalsIgnoreCase("C")
                        || atomName.equalsIgnoreCase("O")
                        || atomName.equalsIgnoreCase("OXT") ) ) {
                
                    motion.transform( secondRes.coords, atomIndex );
                }
            }
        }
        else {
            ArrayList<String> res2PepAtoms = new ArrayList<>();//list of atoms in 2nd res to move
            res2PepAtoms.add("N");

            //secondRes will have either an amide H, or be Pro (CD takes the place of H)
            if(secondRes.template.name.equalsIgnoreCase("PRO"))
                res2PepAtoms.add("CD");
            else
                res2PepAtoms.add("H");


            //now move second-res atoms
            for(String atomName : res2PepAtoms)
                motion.transform( secondRes.coords, secondRes.getAtomIndexByName(atomName) );
        }
    }
    
    
    void applyBackrubLikeMotion(RigidBodyMotion[] pepRots){
        //The backrub motion acts on the backbone with two peptide-plane rotations
        //(both already composed with a primary rotation about the CA0-CA2 axis)
        //This motion can be used to apply both the Backrub and LoopClosureAdjustment
        movePeptidePlane(pepRots[0], 0, true);
        movePeptidePlane(pepRots[1], 1, false);
    }
    
    
    
    public abstract Perturbation copyForNewMolecule(Molecule mol, PerturbationBlock block);
    //copy the perturbation to apply to mol; block is the perturbation block for mol
 
    @Override
    public DOFBlock getBlock(){
        return block;
    }
    
    @Override
    public String getName() {
        return "PERT"+block.allResidues.get(0).getPDBResNumber()+"."+indexInBlock;
    }
}
