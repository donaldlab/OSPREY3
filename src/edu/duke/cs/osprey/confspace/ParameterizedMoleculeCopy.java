/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.confspace;

import edu.duke.cs.osprey.dof.DOFBlock;
import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.structure.Molecule;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;

/**
 *
 * Sometimes we will want to minimize conformations in parallel,
 * so each thread will need a copy of the molecule and its DOFs.  
 * Everything that might change when we change the DOF values is copied deeply.  
 * But the ConfSpace will still have the original DOFs,
 * so we also need a map from original to copied DOFs.  
 * 
 * @author mhall44
 */
public class ParameterizedMoleculeCopy {
    
    private Molecule copiedMol = null;
    private final LinkedHashMap<DegreeOfFreedom,DegreeOfFreedom> original2CopiedDOF = new LinkedHashMap<>();
    
    public ParameterizedMoleculeCopy(ConfSpace confSpace){
        //Make the copy
        //Molecule and DOFs must be copied together because for some DOFs (e.g., strand rotation/translation),
        //the DOF has a bunch of internal state that will be messed up if it gets out of sync with the molecule's movement
        copiedMol = new Molecule(confSpace.m);
        
        ArrayList<DegreeOfFreedom> origDOFs = confSpace.listAllDOFs();

        //map of DOF blocks in confSpace to their copies
        HashMap<DOFBlock,DOFBlock> newBlocks = new HashMap<>();

        
        for(DegreeOfFreedom dof : origDOFs){
            DOFBlock block = dof.getBlock();
            if(block==null) {//standalone DOF
                DegreeOfFreedom copiedDOF = dof.copy();
                copiedDOF.setMolecule(copiedMol);
                original2CopiedDOF.put(dof, copiedDOF);
            }
            else {//DOF is part of a block
                if(!newBlocks.containsKey(block)){//first DOF in its block
                    DOFBlock copiedBlock = block.copyForNewMolecule(copiedMol, original2CopiedDOF);
                    //this will copy all the DOFs in the block at once
                    newBlocks.put(block,copiedBlock);
                }
            }
        }
    }
    
    
    //making without copying (used for test purposes)
    private ParameterizedMoleculeCopy(){}
    
    public static ParameterizedMoleculeCopy makeNoCopy(ConfSpace confSpace){
        ParameterizedMoleculeCopy pmc = new ParameterizedMoleculeCopy();
        pmc.copiedMol = confSpace.m;
        for(DegreeOfFreedom dof : confSpace.listAllDOFs()){
            pmc.original2CopiedDOF.put(dof, dof);
        }
        return pmc;
    }
    
    
    
    public Molecule getCopiedMolecule(){
        return copiedMol;
    }
    
    public DegreeOfFreedom getCopiedDOF(DegreeOfFreedom origDOF){
        return original2CopiedDOF.get(origDOF);
    }
    
}
