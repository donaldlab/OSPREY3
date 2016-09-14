/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.dof;

import edu.duke.cs.osprey.structure.Molecule;
import java.util.LinkedHashMap;

/**
 * 
 * A block of degrees of freedom that are applied together
 * 
 * @author mhall44
 */
public interface DOFBlock {

    //Make a copy of the block that operates in a new molecule, mol.  
    //For each DOF in the block, add an entry (DOF -> copy of DOF) to copiedDOFMap
    public DOFBlock copyForNewMolecule(Molecule mol, LinkedHashMap<DegreeOfFreedom, DegreeOfFreedom> copiedDOFMap);
    
}
