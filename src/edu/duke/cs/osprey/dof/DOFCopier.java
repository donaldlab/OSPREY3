/*
** This file is part of OSPREY.
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation, either version 2 of the License, or
** (at your option) any later version.
** 
** OSPREY is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
*/


package edu.duke.cs.osprey.dof;

import edu.duke.cs.osprey.structure.Molecule;
import java.util.LinkedHashMap;
import java.util.Map;

/**
 *
 * This class copies a set of DOFs to a new molecule, so that different threads can minimize
 * energy in parallel and each can have its own molecule.
 * Some DOF types like dihedrals can be copied one at a time, but others are in blocks
 * that must be handled together.  
 * Some DOFs have internal state, so they must be copied to a molecule
 * with the same sequence and conformation as the one they originally referred to.
 * ParameterizedMoleculeCopy takes care of this (copies molecule and DOFs together).
 * 
 * @author mhall44
 */
/*public class DOFCopier {
    
    
    public DOFCopier(){
        
    }
    
    
    public LinkedHashMap<DegreeOfFreedom,double[]> copyForNewMolecule(
            LinkedHashMap<DegreeOfFreedom,double[]> DOFBounds, Molecule mol){
        //Given a map from DOFs to their bounds (defining a voxel),
        //Create a map from copies of those DOFs in molecule mol to their bounds
        LinkedHashMap<DegreeOfFreedom,double[]> copiedDofs = new LinkedHashMap<>();
        
        //First, we must have copies of all the DOF blocks to their new molecules
        //when we generate these, we will get copies of the perturbations that are in blocks too
        LinkedHashMap<DOFBlock,DOFBlock> newBlocks = new LinkedHashMap<>();
        LinkedHashMap<DegreeOfFreedom,DegreeOfFreedom> copiedDOFMap = new LinkedHashMap<>();
        for(DegreeOfFreedom dof : DOFBounds.keySet()){
            DOFBlock block = dof.getBlock();
            if(block != null){//DOF is part of a block
                if(!newBlocks.containsKey(block)){
                    DOFBlock copiedBlock = block.copyForNewMolecule(mol, copiedDOFMap);
                    newBlocks.put(block,copiedBlock);
                }
            }
        }
        
        for (Map.Entry<DegreeOfFreedom,double[]> entry : DOFBounds.entrySet()) {
            DegreeOfFreedom dof = entry.getKey();
            if(copiedDOFMap.containsKey(dof))//dof has already been copied as part of a block
                dof = copiedDOFMap.get(dof);
            else {//copy it now (this is for a standalone DOF, e.g. sidechain dihedral)
                dof = dof.copy();
                dof.setMolecule(mol);
            }
            copiedDofs.put(dof, entry.getValue());
        }
        return copiedDofs;
    }
    
}*/

