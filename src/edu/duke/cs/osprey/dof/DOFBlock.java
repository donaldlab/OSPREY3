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

