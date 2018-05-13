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

import edu.duke.cs.osprey.tools.RigidBodyMotion;
import edu.duke.cs.osprey.tools.VectorAlgebra;

/**
 *
 * This is a RigidBodyMotion that changes the dihedral defined by some specified atoms
 * 
 * @author mhall44
 */
public class DihedralRotation extends RigidBodyMotion {
    
    
    //the rotation is actually fully defined by the two atoms forming the bond
    //we're rotating around (2nd then 3rd of the 4 atoms defining the dihedral),
    //along with the CHANGE in dihedral angle that we want
    
    public DihedralRotation(double[] atom2Coords, double[] atom3Coords, double dihedralChange ){
        //dihedralChange is in degrees
        super (atom3Coords, VectorAlgebra.subtract(atom3Coords, atom2Coords), dihedralChange, false);
    }
    
    public DihedralRotation(double[] atom2Coords, double[] atom3Coords, double sinDihedralChange,
            double cosDihedralChange){
        super (atom3Coords, VectorAlgebra.subtract(atom3Coords, atom2Coords), sinDihedralChange, cosDihedralChange);
    }
    
}

