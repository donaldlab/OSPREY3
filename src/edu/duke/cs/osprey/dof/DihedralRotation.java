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
