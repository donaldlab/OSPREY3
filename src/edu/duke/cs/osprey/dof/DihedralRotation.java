/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
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
