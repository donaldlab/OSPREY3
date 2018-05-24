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

package edu.duke.cs.osprey.tools;

/**
 *
 * @author mhall44
 */
public class RigidBodyMotion {
    //This class can represent any rigid-body motion
    //subtract center1, apply rotation about the origin, and then add center2
    
    double center1[];
    RotationMatrix rotation;
    double center2[];

    
    public RigidBodyMotion(double[] center1, RotationMatrix rotation, double[] center2) {
        this.center1 = center1;
        this.rotation = rotation;
        this.center2 = center2;
    }
    
    
    
    public RigidBodyMotion(double[] center, double[] axis, double angle, boolean angleInRadians){
        //rotate about center with the specified axis and angle (in degrees or radians as specified)
        center1 = center;
        center2 = center;
        
        rotation = new RotationMatrix(axis[0],axis[1],axis[2],angle,angleInRadians);
    }
    

    public RigidBodyMotion(double[] center, double[] axis, double sinAngle, double cosAngle){
        //rotate about center with the specified axis; sine and cosine of angle given separately
        center1 = center;
        center2 = center;
        
        rotation = new RotationMatrix(axis[0],axis[1],axis[2],sinAngle,cosAngle);
    }
    
    public RigidBodyMotion(double[][] initCoords, double[][] finalCoords){
        //we're given three sets of 3-D coordinates, in an initial and a final state
        //we superimpose them, matching the first pair of coordinates exactly, then the direction for the
        //difference between the second pair exactly
        center1 = initCoords[0].clone();
        center2 = finalCoords[0].clone();
        
        //vectors to superimpose by a rotation
        double uold[] = VectorAlgebra.subtract(initCoords[1], initCoords[0]);
        double unew[] = VectorAlgebra.subtract(finalCoords[1], finalCoords[0]);
        double vold[] = VectorAlgebra.subtract(initCoords[2], initCoords[0]);
        double vnew[] = VectorAlgebra.subtract(finalCoords[2], finalCoords[0]);
        
        rotation = RotationMatrix.getSuperposingRotMatrix(uold, unew, vold, vnew);
    }
    
    
    public void transform(double[] concatCoords){
        //given a bunch of concatenated 3D vectors, apply the transformation to each of them, in place
        int numVectors = concatCoords.length/3;
        
        for(int v=0; v<numVectors; v++)
            transform(concatCoords,v);
    }
    
    public void transform(double[] concatCoords, int index){
        //transform only the vector in concatCoords with the specified index

        for(int a=0; a<3; a++)
            concatCoords[3*index+a] -= center1[a];

        rotation.applyRotation(concatCoords, index);

        for(int a=0; a<3; a++)
            concatCoords[3*index+a] += center2[a];
    }
    
}
