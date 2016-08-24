/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
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
