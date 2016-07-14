/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey;


import org.junit.BeforeClass;
import org.junit.Test;

import edu.duke.cs.osprey.tools.RigidBodyMotion;
import edu.duke.cs.osprey.tools.RotationMatrix;

/**
 *
 * Testing stuff from the tools package
 * 
 * @author mhall44
 */
public class TestTools extends TestBase {
    
    @Test
    public void RigidMotionTest(){
        //Take 3 random points in 3-D, apply a random rigid motion to them, then try to reconstruct 
        //the rigid motion based on superposition
        double points[][] = new double[3][3];
        
        for(int a=0; a<3; a++){
            for(int b=0; b<3; b++)
                points[a][b] = Math.random()-0.5;
        }
        
        RigidBodyMotion motion = new RigidBodyMotion(new double[3]/*] {Math.random(),Math.random(),Math.random()}*/,
                new double[] {Math.random(),Math.random(),Math.random()}, 360*Math.random(), false);
        
        
        double transformedPoints[][] = new double[3][3];
        for(int a=0; a<3; a++){
            System.arraycopy(points[a],0,transformedPoints[a],0,3);
            motion.transform(transformedPoints[a]);
        }
        
        RigidBodyMotion reconstructedMotion = new RigidBodyMotion(points,transformedPoints);
        
        //now apply reconstructed motion to points and see if they match transformedPoints
        //since reconstructedMotion is supposed to be the same as motion
        double error = 0;
        
        for(int a=0; a<3; a++){
            reconstructedMotion.transform(points[a]);
            for(int b=0; b<3; b++){
                double change = points[a][b] - transformedPoints[a][b];
                error += change*change;
            }
        }
        
        //error will be positive and should be very close to 0
        assert(error<1e-16);
    }
    
    
    
    
    @Test
    public void SuperposingRotMatrixTest(){
        //Take 2 random points in 3-D, apply a random rotation to them about the origin, then try to reconstruct 
        //the rotation based on superposition
        //(i.e., rotation-only version of RigidMotionTest)
        double points[][] = new double[2][3];
        
        for(int a=0; a<2; a++){
            for(int b=0; b<3; b++)
                points[a][b] = Math.random()-0.5;
        }
        
        RotationMatrix rot = new RotationMatrix(Math.random(),Math.random(),Math.random(),360*Math.random(),false);
        
        double transformedPoints[][] = new double[2][3];
        for(int a=0; a<2; a++){
            System.arraycopy(points[a],0,transformedPoints[a],0,3);
            rot.applyRotation(transformedPoints[a],0);
        }
        
        RotationMatrix supRot = RotationMatrix.getSuperposingRotMatrix(points[0],transformedPoints[0],
                points[1],transformedPoints[1]);
        
        
        //now apply reconstructed motion to points and see if they match transformedPoints
        //since reconstructedMotion is supposed to be the same as motion
        double error = 0;
        
        for(int a=0; a<2; a++){
            supRot.applyRotation(points[a],0);
            for(int b=0; b<3; b++){
                double change = points[a][b] - transformedPoints[a][b];
                error += change*change;
            }
        }
        
        assert(error<1e-16);
    }
    
    
}
