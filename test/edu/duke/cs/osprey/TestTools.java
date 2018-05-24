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
