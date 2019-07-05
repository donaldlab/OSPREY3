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

package edu.duke.cs.osprey.dof.deeper.perts;

import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.Protractor;
import edu.duke.cs.osprey.tools.RigidBodyMotion;
import edu.duke.cs.osprey.tools.RotationMatrix;
import edu.duke.cs.osprey.tools.VectorAlgebra;
import java.util.ArrayList;

/**
 *
 * A backrub perturbation
 * 
 * @author mhall44
 */
public class Backrub extends Perturbation {
    
    final static double thetaSmallScale = 0.7f; //the scaling factor for the rotation angles for the small rotations

    
    public Backrub(ArrayList<Residue> resDirectlyAffected) {
        super(resDirectlyAffected);
    }
    
    
    @Override
    public boolean doPerturbationMotion(double paramVal) {
        //Apply the perturbation
        //Use an arbitrary param (primary backrub angle in degrees)
        //Don't store rotation matrices or translations

        RigidBodyMotion[] rotations = calcTransRot(paramVal);
        applyBackrubLikeMotion(rotations);
        return true;//we can always do a backrub
    }
    
    
    
    public RigidBodyMotion[] calcTransRot(double param){
        //Calculate the two rigid-body motions that make up
        //the action of the backrub on the backbone
        
        //get coordinates of atoms we need
        //Calpha coordinates first
        double[][] x=new double[3][];
        for(int a=0;a<3;a++){
            x[a] = resDirectlyAffected.get(a).getCoordsByAtomName("CA");
        }

        double[] oldO1Coord = resDirectlyAffected.get(0).getCoordsByAtomName("O");
        double[] oldO2Coord = resDirectlyAffected.get(1).getCoordsByAtomName("O");

        double rotax[]=VectorAlgebra.subtract(x[2],x[0]);//vector from 1st to second Calpha: primary rotation axis

        //Create the primary rotation matrix.  This can be used to rotate about either anchor CA.  
        RotationMatrix rmPrimary = new RotationMatrix(rotax[0],rotax[1],rotax[2],param,false);
        RigidBodyMotion primaryRot = new RigidBodyMotion(x[0], rmPrimary, x[0]);
        
        double newMidCACoord[] = x[1].clone();
        primaryRot.transform(newMidCACoord);
        //Rotated middle CA coordinates
        
        double newO1Coord[] = oldO1Coord.clone();
        primaryRot.transform(newO1Coord);
        double newO2Coord[] = oldO2Coord.clone();
        primaryRot.transform(newO2Coord);
        //Rotated oxygen coordinates

       
        //Compute the small rotation for the first peptide plane, and compute
        //the rotation matrix for it
        double theta = getSmallRotAngle(newO1Coord,newMidCACoord,x[0],oldO1Coord);
        theta *= thetaSmallScale;
        if (Math.signum(theta)==Math.signum(param))
            theta = -theta;
        rotax = VectorAlgebra.subtract( x[1], x[0] );
        RotationMatrix firstSmallRot = new RotationMatrix( rotax[0], rotax[1], rotax[2], theta, false );
        RotationMatrix firstPepRM = firstSmallRot.multiply(rmPrimary);
        RigidBodyMotion firstPepRot = new RigidBodyMotion(x[0], firstPepRM, x[0]);
        //first peptide rotation can be performed about first CA (x[0])
        

        theta = getSmallRotAngle(newO2Coord,newMidCACoord,x[2],oldO2Coord);
        theta *= thetaSmallScale;
        if (Math.signum(theta)==Math.signum(param))
            theta = -theta;
        rotax = VectorAlgebra.subtract( x[2], x[1] );
        
        RotationMatrix secondSmallRot = new RotationMatrix( rotax[0], rotax[1], rotax[2], theta, false );
        RotationMatrix secondPepRM = secondSmallRot.multiply(rmPrimary);
        RigidBodyMotion secondPepRot = new RigidBodyMotion(x[2], secondPepRM, x[2]);
        //second peptide rotation can be performed about last CA (x[2])
        
        return new RigidBodyMotion[] {firstPepRot,secondPepRot};
    }
 
    
    //Get the small rotation angle that will rotate atom pp1 around the axis defined by atoms (pp2,a3), so that pp1 will be as close as possible to atom a4
    //Adapted from Backrubs to use only Atom.coord arrays (needed for the setup of the calculations here)
    private double getSmallRotAngle(double[] pp1, double[] pp2, double[] a3, double[] a4){

        double[] pp3 = projectPointLine(a3, pp2, pp1);
        double[] pp4 = projectPointPlane(a3, pp2, pp3, a4);
        double[] closestPoint = getClosestPoint(pp3,pp1,pp4);

        return Protractor.getAngleDegrees(pp1, pp3, closestPoint);
    }
    
    
    
    
    
    //The next three function are from Backrubs; they now operate on coordinates (double[])
	
    //Returns the projection of atom p3 onto the line between atoms l1 and l2
    public double[] projectPointLine(double[] l1, double[] l2, double[] p3){

            double c[] = new double[3];
            double d12sq = VectorAlgebra.normsq(VectorAlgebra.subtract(l1,l2));
            
            double u = ( (p3[0]-l1[0])*(l2[0]-l1[0]) + (p3[1]-l1[1])*(l2[1]-l1[1]) + (p3[2]-l1[2])*(l2[2]-l1[2]) );
            u = u / d12sq;

            c[0] = (double)(l1[0] + u * (l2[0]-l1[0]));
            c[1] = (double)(l1[1] + u * (l2[1]-l1[1]));
            c[2] = (double)(l1[2] + u * (l2[2]-l1[2]));

            return c;
    }

    //Returns the projection (a pseudo-atom) of atom p4 onto the plane with normal defined by atoms (l1,l2) and a point on that plane c3
    public double[] projectPointPlane(double[] l1, double[] l2, double[] c3, double[] p4){

            double l[] = new double[3];
            for (int i=0; i<3; i++)
                    l[i] = l2[i] - l1[i];

            double d = 0.0f;
            for (int i=0; i<3; i++)
                    d += l[i]*c3[i];

            double t1 = 0.0f;
            for (int i=0; i<3; i++)
                    t1 += l[i]*p4[i];
            t1 -= d;

            double t2 = 0.0f;
            for (int i=0; i<3; i++)
                    t2 += l[i]*l[i];

            double t = t1/t2;

            double r[] = new double[3];
            for (int i=0; i<3; i++)
                    r[i] = p4[i] - t*l[i];

            return r;
    }

    //Compute the closest point on the circle defined by the Atom c (center) and radius (c,p1) to Atom q1 (this assumes c, p1, and q1 are coplanar)
    public double[] getClosestPoint (double[] c, double[] p1, double[] q1){
            
            double r = VectorAlgebra.distance(c, p1);
        
            double t[] = new double[3];
            for (int i=0; i<3; i++)
                    t[i] = q1[i] - c[i];

            double d = 0.0f;
            for (int i=0; i<3; i++)
                    d += t[i]*t[i];
            d = (double)Math.sqrt(d);

            double a[] = new double[3];
            for (int i=0; i<3; i++)
                    a[i] = c[i] + r*(t[i]/d);

            return a;
    }
    
    @Override
    public Perturbation copyForNewMolecule(Molecule mol, PerturbationBlock block){
        Backrub br = new Backrub(Residue.equivalentInMolec(resDirectlyAffected, mol));
        br.curParamVal = curParamVal;
        br.indexInBlock = indexInBlock;
        br.block = block;
        return br;
    }
    
}
