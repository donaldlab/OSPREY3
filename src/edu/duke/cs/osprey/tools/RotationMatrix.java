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

import java.io.Serializable;

/**x
 *
 * @author mhall44
 */
public class RotationMatrix implements Serializable {
    
	private static final long serialVersionUID = 2915420374293280379L;
	double[][] matrix;
    
    
    public RotationMatrix(double[][] mtx){
        matrix = mtx;
    }
    
    // This function constructs a rotation matrix from a rotation in
	//  axis-angle notation.  (fx,fy,fz) is the axis (not necessarily normalized)
    public RotationMatrix (double fx, double fy, double fz, double angle, boolean angleInRadians) {//if not radians, then degrees
        this(
            fx, fy, fz,
            angleInRadians ? Math.sin(angle) : Math.sin(Math.toRadians(angle)),
            angleInRadians ? Math.cos(angle) : Math.cos(Math.toRadians(angle))
        );
    }
        
        
    // rotation matrix from axis (fx,fy,fz) angle (sin,cos) representation
    public RotationMatrix (double fx, double fy, double fz, double sin, double cos) {

        /* this implementation is numerically unstable
            let's try something simpler
        
        // First convert the axisangle to a quaternion
        //use half-angle formulae
        double[] halfAngSC = Protractor.getHalfAngleSinCos(sinTheta,cosTheta);
        double sin_a = halfAngSC[0];
        double cos_a = halfAngSC[1];
        double tmp = (double) Math.sqrt(fx*fx + fy*fy + fz*fz);
        double qx = fx / tmp * sin_a;
        double qy = fy / tmp * sin_a;
        double qz = fz / tmp * sin_a;
        double qw = cos_a;
        tmp = (double) Math.sqrt(qx*qx + qy*qy + qz*qz + qw*qw);
        qx /= tmp;
        qy /= tmp;
        qz /= tmp;
        qw /= tmp;
        double xx = qx * qx;
        double xy = qx * qy;
        double xz = qx * qz;
        double xw = qx * qw;

        double yy = qy * qy;
        double yz = qy * qz;
        double yw = qy * qw;

        double zz = qz * qz;
        double zw = qz * qw;

        matrix = new double[3][3];
        
        matrix[0][0] = 1 - 2 * (yy + zz);
        matrix[0][1] = 2 * (xy - zw);
        matrix[0][2] = 2 * (xz + yw);

        matrix[1][0] = 2 * (xy + zw);
        matrix[1][1] = 1 - 2 * (xx + zz);
        matrix[1][2] = 2 * (yz - xw);

        matrix[2][0] = 2 * (xz - yw);
        matrix[2][1] = 2 * (yz + xw);
        matrix[2][2] = 1 - 2 * (xx + yy);
        */
        
        // normalize f
        double len = Math.sqrt(fx*fx + fy*fy + fz*fz);
        double ux = fx/len;
        double uy = fy/len;
        double uz = fz/len;
        
        double omcos = 1 - cos;
        double uxyomcos = ux*uy*omcos;
        double uxzomcos = ux*uz*omcos;
        double uyzomcos = uy*uz*omcos;
        double uxsin = ux*sin;
        double uysin = uy*sin;
        double uzsin = uz*sin;
        
        matrix = new double[3][3];
        
        matrix[0][0] = cos + ux*ux*omcos;
        matrix[0][1] = uxyomcos - uzsin;
        matrix[0][2] = uxzomcos + uysin;
        
        matrix[1][0] = uxyomcos + uzsin;
        matrix[1][1] = cos + uy*uy*omcos;
        matrix[1][2] = uyzomcos - uxsin;
        
        matrix[2][0] = uxzomcos - uysin;
        matrix[2][1] = uyzomcos + uxsin;
        matrix[2][2] = cos + uz*uz*omcos;
    }

    public RotationMatrix multiply(RotationMatrix rotation2){
        return new RotationMatrix(multiplyMatrices(matrix,rotation2.matrix));
    }

    private static double[][] multiplyMatrices(double[][] M1, double[][] M2){
        double[][] ans=new double[3][3];
        for(int a=0;a<3;a++){
            for(int b=0;b<3;b++){
                ans[a][b]=0;
                for(int c=0;c<3;c++){
                    ans[a][b]+=M1[a][c]*M2[c][b];
                }
            }
        }
        return ans;
    }


    public RotationMatrix transpose(){
        double[][] ans=new double[3][3];
        for(int a=0;a<3;a++){
            for(int b=0;b<3;b++){
                ans[a][b] = matrix[b][a];
            }
        }
        return new RotationMatrix(ans);
    }
    
    
    
    
    
    public void applyRotation(double[] x, int index){
        //Given the concatenated 3-D coordinates of a bunch of atoms,
        //apply rotation to the atom with the given index
        //(in place)

        double val;
        double ans[] = new double[3];
        for(int a=0;a<3;a++){
            val=0;
            for(int b=0;b<3;b++){
                val+=matrix[a][b]*x[3*index+b];
            }
            ans[a]=val;
        }
        
        System.arraycopy(ans,0,x,3*index,3);
    }
    
    
    public double[] rotateVector(double[] vec){//Apply rotation matrix to vector vec, i.e. compute matrix*vec
        //return ans, keeping vec as before
        double ans[]=new double[3];
        double val;
        for(int a=0;a<3;a++){
            val=0;
            for(int b=0;b<3;b++){
                val+=matrix[a][b]*vec[b];
            }
            ans[a]=val;
        }
        return ans;
    }
    
    public double[] unrotateVector(double[] vec){//Reverse rotation on vector vec,
        //i.e. compute inv(matrix)*vec=transpose(rm)*vec
        //return ans, keeping vec as before
        double ans[]=new double[3];
        double val;
        for(int a=0;a<3;a++){
            val=0;
            for(int b=0;b<3;b++){
                val+=matrix[b][a]*vec[b];
            }
            ans[a]=val;
        }
        return ans;
    }
    
    
    public static RotationMatrix identity(){
        //identity matrix
        double M[][] = new double[3][3];

        for(int a=0;a<3;a++){
            for(int b=0;b<3;b++){
                if(a==b)
                    M[a][b] = 1;
                else
                    M[a][b] = 0;
            }
        }
        
        return new RotationMatrix(M);
    }
    
    
    public static RotationMatrix getSuperposingRotMatrix(double uold[], double unew[], double vold[], double vnew[]){
        //Returns a rotation matrix that rotates vector uold to point in the direction of unew, and vold to point in the direction of vnew
        //Relies on uold-vold and unew-vnew angles being basically the same (no rotation exactly satisfies the requirements if they are different)
        //If the angles are close to equal the rotation matrix will exactly superimpose uold, unew and be close to right for vold, vnew

        double th;

        //First create a matrix rotating uold to unew
        //Axis for this will be uold X unew (sign for this will give an angle < 180 degrees)
        RotationMatrix mtx1;
        double axis1[] = VectorAlgebra.cross(uold, unew);
        if( VectorAlgebra.norm(axis1) == 0 ){//uold and unew are collinear
            if( VectorAlgebra.dot(uold,unew) > 0 )//uold and unew are already superimposed
                mtx1 = identity();
            else{//Need a 180-degree rotation.
                double normal[] = VectorAlgebra.getPerpendicular(uold);
                mtx1 = new RotationMatrix( normal[0], normal[1], normal[2], 180, false );
            }
        }
        else{
            th = Protractor.getAngleRadians(uold,unew);//angle of this first rotation (radians)
            mtx1 = new RotationMatrix( axis1[0], axis1[1], axis1[2], th, true );
        }

        //Now create a matrix to rotate about unew, rotating mtx1*vold to vnew
        //The angle assumption comes in here, because we are using unew as the axis
        //So we need vnew-unew angle = applyMatrix(mtx1,vold)-unew angle (which = vold-uold angle)
        RotationMatrix mtx2;
        double w1[] = VectorAlgebra.perpendicularComponent( mtx1.rotateVector(vold), unew );
        double w2[] = VectorAlgebra.perpendicularComponent( vnew, unew );
        th = Protractor.getAngleRadians(w1,w2);//angle of second rotation (radians)
        if( VectorAlgebra.dot( VectorAlgebra.cross(w1,w2), unew) > 0)
            mtx2 = new RotationMatrix( unew[0], unew[1], unew[2], th, true );
        else
            mtx2 = new RotationMatrix( unew[0], unew[1], unew[2], -th, true );

        return mtx2.multiply(mtx1);//Return the product of these two rotations
    }

    public RotationMatrix copy() {//deep copy
        double mtx[][] = new double[matrix.length][];
        for(int a=0; a<matrix.length; a++)
            mtx[a] = matrix[a].clone();
        return new RotationMatrix(mtx);
    }    
}
