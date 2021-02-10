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
 * Shear perturbation: small perturbation appropriate for helices
 * 
 * @author mhall44
 */
public class Shear extends Perturbation {

    public Shear(ArrayList<Residue> resDirectlyAffected) {
        super(resDirectlyAffected);
    }
    
    
    @Override
    public boolean doPerturbationMotion(double paramVal) {
        
        RigidBodyMotion[] pepMotions = calcMotions(paramVal);
        
        movePeptidePlane(pepMotions[0], 0, true);
        movePeptidePlane(pepMotions[1], 1, true);
        movePeptidePlane(pepMotions[2], 2, false);
        
        //pepMotions[a] is rm rotation about about CA (orig config), followed by tr
        
        return true;//we can always do a shear
    }
    
    
    
    RigidBodyMotion[] calcMotions(double paramVal){
        
        //the motions will consist of rotations (about CA's) followed by translations
        //(for the latter two)
        //each motion is for a particular peptide plane
        RotationMatrix rm[] = new RotationMatrix[3];
        double tr[][] = new double[2][];
        
        
        double x[][]=new double[4][3];//Calpha coordinates
        for(int a=0;a<4;a++){
            x[a] = resDirectlyAffected.get(a).getCoordsByAtomName("CA");
        }
        
        Residue midRes = resDirectlyAffected.get(1);
        double midOCoord[] = midRes.getCoordsByAtomName("O");
        double midCCoord[] = midRes.getCoordsByAtomName("C");
        double w[] = VectorAlgebra.subtract( midOCoord, midCCoord );
        
        
        double x01[] = VectorAlgebra.subtract(x[1],x[0]);//vector from 1st to second Calpha
        double x12[] = VectorAlgebra.subtract(x[2],x[1]);
        double x23[] = VectorAlgebra.subtract(x[3],x[2]);

        double rotax1[]=VectorAlgebra.cross(x01,x12);//First peptide and sidechain
        
        rm[0] = new RotationMatrix(rotax1[0],rotax1[1],rotax1[2],paramVal,false);
        
        tr[0]=VectorAlgebra.subtract( rm[0].rotateVector(x01), x01 );

        double y1[]=VectorAlgebra.add(x[1],tr[0]);//Last peptide and sidechain 2.  y1=new 2nd CA coordinates
        double y13[]=VectorAlgebra.subtract(x[3],y1);
        double a=VectorAlgebra.norm(x23);
        double b=VectorAlgebra.norm(y13);
        double d=VectorAlgebra.norm(x12);
        double beta=(double)(Math.acos(VectorAlgebra.dot(x23,y13)/(a*b)) - Math.acos((b*b+a*a-d*d)/(2*a*b)));//Second rotation angle
        double rotax2[]=VectorAlgebra.cross(x23,y13);
        
        rm[2] = new RotationMatrix(rotax2[0], rotax2[1], rotax2[2], beta, true);
        tr[1] = VectorAlgebra.subtract( x23, rm[2].rotateVector(x23) );
        

        double y2[]=VectorAlgebra.add(x[2],tr[1]);//Calculating middle peptide rotation, rm[1].  y2=new 3rd CA coordinates
        double y12[]=VectorAlgebra.subtract(y2,y1);
        double theta=Protractor.getAngleRadians(x12,y12);
        double srotax[]=VectorAlgebra.cross(x12,y12);//Rotation axis to superimpose

        if( VectorAlgebra.norm(srotax) == 0 )//This will happen if x12, y12 are the same.  No feasible shear will have x12 and y12 in opposite directions
            rm[1] = RotationMatrix.identity();
        else
            rm[1] = new RotationMatrix(srotax[0],srotax[1],srotax[2],theta,true);//Rotate middle peptide to superimpose Calphas
        double u[]=rm[1].rotateVector(w);
        
        
        //Additional rotation to correct carbonyl orientation
        double vhat[]=VectorAlgebra.scale(y12,1/VectorAlgebra.norm(y12));//unit vector along line from 2nd to 3rd Calpha
        double alpha=(double)Math.atan2( VectorAlgebra.dot(VectorAlgebra.cross(u,w),vhat), VectorAlgebra.dot(w,u)-VectorAlgebra.dot(w,vhat)*VectorAlgebra.dot(u,vhat));
        RotationMatrix rmalpha = new RotationMatrix(vhat[0],vhat[1],vhat[2],alpha,true);
        
        rm[1] = rmalpha.multiply(rm[1]);        
        
        //OK now generate the full motions
        RigidBodyMotion[] ans = new RigidBodyMotion[3];
        ans[0] = new RigidBodyMotion(x[0], rm[0], x[0]);//rotate about first CA
        //the other two motions have a translation
        for(int m=1; m<3; m++){
            //Rotate about x[m], then translate by tr[m-1]
            double center2[] = VectorAlgebra.add(x[m], tr[m-1]);
            ans[m] = new RigidBodyMotion(x[m], rm[m], center2);
        }
        
        return ans;
    }
    
    
    @Override
    public Perturbation copyForNewMolecule(Molecule mol, PerturbationBlock block){
        Shear s = new Shear(Residue.equivalentInMolec(resDirectlyAffected, mol));
        s.curParamVal = curParamVal;
        s.indexInBlock = indexInBlock;
        s.block = block;
        return s;
    }
    
}
