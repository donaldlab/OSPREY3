/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.tools;

import static edu.duke.cs.osprey.tools.VectorAlgebra.cross;
import static edu.duke.cs.osprey.tools.VectorAlgebra.dot;
import static edu.duke.cs.osprey.tools.VectorAlgebra.norm;
import static edu.duke.cs.osprey.tools.VectorAlgebra.normsq;
import static edu.duke.cs.osprey.tools.VectorAlgebra.perpendicularComponent;
import static edu.duke.cs.osprey.tools.VectorAlgebra.subtract;
import edu.duke.cs.osprey.structure.Residue;

/**
 *
 * @author mhall44
 */
public class Protractor {
    //Measure angles
    
    
    
    public static double getAngleRadians(double vec1[], double vec2[]){//Get the angle, in radians, between two vectors

        double costh = dot(vec1,vec2) / ( norm(vec1) * norm(vec2) );
        if( costh > 1)//It might be slightly over due to numerical error...this means the angle is basically 0
            costh = 1;
        else if(costh < -1 )
            costh = -1;

        return (double) Math.acos( costh );
    }

    public static double getAngleRadians(double A[], double B[], double C[]){//Get the angle ABC
        double BA[] = subtract(A,B);
        double BC[] = subtract(C,B);
        return getAngleRadians(BA,BC);
    }


    public static double getLongitudeRadians(double pt[], double northPole[], double london[]){
        //get the longitude of pt, given the pole direction and something at 0 longitude
        //(all vectors of arbitrary length, relative to origin/center of earth)
        //returns with range (-pi,pi) (in radians)
        double[] ptComp = perpendicularComponent(pt,northPole);

        if(normsq(ptComp)<1e-14)//we're at a pole...define this as longitude 0
            return 0;

        double[] londonComp = perpendicularComponent(london,northPole);
        double ans = getAngleRadians(ptComp,londonComp);
        if( dot(cross(london,pt),northPole) < 0 )
            return -ans;
        else
            return ans;
    }
    
    
    
    public static double measureDihedral(double[][] coords){
        //given 3D coords for four atoms, return their dihedral (standard sign convention; in degrees)
    
	// This was not written by me, but I have checked it
	// If all 4 atoms lie in a plane and the first and fourth
	//  atoms are trans then 180 is returned, if they are cis
	//  then 0 is returned.
	// Between these two extremes, the angle of the right
	//  handed rotation where the axis is the vector from
	//  atom2 to atom3 (ie thumb points to atom3) is returned.
	// The returned angle is between -180-epsilon .. +180
        
        double xij, yij, zij;
        double xkj, ykj, zkj;
        double xkl, ykl, zkl;
        double dx, dy, dz;
        double gx, gy, gz;
        double bi, bk;
        double ct, d, ap, app, bibk;

        xij = coords[0][0] - coords[1][0];
        yij = coords[0][1] - coords[1][1];
        zij = coords[0][2] - coords[1][2];
        xkj = coords[2][0] - coords[1][0];
        ykj = coords[2][1] - coords[1][1];
        zkj = coords[2][2] - coords[1][2];
        xkl = coords[2][0] - coords[3][0];
        ykl = coords[2][1] - coords[3][1];
        zkl = coords[2][2] - coords[3][2];

                    // d = ij cross kj
                    // g = kl cross kj
        dx = yij * zkj - zij * ykj;
        dy = zij * xkj - xij * zkj;
        dz = xij * ykj - yij * xkj;
        gx = zkj * ykl - ykj * zkl;
        gy = xkj * zkl - zkj * xkl;
        gz = ykj * xkl - xkj * ykl;

        bi = dx * dx + dy * dy + dz * dz;  // magnitude of d
        bk = gx * gx + gy * gy + gz * gz;  // magnitude of g
        ct = dx * gx + dy * gy + dz * gz;  // d dot g
                    bibk = bi * bk;
                    if (bibk < 1.0e-6)	
                            return 0;
        ct = ct / Math.sqrt(bibk);
        if(ct < -1.0)
          ct = -1.0;
        else if(ct > 1.0)
          ct = 1.0;

        ap = Math.acos(ct);
        d  = xkj*(dz*gy-dy*gz) + ykj*(dx*gz-dz*gx) + zkj*(dy*gx-dx*gy);
        if(d < 0.0)
          ap = -ap;
        ap = Math.PI - ap;
        app = 180.0 * ap / Math.PI;
        if(app > 180.0)
          app = app - 360.0;
        return(app);
    }
    
    
    
    
        //return (sin(theta/2),cos(theta/2))
    //meant to be quick, without inverse trig
    public static double[] getHalfAngleSinCos(double sinTheta, double cosTheta){
        //first get the right absolute values
        double s2 = (1-cosTheta)/2;
        double c2 = (1+cosTheta)/2;
        
        //assuming negative values of these (cos outside [-1,1]) are just slight numerical errors
        if(s2<0)
            s2=0;
        if(c2<0)
            c2=0;
        
        double ans[] = new double[] {Math.sqrt(s2),Math.sqrt(c2)};
        //now flip signs as needed
        //we are assuming theta is between 0 and 360
        //though this is intended for use in making rotation matrices via quaternions,
        //so actually (sin,cos) and (-sin,-cos) are equivalent
        if(sinTheta<0)//theta/2 in second quadrant.  Else it's in the first.  
            ans[1] *= -1;
        return ans;
    }
    
    /** Returns {phi,psi} for the residue.
     * 
     * @param res: the residue to compute phi and psi
     * @return an array where the first field is phi and the second field is psi
     */
    public static double[] getPhiPsi(Residue res){

        double ans[] = new double[2];

        //Get coordinates of relevant atoms
        //return null (undefined) if can't find one or more atoms
        if(res.indexInMolecule==0 || res.indexInMolecule==res.molec.residues.size()-1)//first or last res
            return null;

        Residue prevRes = res.molec.residues.get(res.indexInMolecule-1);
        Residue nextRes = res.molec.residues.get(res.indexInMolecule+1);


        double[] CLast = prevRes.getCoordsByAtomName("C");
        double[] NCur = res.getCoordsByAtomName("N");
        double[] CACur = res.getCoordsByAtomName("CA");
        double[] CCur = res.getCoordsByAtomName("C");
        double[] NNext = nextRes.getCoordsByAtomName("N");

        if ( CLast==null || NCur==null || CACur==null || CCur==null || NNext==null )
            return null;//atom not found

        ans[0] = Protractor.measureDihedral( new double[][] {CLast,NCur,CACur,CCur} );//phi
        ans[1] = Protractor.measureDihedral( new double[][] {NCur,CACur,CCur,NNext} );//psi

        return ans;
    }
}
