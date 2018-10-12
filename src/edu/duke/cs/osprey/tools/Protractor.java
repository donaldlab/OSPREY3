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

        return Math.acos( costh );
    }

	public static double getAngleDegrees(double vec1[], double vec2[]) {
    	return Math.toDegrees(getAngleRadians(vec1, vec2));
	}

    
    public static double getAngleDegrees(double A[], double B[], double C[]){
        return Math.toDegrees(getAngleRadians(A,B,C));
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

    public static double measureBondAngle(double[] coords, int[] indices) {
    	return measureBondAngle(coords, indices[0], indices[1], indices[2]);
	}

	public static double measureBondAngle(double[] coords, int a, int b, int c) {
    	return measureBondAngle(coords, a, coords, b, coords, c);
	}

	public static double measureBondAngle(double[][] coords) {
    	return measureBondAngle(coords[0], 0, coords[1], 0, coords[2], 0);
	}

	/** measure the bond angle between three atoms, A, B, C, in degrees */
	public static double measureBondAngle(double[] acoords, int aindex, double[] bcoords, int bindex, double[] ccoords, int cindex) {

		int a3 = aindex*3;
		int b3 = bindex*3;
		int c3 = cindex*3;

    	// ba = a - b
    	double bax = acoords[a3    ] - bcoords[b3    ];
    	double bay = acoords[a3 + 1] - bcoords[b3 + 1];
    	double baz = acoords[a3 + 2] - bcoords[b3 + 2];

    	// bc = c - b
		double bcx = ccoords[c3    ] - bcoords[b3    ];
		double bcy = ccoords[c3 + 1] - bcoords[b3 + 1];
		double bcz = ccoords[c3 + 2] - bcoords[b3 + 2];

		// |ba||bc|cos = ba . bc
		// cos = (ba . bc)/|ba|/|bc|
		double cos = bax*bcx + bay*bcy + baz*bcz;
		cos /= Math.sqrt(bax*bax + bay*bay + baz*baz);
		cos /= Math.sqrt(bcx*bcx + bcy*bcy + bcz*bcz);

		// convert to degrees (clamp to [-1,1] to avoid numerical error)
		if (cos <= -1) {
			return -180;
		} else if (cos >= 1) {
			return 0;
		} else {
			return Math.toDegrees(Math.acos(cos));
		}
	}
    
    public static double measureDihedral(double[] coords, int[] indices) {
    	return measureDihedral(coords, indices[0], indices[1], indices[2], indices[3]);
    }
    
    public static double measureDihedral(double[] coords, int a, int b, int c, int d) {
    	return measureDihedral(coords, a, coords, b, coords, c, coords, d);
    }
    
    public static double measureDihedral(double[][] coords) {
    	return measureDihedral(coords[0], 0, coords[1], 0, coords[2], 0, coords[3], 0);
    }
    
    public static double[] measureDihedralSinCos(double[][] coords) {
    	return measureDihedralSinCos(coords[0], 0, coords[1], 0, coords[2], 0, coords[3], 0);
    }
    
    //given 3D coords for four atoms, return their dihedral (standard sign convention; in degrees)
    public static double measureDihedral(double[] acoords, int aindex, double[] bcoords, int bindex, double[] ccoords, int cindex, double[] dcoords, int dindex) {
    
        double[] sincos = measureDihedralSinCos(acoords, aindex, bcoords, bindex, ccoords, cindex, dcoords, dindex);
        
        // compute theta from sin and cos
        double angleRadians = Math.atan2(sincos[0], sincos[1]);
        return Protractor.normalizeDegrees(Math.toDegrees(angleRadians));
    }

	public static double[] measureDihedralSinCos(double[] coords, int aindex, int bindex, int cindex, int dindex) {
		return measureDihedralSinCos(coords, aindex, coords, bindex, coords, cindex, coords, dindex);
	}

	public static double[] measureDihedralSinCos(double[] acoords, int aindex, double[] bcoords, int bindex, double[] ccoords, int cindex, double[] dcoords, int dindex) {
        //This version returns the {sine,cosine} of the dihedral
        
        // This was not written by me, but I have checked it
        // If all 4 atoms lie in a plane and the first and fourth
        //  atoms are trans then 180 is returned, if they are cis
        //  then 0 is returned.
        // Between these two extremes, the angle of the right
        //  handed rotation where the axis is the vector from
        //  atom2 to atom3 (ie thumb points to atom3) is returned.
        // The returned angle is between -180-epsilon .. +180
        
        // Jeff: I rewrote this, and fixed the numerical stability issues
        
        int a3 = aindex*3;
        int b3 = bindex*3;
        int c3 = cindex*3;
        int d3 = dindex*3;
        
        // let ba = a-b
        double bax = acoords[a3    ] - bcoords[b3    ];
        double bay = acoords[a3 + 1] - bcoords[b3 + 1];
        double baz = acoords[a3 + 2] - bcoords[b3 + 2];
        
        // let bc = c-b
        double bcx = ccoords[c3    ] - bcoords[b3    ];
        double bcy = ccoords[c3 + 1] - bcoords[b3 + 1];
        double bcz = ccoords[c3 + 2] - bcoords[b3 + 2];
        
        // let dc = c-d
        double dcx = ccoords[c3    ] - dcoords[d3    ];
        double dcy = ccoords[c3 + 1] - dcoords[d3 + 1];
        double dcz = ccoords[c3 + 2] - dcoords[d3 + 2];

        // let d = ba cross bc
        double dx = bay*bcz - baz*bcy;
        double dy = baz*bcx - bax*bcz;
        double dz = bax*bcy - bay*bcx;
        
        // let g = bc cross dc
        double gx = bcz*dcy - bcy*dcz;
        double gy = bcx*dcz - bcz*dcx;
        double gz = bcy*dcx - bcx*dcy;
        
        // cos(theta) = (d dot g)/(|d|*|g|)
        double bi = dx*dx + dy*dy + dz*dz;
        double bk = gx*gx + gy*gy + gz*gz;
        double norm = Math.sqrt(bi*bk);
        double cos = (dx*gx + dy*gy + dz*gz)/norm;
        
        /* NOTE: don't do this, this is numerically unstable around 180
        // sin(theta) follows the unit circle
        double sin = 0;
        if (cos > -1 && cos < 1) {
            sin = Math.sqrt(1 - cos*cos);
        }
        */
        
        // let v = g cross d
        double vx = gy*dz - gz*dy;
        double vy = gz*dx - gx*dz;
        double vz = gx*dy - gy*dx;
        
        // sin(theta) = |v|/(|d|*|g|)
        double sin = Math.sqrt(vx*vx + vy*vy + vz*vz)/norm;
        
        // bc dot v
        double bcdotv  = bcx*vx + bcy*vy + bcz*vz;
        
        // adjust the signs to move sin,cos out of the fist quadrant
        // and to match dihedral angle conventions
        cos = -cos;
        if (bcdotv < 0.0) {
            sin = -sin;
        }
        
        return new double[] { sin, cos };
    }
    
    //return (sin(theta/2),cos(theta/2))
    //meant to be quick, without inverse trig
    /**
     * this is numerically unstable, don't use it!
     * (there's a singularity at -180)
     */
    @Deprecated
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
        //return nan (undefined) if can't find one or more atoms
        if(res.indexInMolecule==0 || res.indexInMolecule==res.molec.residues.size()-1)//first or last res
            return new double[] {Double.NaN, Double.NaN};

        Residue prevRes = res.molec.residues.get(res.indexInMolecule-1);
        Residue nextRes = res.molec.residues.get(res.indexInMolecule+1);


        double[] CLast = prevRes.getCoordsByAtomName("C");
        double[] NCur = res.getCoordsByAtomName("N");
        double[] CACur = res.getCoordsByAtomName("CA");
        double[] CCur = res.getCoordsByAtomName("C");
        double[] NNext = nextRes.getCoordsByAtomName("N");

        if ( CLast==null || NCur==null || CACur==null || CCur==null || NNext==null )
            return new double[] {Double.NaN, Double.NaN};//atom not found

        ans[0] = Protractor.measureDihedral( new double[][] {CLast,NCur,CACur,CCur} );//phi
        ans[1] = Protractor.measureDihedral( new double[][] {NCur,CACur,CCur,NNext} );//psi

        return ans;
    }
    
    public static double normalizeDegrees(double angleDegrees) {
    	if (!Double.isFinite(angleDegrees)) {
    		throw new IllegalArgumentException("can't normalize angle: " + angleDegrees);
		}
        while (angleDegrees <= -180) {
            angleDegrees += 360;
        }
        while (angleDegrees > 180) {
            angleDegrees -= 360;
        }
        assert (angleDegrees > -180);
        assert (angleDegrees <= 180);
        return angleDegrees;
    }

    /** returns the smallest distance between a and b */
    public static double getDistDegrees(double a, double b) {
    	return Math.abs(getDeltaDegrees(a, b));
	}

	/** returns the smallest directed distance from a to b */
	public static double getDeltaDegrees(double a, double b) {
		return normalizeDegrees(b - a);
	}

	/**
	 * returns the degrees (relative to 0) in a frame defined by the axis and zero vectors
	 * the axis and zero vectors are assumed to be perpendicular
	 * v is the query vector, assumed to be perpendicular to axis
	 * angles are positive in the ccw direction when looking against the axis
	 */
	public static double getFrameDegrees(double[] axis, double[] zero, double[] v) {

		// the angle magnitude is the easy part
		double degrees = getAngleDegrees(v, zero);

		// but is it positive or negative?
		if (VectorAlgebra.dot(VectorAlgebra.cross(axis, zero), v) < 0) {
			degrees = -degrees;
		}

		return degrees;
	}

	public static double measureBondLength(double[] coords, int aindex, int bindex) {
		return measureBondLength(coords, aindex, coords, bindex);
	}

	public static double measureBondLength(double[] acoords, int aindex, double[] bcoords, int bindex) {

		int a3 = aindex*3;
		int b3 = bindex*3;

		double dx = acoords[a3    ] - bcoords[b3    ];
		double dy = acoords[a3 + 1] - bcoords[b3 + 1];
		double dz = acoords[a3 + 2] - bcoords[b3 + 2];

		return Math.sqrt(dx*dx + dy*dy + dz*dz);
	}


	/**
	 * eg the bond geometry around a Ca atom:
	 * a = N
	 * b = Ca
	 * c = C
	 * d = Cb
	 *
	 * center is b
	 * inPlaneAxis is halfway between ba and bc, negated and normalized
	 * outOfPlaneAxis is perpendicular to the plane of a-b-c, normalized
 	 */
	public static class TetrahedralGeometry {

		public static class Angles {

			public final double inPlaneDegrees;
			public final double outOfPlaneDegrees;

			public Angles(double inPlaneDegrees, double outOfPlaneDegrees) {
				this.inPlaneDegrees = inPlaneDegrees;
				this.outOfPlaneDegrees = outOfPlaneDegrees;
			}
		}


		public double[] ba = new double[3];
		public double[] bc = new double[3];
		public double[] bd = new double[3];

		public double[] center = new double[3];
		public double[] inPlaneAxis = new double[3];
		public double[] outOfPlaneAxis = new double[3];

		public double[] baflat = new double[3];
		public double[] bcflat = new double[3];
		public double[] bdflat = new double[3];

		public double inPlaneDegrees = Double.NaN;
		public double outOfPlaneDegrees = Double.NaN;

		public void update(double[] coords, int aindex, int bindex, int cindex, int dindex) {
			update(coords, aindex, coords, bindex, coords, cindex, coords, dindex);
		}

		public void update(double[] acoords, int aindex, double[] bcoords, int bindex, double[] ccoords, int cindex, double[] dcoords, int dindex) {

			int a3 = aindex*3;
			int b3 = bindex*3;
			int c3 = cindex*3;
			int d3 = dindex*3;

			// center is just b
			center[0] = bcoords[b3    ];
			center[1] = bcoords[b3 + 1];
			center[2] = bcoords[b3 + 2];

			// ba = a - b
			ba[0] = acoords[a3    ] - center[0];
			ba[1] = acoords[a3 + 1] - center[1];
			ba[2] = acoords[a3 + 2] - center[2];

			// bc = c - b
			bc[0] = ccoords[c3    ] - center[0];
			bc[1] = ccoords[c3 + 1] - center[1];
			bc[2] = ccoords[c3 + 2] - center[2];

			// bd = d - b
			bd[0] = dcoords[d3    ] - center[0];
			bd[1] = dcoords[d3 + 1] - center[1];
			bd[2] = dcoords[d3 + 2] - center[2];

			// the out-of-plane axis is ba x bc
			VectorAlgebra.copy(VectorAlgebra.cross(ba, bc), outOfPlaneAxis);
			VectorAlgebra.normalizeInPlace(outOfPlaneAxis);

			VectorAlgebra.copy(VectorAlgebra.perpendicularComponent(ba, outOfPlaneAxis), baflat);
			VectorAlgebra.copy(VectorAlgebra.perpendicularComponent(bc, outOfPlaneAxis), bcflat);
			VectorAlgebra.copy(VectorAlgebra.perpendicularComponent(bd, outOfPlaneAxis), bdflat);

			// the in-plane axis is halfway between baflat and bcflat, but negated
			VectorAlgebra.copy(VectorAlgebra.normalize(baflat), inPlaneAxis);
			VectorAlgebra.addInPlace(inPlaneAxis, VectorAlgebra.normalize(bcflat));
			VectorAlgebra.scaleInPlace(inPlaneAxis, -0.5);
			VectorAlgebra.normalizeInPlace(inPlaneAxis);

			// update tetrahedral angles
			inPlaneDegrees = Protractor.getFrameDegrees(outOfPlaneAxis, inPlaneAxis, bdflat);
			outOfPlaneDegrees = Protractor.getAngleDegrees(bd, outOfPlaneAxis);

			// map out-of-plane angle from [0,180] to [-90,90]
			outOfPlaneDegrees -= 90;
		}

		public Angles getAngles() {
			return new Angles(inPlaneDegrees, outOfPlaneDegrees);
		}
	}
}
