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

package edu.duke.cs.osprey.bbfree;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.VectorAlgebra;
import java.io.Serializable;

/**
 *
 * This is a linear model of a peptide plane, allowing us to get C', H, and O
 * coordinates as a linear function of the N and CA coordinates (CA1, N, CA2)
 * This model can be used even if the peptide plane is translated and rotated
 * 
 * @author mhall44
 */
public class PepPlaneLinModel implements Serializable {
    
    //OK each of C', H, and O will be expressed as CA1 + 
    //a linear combination of (N-CA1), (CA2-CA1), and (N-CA1) X (CA2-CA1).  Coeffs here
    private DoubleMatrix1D CCoeffs;
    private DoubleMatrix1D HCoeffs;
    private DoubleMatrix1D OCoeffs;
    
    
    public PepPlaneLinModel(Residue res1, Residue res2){
        //construct from the current conformations
        //of the residue involved in this peptide plane
        
        if(res1.template.name.equalsIgnoreCase("PRO") || res2.template.name.equalsIgnoreCase("PRO"))
            throw new RuntimeException("ERROR: CATS doesn't handle proline");
        
        double CA1[] = res1.getCoordsByAtomName("CA");
        double CA2[] = res2.getCoordsByAtomName("CA");
        double N[] = res2.getCoordsByAtomName("N");
        double C[] = res1.getCoordsByAtomName("C");
        double O[] = res1.getCoordsByAtomName("O");
        double H[] = res2.getCoordsByAtomName("H");
        
        DoubleMatrix2D M = makeAxisMatrix(CA1,N,CA2);
        DoubleMatrix2D vec = DoubleFactory2D.dense.make(3,1);
        
        vec.viewColumn(0).assign(VectorAlgebra.subtract(C, CA1));
        CCoeffs = Algebra.DEFAULT.solve(M, vec).viewColumn(0);
        
        vec.viewColumn(0).assign(VectorAlgebra.subtract(O, CA1));
        OCoeffs = Algebra.DEFAULT.solve(M, vec).viewColumn(0);
        
        vec.viewColumn(0).assign(VectorAlgebra.subtract(H, CA1));
        HCoeffs = Algebra.DEFAULT.solve(M, vec).viewColumn(0);
    }
    
    
    private DoubleMatrix2D makeAxisMatrix(double[] CA1, double[] N, double CA2[]){
        //matrix of the axes used to expand our atom coordinates
        double axis1[] = VectorAlgebra.subtract(N, CA1);
        double axis2[] = VectorAlgebra.subtract(CA2, CA1);
        double axis3[] = VectorAlgebra.cross(axis1, axis2);
                
        DoubleMatrix2D M = DoubleFactory2D.dense.make(new double[][] {axis1,axis2,axis3});
        return Algebra.DEFAULT.transpose(M);
    }
    
    //calculating atom coordinates from our main atom coordinates
    double[] calcCCoords(double CA1[], double N[], double CA2[], boolean projIntoPlane){
        if(projIntoPlane)
            return calcCoordsProj(CA1,N,CA2,CCoeffs.viewPart(0,2));//skip out-of-plane axis
        else
            return calcCoords(CA1,N,CA2,CCoeffs);
    }
    
    //O and H won't need to be projected into the plane since we don't use constraints on them
    double[] calcOCoords(double CA1[], double N[], double CA2[]){
        return calcCoords(CA1,N,CA2,OCoeffs);
    }
    
    double[] calcHCoords(double CA1[], double N[], double CA2[]){
        return calcCoords(CA1,N,CA2,HCoeffs);
    }
    
    private double[] calcCoords(double CA1[], double N[], double CA2[], DoubleMatrix1D coeffs){
        DoubleMatrix2D M = makeAxisMatrix(CA1,N,CA2);
        return VectorAlgebra.add( CA1, Algebra.DEFAULT.mult(M, coeffs).toArray() );
    }
    
    //calculate the plane projection of the desired atom's coords.  Coeffs for in-plane axes
    private double[] calcCoordsProj(double CA1[], double N[], double CA2[], DoubleMatrix1D coeffs){
        DoubleMatrix2D M = makeAxisMatrix(CA1,N,CA2).viewPart(0, 0, 3, 2);
        return VectorAlgebra.add( CA1, Algebra.DEFAULT.mult(M, coeffs).toArray() );
    }
    
    
    public double[] getProjCAtomCoeffs(){
        //get coefficients for the plane projection of C
        //as a linear combination of CA1, N, and CA2
        double ans[] = new double[3];
        ans[0] = 1 - CCoeffs.get(0) - CCoeffs.get(1);
        ans[1] = CCoeffs.get(0);
        ans[2] = CCoeffs.get(1);
        return ans;
    }
}
