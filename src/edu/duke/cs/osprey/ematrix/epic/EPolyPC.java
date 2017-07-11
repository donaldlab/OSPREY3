package edu.duke.cs.osprey.ematrix.epic;

/*
	This file is part of OSPREY.

	OSPREY Protein Redesign Software Version 2.1 beta
	Copyright (C) 2001-2012 Bruce Donald Lab, Duke University

	OSPREY is free software: you can redistribute it and/or modify
	it under the terms of the GNU Lesser General Public License as
	published by the Free Software Foundation, either version 3 of
	the License, or (at your option) any later version.

	OSPREY is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
	GNU Lesser General Public License for more details.

	You should have received a copy of the GNU Lesser General Public
	License along with this library; if not, see:
	      <http://www.gnu.org/licenses/>.

	There are additional restrictions imposed on the use and distribution
	of this open-source code, including: (A) this header must be included
	in any modification or extension of the code; (B) you are required to
	cite our papers in any publications that use this code. The citation
	for the various different modules of our software, together with a
	complete list of requirements and restrictions are found in the
	document license.pdf enclosed with this distribution.

	Contact Info:
			Bruce Donald
			Duke University
			Department of Computer Science
			Levine Science Research Center (LSRC)
			Durham
			NC 27708-0129
			USA
			e-mail:   www.cs.duke.edu/brd/

	<signature of Bruce Donald>, Mar 1, 2012
	Bruce Donald, Professor of Computer Science
 */

///////////////////////////////////////////////////////////////////////////////////////////////
//	EPolyPC.java
//
//	Version:           2.1 beta
//
//
//	  authors:
// 	  initials    name                 organization                email
//	 ---------   -----------------    ------------------------    ----------------------------
//	  MAH           Mark A. Hallen	  Duke University               mah43@duke.edu
///////////////////////////////////////////////////////////////////////////////////////////////


import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import cern.colt.matrix.linalg.EigenvalueDecomposition;


public class EPolyPC extends EPoly {
    //This is a version of EPoly with a series expanded in the principal components
    //of a quadratic template fit; higher-order terms are based only on the
    //more principal of these components

    int fullOrder;//full terms
    //then we go up to PCOrder in PC-only terms
    int PCOrder;

    boolean isPC[];//which of the (local coordinate) DOFs are PCs
    
    DoubleMatrix2D axisCoeffs;//Coefficients of the new axis vectors (as linear combinations of the DOFs in DOFNums)
    //rows of axisCoeffs are the coefficients for a given vector 
    //(the eigenvectors of the Hessian at the minimum point)
    DoubleMatrix2D invAxisCoeffs;//Inverse of axisCoeffs


    
    
    public EPolyPC( EPoly template, int fullOrder, int PCOrder, double PCFac) {
        //set up the principal component basis and the DOF information
        //coeffs will be added later
        //since we will be using this EPolyPC object when performing the fitting
        
        super(template.numDOFs,template.DOFmax,template.DOFmin,template.center,
                template.minE,null,fullOrder,template.DOFNames);
        
        //get the coordinate transformation into the eigenbasis of the template's Hessian
        //so we can use the high-eigenvalue directions as our principcal components
        DoubleMatrix2D hess = SeriesFitter.getHessian(template.coeffs,numDOFs,false);
        EigenvalueDecomposition edec = new EigenvalueDecomposition(hess);

        DoubleMatrix1D eigVals = edec.getRealEigenvalues();
        DoubleMatrix2D eigVecs = edec.getV();//The columns of eigVec are the eigenvectors

        invAxisCoeffs = eigVecs;//axisCoeffs' inverse is its transpose, since it's orthogonal
        axisCoeffs = Algebra.DEFAULT.transpose(eigVecs);
        

        this.fullOrder = fullOrder;
        this.PCOrder = PCOrder;


        //Now figure out what PC's to use
        //Let's try the criterion to use all PCs
        //whose eigenvalues' absolute values are within a factor of PCFac
        //of the biggest

        double maxEigVal = 0;
        for(int n=0; n<numDOFs; n++)
            maxEigVal = Math.max( Math.abs(eigVals.get(n)), maxEigVal );

        isPC = new boolean[numDOFs];

        for(int n=0; n<numDOFs; n++){
            if( Math.abs(eigVals.get(n)) >= maxEigVal*PCFac )
                isPC[n] = true;
        }
    }


    //evaluate() works as in EPoly
    //except the series is evaluated differently...
    @Override
    double evalSeries(DoubleMatrix1D z){
        //evaluate the actual series
        //(function of relative coordinates)
        DoubleMatrix1D y = toPCBasis(z);
        return SeriesFitter.evalSeries(coeffs, y, numDOFs, false, fullOrder, PCOrder, isPC);
    }
    
    
    //conversion from relative coordinates in usual DOF basis to eigenbasis of template Hessian
    DoubleMatrix1D toPCBasis(DoubleMatrix1D z){
        return axisCoeffs.zMult(z, DoubleFactory1D.dense.make(numDOFs));
    }

}
