/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.minimization;

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
//	GenCoord.java
//
//	Version:           2.1 beta
//
//
//	  authors:
// 	  initials    name                 organization                email
//	 ---------   -----------------    ------------------------    ----------------------------
//	  MAH           Mark A. Hallen	  Duke University               mah43@duke.edu
///////////////////////////////////////////////////////////////////////////////////////////////

import cern.colt.matrix.DoubleMatrix1D;
import java.io.Serializable;


@SuppressWarnings("serial")
public class GenCoord implements Serializable {
    //This class represents a generalized coordinate
    //it facilitates calculation of the value of a generalized coordinate from a set of "regular" coordinates like dihedrals
    //used in the CCDMinimizer

    int type;

    //possible types
    static final int REGULAR=0;//Just a regular coordinate
    static final int SUMSQ=1;//Square root of sum of squares of regular coordinates
    static final int LINCOMB=2;//Linear combination of regular coordinates


    double coeffs[] = null;//coefficients for linear combination

    public GenCoord(){//Regular coordinate by default
        type = REGULAR;
    }

    public GenCoord(double c[]){//Linear-combination coordinate
        type = LINCOMB;
        coeffs = c;
    }

    public double eval( DoubleMatrix1D x, int[] indices ) {
        //The "regular" coordinates to be used have the indicated indices in x

        switch(type){
            case REGULAR:
                return x.get(indices[0]);
            case SUMSQ:
                double sum=0;
                for(int q:indices)
                    sum += Math.pow(x.get(q), 2);
                return Math.sqrt(sum);
            case LINCOMB:
                double ans=0;
                for(int a=0; a<coeffs.length; a++)
                    ans += coeffs[a]*x.get(indices[a]);
                return ans;
            default:
                throw new Error("Unrecognized generalized coordinate type: "+type);
        }
    }

    public double getNearestInRangeDOFVal( double startVal, double min, double max,
            DoubleMatrix1D x, int dof, int[] indices ) {
        //Return the nearest value of DOF #dof in x to startVal that will make this GenCoord (operating on
        //the given indices of x) return a value in the range[min,max]
        //If there is no such value, NaN is returned
        switch(type){
            
            case REGULAR:
                if( startVal < min )
                    return min;
                else if(startVal > max )
                    return max;
                else//startVal is in range already
                    return startVal;

            case SUMSQ:
                double sum=0;
                for(int q:indices)
                    sum += Math.pow(x.get(q), 2);
                double minsq = Math.pow(min, 2);
                if(sum<minsq){//Increase the absolute value of x.get(dof) until it makes sum equal to min^2
                    double absAns = Math.sqrt( Math.pow(min, 2) - (sum-Math.pow(startVal,2)) );
                    if( startVal > 0 )
                        return absAns;
                    else
                        return -absAns;
                }
                double maxsq = Math.pow(max, 2);
                if(sum>maxsq){//Decrease the absolute value of x.get(dof) until it makes sum equal to min^2
                    double absAns = Math.sqrt( Math.pow(max, 2) - (sum-Math.pow(startVal,2)) );
                    if( startVal > 0 )
                        return absAns;
                    else
                        return -absAns;
                }
                else//In range already
                    return startVal;

            case LINCOMB:
                double ans=0;
                int DOFLocalInd = -1;//Index of dof in the "indices" array
                for(int a=0; a<coeffs.length; a++){
                    if( indices[a] == dof ){
                        ans += coeffs[a]*startVal;
                        DOFLocalInd = a;
                    }
                    else
                        ans += coeffs[a]*x.get(indices[a]);
                }
                if( ans < min ){
                    if(coeffs[DOFLocalInd] == 0)
                        return Double.NaN;
                    return startVal + (min-ans)/coeffs[DOFLocalInd];
                }
                else if ( ans > max ){
                    if(coeffs[DOFLocalInd] == 0)
                        return Double.NaN;
                    return startVal + (max-ans)/coeffs[DOFLocalInd];
                }
                else//In range already
                    return startVal;
                
            default:
                throw new Error("Unrecognized generalized coordinate type: "+type);
        }
    }

    public double constrOpt(DoubleMatrix1D constraints[], int indices[], boolean useMin){
         //constraints[0] is DOF minima and constraints[1] is DOF maxima
        //Return the minimum or maximum (according to useMin) of this GenCoord given the
        //constraints on the regular DOFs
        //(indices in the constraint vectors given)
        switch(type){

            case REGULAR:
                if(useMin)
                    return constraints[0].get(indices[0]);
                else
                    return constraints[1].get(indices[0]);

            case SUMSQ:
                double sum=0;
                for(int q:indices){
                    double sq1 = Math.pow( constraints[0].get(q), 2);
                    double sq2 = Math.pow( constraints[1].get(q), 2);
                    if(useMin)
                        sum += Math.min(sq1,sq2);
                    else
                        sum += Math.max(sq1,sq2);
                }
                return Math.sqrt(sum);

            case LINCOMB:
                double ans=0;
                for(int a=0; a<coeffs.length; a++){
                    double val1 = coeffs[a]*constraints[0].get(indices[a]);
                    double val2 = coeffs[a]*constraints[1].get(indices[a]);
                    if(useMin)
                        ans += Math.min(val1,val2);
                    else
                        ans += Math.max(val1,val2);
                }
                return ans;
                
            default:
                throw new Error("Unrecognized generalized coordinate type: "+type);
        }
    }
    

    public double getMeshWidth(double DOFMeshWidths[]){
        //Return a good mesh width for a FuncLBTerm depending on this GenCoord
        //given mesh widths for the DOFs that this GenCoord operates on
        switch(type){

            case REGULAR://Just use the DOF mesh width
                return DOFMeshWidths[0];

            case LINCOMB://Return an average of the DOF mesh widths, weighted by the absolute values of their coefficients
                double numerator = 0, denominator = 0;
                for(int a=0; a<coeffs.length; a++){
                    numerator += Math.abs(coeffs[a])*DOFMeshWidths[a];
                    denominator += Math.abs(coeffs[a]);
                }
                return numerator/denominator;

            default:
                throw new Error("Unrecognized generalized coordinate type for getMeshWidth: "+type);
        }
    }

}
