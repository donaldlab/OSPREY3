/*
** This file is part of OSPREY.
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation, either version 2 of the License, or
** (at your option) any later version.
** 
** OSPREY is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
*/


package edu.duke.cs.osprey.bbfree;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.dof.DOFBlock;
import edu.duke.cs.osprey.dof.DegreeOfFreedom;

/**
 *
 * This is a DOF that's part of a BBFreeBlock
 * 
 * @author mhall44
 */
public class BBFreeDOF extends DegreeOfFreedom {
    
    DoubleMatrix1D coeffs;//This free DOF is a linear combination of the full DOFs in its block
    BBFreeBlock block;
    int indexInBlock;

    public BBFreeDOF(DoubleMatrix1D coeffs, BBFreeBlock block, int indexInBlock) {
        this.coeffs = coeffs;
        this.block = block;
        this.indexInBlock = indexInBlock;
    }

    @Override
    public void apply(double paramVal) {
        double[] newDOFVals = block.curFreeDOFVals.clone();
        newDOFVals[indexInBlock] = paramVal;
        block.setDOFs(DoubleFactory1D.dense.make(newDOFVals));
    }
    
    
    public double evalAtFullDOFs(DoubleMatrix1D fullDOFVals){
        return fullDOFVals.zDotProduct(coeffs);
    }
    
    
    @Override
    public DOFBlock getBlock(){
        return block;
    }
    
    @Override
    public String getName() {
        return "CATS"+block.residues.get(0).getPDBResNumber()+"."+indexInBlock;
    }
}

