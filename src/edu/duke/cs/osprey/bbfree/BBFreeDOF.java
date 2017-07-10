/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
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
