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
