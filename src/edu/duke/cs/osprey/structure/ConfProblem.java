/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.structure;

import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import java.io.Serializable;

/**
 *
 * Some conformational degrees of freedom (loop closure adjustment, proline pucker)
 * have parameter values that, for certain starting conformations, are 
 * geometrically impossible.  This object indicates a broken conformation resulting form this.
 * It is applied by the DegreeOfFreedom introducing the problem,
 * and can only be removed by that DegreeOfFreedom (when applied in a geometrically possible way).
 * It will be treated in a forcefield as causing infinite energy for 
 * interactions involving the broken residue
 * 
 * @author mhall44
 */
public class ConfProblem implements Serializable {
    
    DegreeOfFreedom problemDOF;
    Residue brokenResidue;

    
    public ConfProblem(DegreeOfFreedom problemDOF, Residue brokenResidue) {
        //Generate the ConfProblem and attach it to its residue
        this.problemDOF = problemDOF;
        this.brokenResidue = brokenResidue;
        brokenResidue.confProblems.add(this);
    }
    
    
    public void removeFromRes(){
        //the problem is resolved...remove it from its residue
        brokenResidue.confProblems.remove(this);
    }
    
    public Residue getBrokenResidue(){
        return brokenResidue;
    }
}
