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
