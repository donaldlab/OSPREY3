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

package edu.duke.cs.osprey.dof;

import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.Residue;
import java.io.Serializable;
import java.util.HashMap;
import java.util.List;

/**
 *
 * @author mhall44
 */
public abstract class DegreeOfFreedom implements Serializable {
    //This class represents a conformational degree of freedom
    //It is used to abstract out flexibility: for purposes of all the conformational search machinery,
    //the conformation can be represented as a vector of degree-of-freedom values
    //then implementations of this class determine what that vector means conformationally
    //by applying appropriate changes to molec. 
    
    //Let's say molecule, once created, can only be changed by DegreeOfFreedom.apply!
    
    private static final long serialVersionUID = 3348198591978172994L;
    
    double curVal;
    
    public abstract void apply(double paramVal);//apply the given parameter value
    //(some degrees of freedom may have convenience methods to call this, e.g. mutation called by aa type)
    
    public double getCurVal() { return curVal; }
    
    
    //If this DegreeOfFreedom moves only a single residue, return that residue
    //Otherwise return null
    //Used in setting up partial energy functions (see MultiTermEnergyFunction)
    public Residue getResidue() { return null; }
    
    // enables parallel molecule manipulation without data races
    // these two methods are only implemented for perturbations that aren't part of a block
    // (DOFBlock.copyForNewMolecule handle these operations in that case)
    public DegreeOfFreedom copy() {
        throw new UnsupportedOperationException("unsupported by " + getClass().getName());
    }
    public void setMolecule(Molecule val) {
        throw new UnsupportedOperationException("unsupported by " + getClass().getName());
    }
    
    public abstract DOFBlock getBlock();//return the DOF block for a DOF (return null if none)

    
    public abstract String getName();//make a name for this DOF
    //should suffice to uniquely identify the DOF among DOFs that appear in a given system
    //but should be the same for equivalent DOFs in different copies of a system
    
    public static HashMap<String,Integer> nameToIndexMap(List<DegreeOfFreedom> dofs){
        //map names to indices in the list dofs
        HashMap<String,Integer> ans = new HashMap<>();
        for(int index=0; index<dofs.size(); index++)
            ans.put(dofs.get(index).getName(), index);
        return ans;
    }
}
