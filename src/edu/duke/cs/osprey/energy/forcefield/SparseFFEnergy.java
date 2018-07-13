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

package edu.duke.cs.osprey.energy.forcefield;

import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.restypes.ResidueTemplate;
import edu.duke.cs.osprey.structure.Atom;
import edu.duke.cs.osprey.structure.Residue;
import java.util.ArrayList;

/**
 * A forcefield energy consisting of only the electrostatics and VDW energies for a selected set of residue pairs
 * Involves only 1 pairs of residues (pairwise interactions) or a single residue (intra-residue energies)
 * Used for SAPE
 * 
 * @author mhall44
 */
public class SparseFFEnergy implements EnergyFunction {
    
    Residue res1, res2;//will be the same for the intra case
    ResidueTemplate templ1=null, templ2=null;//store these so we can tell if there was a mutation
    
    ForcefieldParams ffParams;
    ForcefieldEnergy ffEnergy;
    
    
    public SparseFFEnergy (ArrayList<Atom[]> interactingAtomPairs, ForcefieldParams ffParams){
        
        res1 = interactingAtomPairs.get(0)[0].res;
        res2 = interactingAtomPairs.get(0)[1].res;
        this.ffParams = ffParams;
        
        templ1 = res1.template;
        templ2 = res2.template;
        
        ffEnergy = new ForcefieldEnergy(interactingAtomPairs, ffParams);//sparse ForcefieldEnergy object
    }

    @Override
    public double getEnergy() {
        
        //Mutations mess up our listing of atom pairs, so we cannot use a SparseFFEnergy after mutating any residues involved
        if(res1.template!=templ1 || res2.template!=templ2){
            throw new RuntimeException("ERROR: Cannot calculate a SparseFFEnergy on residues mutated since the SparseFFEnergy was initialized");
        }
        
        if( ! (res1.confProblems.isEmpty() && res2.confProblems.isEmpty()) )
            return Double.POSITIVE_INFINITY;//conformation geometrically impossible
        
        return ffEnergy.calculateTotalEnergy();
    }
    
    
}
