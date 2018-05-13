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

