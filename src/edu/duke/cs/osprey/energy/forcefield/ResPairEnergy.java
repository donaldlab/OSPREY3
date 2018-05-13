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

import edu.duke.cs.osprey.restypes.ResidueTemplate;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.structure.Residue;

/**
 *
 * @author mhall44
 */
public class ResPairEnergy implements EnergyFunction {
    //interaction energy between two residues, as modeled by a forcefield
    
    Residue res1, res2;
    
    ResidueTemplate templ1=null, templ2=null;//store these so we can tell if there was a mutation
    
    ForcefieldParams ffParams;
    ForcefieldEnergy ffEnergy;
    
    
    public ResPairEnergy(Residue res1, Residue res2, ForcefieldParams ffParams){
        this.res1 = res1;
        this.res2 = res2;
        this.ffParams = ffParams;
        
        templ1 = res1.template;
        templ2 = res2.template;
        
        initFFE();
    }

    @Override
    public double getEnergy() {
        
        //we'll need to re-initialize the ffe if our residues have been mutated
        if(res1.template!=templ1 || res2.template!=templ2)
            initFFE();
        
        if( ! (res1.confProblems.isEmpty() && res2.confProblems.isEmpty()) )
            return Double.POSITIVE_INFINITY;//conformation geometrically impossible
        
        return ffEnergy.calculateTotalEnergy();
    }
    
    
    void initFFE(){
        ffEnergy = new ForcefieldEnergy(false,res1.atoms,res2.atoms,ffParams);
        templ1 = res1.template;
        templ2 = res2.template;
    }

    public Residue getRes1() {
        return res1;
    }

    public Residue getRes2() {
        return res2;
    }
    
    
    public ForcefieldParams getFFParams() {
        return ffParams;
    }
    
    public ForcefieldEnergy getFFEnergy() {
    	return ffEnergy;
    }
    
    
}

