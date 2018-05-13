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
import edu.duke.cs.osprey.structure.Residue;

/**
 *
 * @author mhall44
 */
public class SingleResEnergy implements EnergyFunction {
    //internal energy as a single residue, as modeled by a forcefield
    //very similar to ResPairEnergy
    
    Residue res;
    
    ResidueTemplate templ;//store this so we can tell if there was a mutation
    
    ForcefieldParams ffParams;
    ForcefieldEnergy ffEnergy;
    
    
    public SingleResEnergy(Residue res, ForcefieldParams ffParams){
        this.res = res;
        this.ffParams = ffParams;
        
        templ = res.template;
        
        initFFE();
    }

    
    @Override
    public double getEnergy() {
        
        //we'll need to re-initialize the ffe if our residues have been mutated
        if(res.template!=templ)
            initFFE();
        
        if( ! res.confProblems.isEmpty() )
            return Double.POSITIVE_INFINITY;//conformation geometrically impossible
        
        return ffEnergy.calculateTotalEnergy();
    }
    
    
    void initFFE(){
        ffEnergy = new ForcefieldEnergy(true,res.atoms,res.atoms,ffParams);
        templ = res.template;
    }

    public Residue getRes() {
        return res;
    }

    public ForcefieldParams getFFParams() {
        return ffParams;
    }
    
    public ForcefieldEnergy getFFEnergy() {
    	return ffEnergy;
    }
    
    
}

