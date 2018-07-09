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
