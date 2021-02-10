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
