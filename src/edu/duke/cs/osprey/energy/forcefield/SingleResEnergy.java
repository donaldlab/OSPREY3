/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
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
