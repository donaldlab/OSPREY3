/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.confspace;

import edu.duke.cs.osprey.minimization.ObjectiveFunction.DofBounds;

/**
 *
 * A parametric molecule with bounds on its degrees of freedom
 * SimpleConfSpace generates a ParametricMolecule for a particular conf (RC list),
 * with the residue types and starting conformation implied by that conf
 * It is helpful to bundle it with the conf's DOF bounds, especially for DEEPer and CATS
 * 
 * @author mhall44
 */
public class BoundedParametricMolecule {
    
    public final ParametricMolecule pmol;
    public final DofBounds dofBounds;

    public BoundedParametricMolecule(ParametricMolecule pmol, DofBounds dofBounds) {
        this.pmol = pmol;
        this.dofBounds = dofBounds;
    }
    
    
    
}
