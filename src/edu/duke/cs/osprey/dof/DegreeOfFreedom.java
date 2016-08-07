/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.dof;

import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.Residue;
import java.io.Serializable;

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
    public DegreeOfFreedom copy() {
        // TODO: once all subclasses implement this, change to abstract method
        throw new UnsupportedOperationException();
    }
    public void setMolecule(Molecule val) {
        // TODO: once all subclasses implement this, change to abstract method
        throw new UnsupportedOperationException();
    }
}
