/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.dof;

import edu.duke.cs.osprey.structure.Molecule;

/**
 *
 * @author mhall44
 */
public abstract class DegreeOfFreedom {
    //This class represents a conformational degree of freedom
    //It is used to abstract out flexibility: for purposes of all the conformational search machinery,
    //the conformation can be represented as a vector of degree-of-freedom values
    //then implementations of this class determine what that vector means conformationally
    //by applying appropriate changes to molec. 
    
    //Let's say molecule, once created, can only be changed by DegreeOfFreedom.apply!
    
    
    /*Molecule molec;
    double curVal;*/
    
    public abstract void apply(double paramVal);//apply the given parameter value
    //(some degrees of freedom may have convenience methods to call this, e.g. mutation called by aa type)
    
}
