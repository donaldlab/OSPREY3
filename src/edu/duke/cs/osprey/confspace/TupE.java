/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.confspace;
/**
 *
 * Like ConfPair but with an (RCTuple,energy) pair
 * 
 * @author mhall44
 */
public class TupE implements Comparable<TupE> {
    
    public RCTuple tup;
    public double E;

    public TupE(RCTuple tup, double E) {
        this.tup = tup;
        this.E = E;
    }
    
    
	
    @Override
    public int compareTo(TupE o) {
        return Double.compare(E,o.E);
    }
    
}
