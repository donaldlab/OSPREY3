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
public class TupE implements Comparable {
    
    RCTuple tup;
    double E;

    public TupE(RCTuple tup, double E) {
        this.tup = tup;
        this.E = E;
    }
    
    
	
    @Override
    public int compareTo(Object o) throws ClassCastException {
            // TODO Auto-generated method stub
            if(!(o instanceof TupE))
                    throw new ClassCastException("Another tupE was expected.");
            double otherE = ((TupE) o).E;
            //NOTE THIS IS NORMAL ORDERING FOR ENERGIES. BACKWARDS FROM CONFPAIR
            return Double.valueOf(E).compareTo(otherE);
    }
    
}
