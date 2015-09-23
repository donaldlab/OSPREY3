/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.confspace;

/**
 *
 * @author hmn5
 */
public class SuperTupE implements Comparable {

    SuperRCTuple tup;
    double E;

    public SuperTupE(SuperRCTuple tup, double E) {
        this.tup = tup;
        this.E = E;
    }

    @Override
    public int compareTo(Object o) throws ClassCastException {
        // TODO Auto-generated method stub
        if (!(o instanceof SuperTupE)) {
            throw new ClassCastException("Another tupE was expected.");
        }
        double otherE = ((SuperTupE) o).E;
        //NOTE THIS IS NORMAL ORDERING FOR ENERGIES. BACKWARDS FROM CONFPAIR
        return Double.valueOf(E).compareTo(otherE);
    }
}
