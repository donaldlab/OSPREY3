/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.astar.comets;

import java.io.Serializable;
import java.util.StringTokenizer;

/**
 *
 * Linear multistate energy
 * 
 * @author mhall44
 */
public class LME implements Serializable {
    //this is a function of sequence
    //it is a linear function of the GMEC energies for this sequence in all the states
    //just constTerm + sum_s coeffs_s * GMEC_E_for_state_s

    double[] coeffs;
    double constTerm;

    public LME(String s, int numStates){
        //parse from a string listing coefficients in order
        StringTokenizer st = new StringTokenizer(s);
        if(st.countTokens()!=numStates+1){
            throw new RuntimeException("ERROR: SeqTree has "+numStates+" states but GMECLinFunc "
                    + "specified with "+st.countTokens()+" coefficients: "+s);
        }

        coeffs = new double[numStates];
        for(int state=0; state<numStates; state++)
            coeffs[state] = Double.valueOf( st.nextToken() );

        constTerm = Double.valueOf( st.nextToken() );
    }

    public double[] getCoeffs() {
        return coeffs;
    }

    public double getConstTerm() {
        return constTerm;
    }
    
    
    public double eval(double[] stateGMECVals){
        //Given the optimized conformational energies for each state, evaluate the LME
        
        if(stateGMECVals.length != coeffs.length)
            throw new RuntimeException("ERROR: Wrong number of state GMECs");
        
        double ans = 0;
        
        for(int c=0; c<coeffs.length; c++){
            if(coeffs[c] != 0)//specified explicitly to deal with infinite state GMECs
                ans += coeffs[c] * stateGMECVals[c];
        }
        
        ans += constTerm;
        
        if(Double.isNaN(ans))//both positive and negative infinities --> sequence is impossible
            return Double.POSITIVE_INFINITY;
        
        return ans;
    }
    
    
}
