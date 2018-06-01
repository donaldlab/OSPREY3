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
