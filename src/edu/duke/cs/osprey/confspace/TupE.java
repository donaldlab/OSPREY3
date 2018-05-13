/*
** This file is part of OSPREY.
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation, either version 2 of the License, or
** (at your option) any later version.
** 
** OSPREY is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
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

