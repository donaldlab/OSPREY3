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

package edu.duke.cs.osprey.tupexp;

/**
 *
 * @author mhall44
 */
///////////////////////////////////////////////////////////////////////////////////////////////
//	ConfPair.java
//
//	Version:           2.1 beta
//
//
//	  authors:
// 	  initials    name                 organization                email
//	---------   -----------------    ------------------------    ----------------------------
//     KER        Kyle E. Roberts       Duke University         ker17@duke.edu


// This class is used to store conformations.  It is used in OSPREY to store the top conformations.  
public class ConfPair implements Comparable{
	int[] conf;
	//minE: 0 unMinE: 1
	double[] energy;
	public ConfPair(int[] conformation, double[] e){
		conf = new int[conformation.length];
		for(int i=0; i<conformation.length;i++)
			conf[i] = conformation[i];
		energy = new double[e.length];
		for(int i=0; i<e.length;i++)
			energy[i] = e[i];
		
	}
	
	@Override
	public int compareTo(Object o) throws ClassCastException {
		// TODO Auto-generated method stub
		if(!(o instanceof ConfPair))
			throw new ClassCastException("Another confPair was expected.");
		double otherE = ((ConfPair) o).energy[0];
		if(otherE >= energy[0])
			return 1;
		else
			return -1;
		
	}

}


