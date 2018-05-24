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

package edu.duke.cs.osprey.energy.forcefield;

import java.util.Arrays;

import edu.duke.cs.osprey.structure.Residue;
import java.io.Serializable;

public class CoordsAndCharges implements Serializable {

	private static final long serialVersionUID = 3611517841418160913L;
	
	// combine coords and charges into one array to be cpu cache friendly
	// layout: x, y, z, charge
	public double[] data;
	public int res1Start;
	public int res2Start;
	public Residue res1;
	public Residue res2;
	
	public CoordsAndCharges(Residue res1, Residue res2) {
		boolean isInternal = res1 == res2;
		if (isInternal) {
			this.data = new double[res1.atoms.size()*4];
			this.res1Start = 0;
			this.res2Start = 0;
		} else {
			this.data = new double[(res1.atoms.size() + res2.atoms.size())*4];
			this.res1Start = 0;
			this.res2Start = res1.atoms.size()*4;
		}
		this.res1 = res1;
		this.res2 = res2;
		
		// intialize data to something bogus so the first update() catches changes correctly
		Arrays.fill(data, Double.NaN);
		
		// write the charges now, since they don't change over time
		updateCharges(res1Start, res1);
		updateCharges(res2Start, res2);
	}
	
	private void updateCharges(int startIndex, Residue res) {
		for (int i=0; i<res.atoms.size(); i++) {
			data[startIndex + i*4 + 3] = res.atoms.get(i).charge;
		}
	}
	
	public void updateCoords() {
		updateCoords(res1Start, res1.coords);
		updateCoords(res2Start, res2.coords);
	}
	
	private void updateCoords(int startIndex, double coords[]) {
		int ix4 = startIndex - 4;
		int ix3 = -3;
		int n = coords.length/3;
		for (int i=0; i<n; i++) {
			ix4 += 4;
			ix3 += 3;
			data[ix4] = coords[ix3];
			data[ix4 + 1] = coords[ix3 + 1];
			data[ix4 + 2] = coords[ix3 + 2];
		}
	}
}
