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
