package edu.duke.cs.osprey.energy.forcefield;

import java.util.Arrays;

import edu.duke.cs.osprey.structure.Residue;
import java.io.Serializable;

public class AtomCache implements Serializable {

	// combine coords and charges into one array to be cpu cache friendly
	// layout: x, y, z, charge
	public double[] data;
	public int res1Start;
	public int res2Start;
	public Residue res1;
	public Residue res2;
	
	public AtomCache(Residue res1, Residue res2) {
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
	
	public boolean updateCoords() {
		boolean isChanged = false;
		isChanged |= updateCoords(res1Start, res1.coords);
		isChanged |= updateCoords(res2Start, res2.coords);
		
		/* DEBUG: make sure there's no NaNs left
		for (int i=0; i<data.length; i++) {
			if (Double.isNaN(data[i]) || Double.isInfinite(data[i])) {
				throw new Error("AtomCache Updated missed values\n\t" + Arrays.toString(data));
			}
		}
		*/
		
		return isChanged;
	}
	
	private boolean updateCoords(int startIndex, double coords[]) {
		boolean isChanged = false;
		int ix4 = startIndex - 4;
		int ix3 = -3;
		int n = coords.length/3;
		for (int i=0; i<n; i++) {
			ix4 += 4;
			ix3 += 3;
			isChanged |= update(data, ix4 + 0, coords[ix3 + 0]);
			isChanged |= update(data, ix4 + 1, coords[ix3 + 1]);
			isChanged |= update(data, ix4 + 2, coords[ix3 + 2]);
		}
		return isChanged;
	}
	
	private boolean update(double[] out, int i, double val) {
		if (val == out[i]) {
			return false;
		}
		out[i] = val;
		return true;
	}
}
