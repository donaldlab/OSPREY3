package edu.duke.cs.osprey.multistatekstar;

import java.io.Serializable;
import edu.duke.cs.osprey.structure.Residue;

/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 * 
 */

@SuppressWarnings("serial")
public class ResidueTermini implements Serializable {

	public int state;
	public int lBound;
	public int uBound;
	
	public ResidueTermini(int state, int lBound, int uBound) {
		if(lBound > uBound) throw new RuntimeException("ERROR: lBound: "+lBound+" must be <= uBound: "+uBound);
		this.state = state;
		this.lBound = lBound;
		this.uBound = uBound;
	}
	
	public boolean contains(Residue res) {
		int pdbResNum = Integer.parseInt(res.getPDBResNumber());
		return contains(pdbResNum);
	}
	
	public boolean contains(int pdbResNum) {
		return pdbResNum >= lBound && pdbResNum <= uBound ? true : false;
	}
	
	public int[] toIntArray() {
		int[] ans = new int[2];
		ans[0] = lBound;
		ans[1] = uBound;
		return ans;
	}
	
	public String[] toStringArray() {
		String[] ans = new String[2];
		ans[0] = String.valueOf(lBound);
		ans[1] = String.valueOf(uBound);
		return ans;
	}
}
