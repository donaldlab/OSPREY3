package edu.duke.cs.osprey.multistatekstar;

import java.io.Serializable;
import edu.duke.cs.osprey.structure.Residue;

@SuppressWarnings("serial")
public class ResidueTermini implements Serializable {

	public int state;
	public int lBound;
	public int uBound;
	
	public ResidueTermini(int state, int begin, int end) {
		if(begin > end) throw new RuntimeException("ERROR: begin: "+begin+" must be <= end: "+end);
		this.state = state;
		this.lBound = begin;
		this.uBound = end;
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
