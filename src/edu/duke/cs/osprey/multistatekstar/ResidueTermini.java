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
	public String lBound;
	public String uBound;
	//MH: changing the bounds from int to string to accommodate numbers with chain ID
	//(or Kabat numbering)
	//For backwards compatibility, these strings can omit the chain ID,
	//and comparison will be made based on number or chain ID comparisons,
	//not on sequence in an input molecule

	
	public ResidueTermini(int state, String lBound, String uBound) {
		if(lBound.compareTo(uBound)>0) 
			throw new RuntimeException("ERROR: lBound: "+lBound+" must be <= uBound: "+uBound);
		this.state = state;
		this.lBound = lBound;
		this.uBound = uBound;
	}
	
	public boolean contains(Residue res) {
		return contains(res.getPDBResNumber());
	}
	
	public boolean contains(String pdbResNum) {
		return compareResNums(pdbResNum,lBound)>-1 && compareResNums(pdbResNum,uBound)<1;
	}
		
		
	public static int compareResNums(String resNum1, String resNum2) {
		if(hasChainID(resNum1) && hasChainID(resNum2)) {
			Character chainID1 = resNum1.charAt(0);
			int compareChainID = chainID1.compareTo(resNum2.charAt(0));
			if(compareChainID!=0)
				return compareChainID;
		}
		//either chain ID's not provided for both, or they are the same
		//compare numbers
		int compareNumbers = extractPureNumber(resNum1).compareTo(extractPureNumber(resNum2));
		if(compareNumbers!=0)
			return compareNumbers;
		//if we get here either they're the same, or the difference is the Kabat letter
		//If only one has a letter, 123 < 123A
		Character lastChar1 = resNum1.charAt(resNum1.length()-1);
		Character lastChar2 = resNum2.charAt(resNum2.length()-1);
		if(Character.isAlphabetic(lastChar1)) {
			return Character.isAlphabetic(lastChar2) ? lastChar1.compareTo(lastChar2) : 1;
		}
		else {
			return Character.isAlphabetic(lastChar2) ? -1 : 0;
		}
	}
	
	private static boolean hasChainID(String resNum) {
		return Character.isAlphabetic(resNum.charAt(0));
	}
	
	private static Integer extractPureNumber(String resNum) {
		String fullResNum = resNum;
		if(hasChainID(resNum))
			resNum = resNum.substring(1);
		if(Character.isAlphabetic(resNum.charAt(resNum.length()-1)))//Strip Kabat letter
			resNum = resNum.substring(0, resNum.length()-1);
		try {
			return Integer.parseInt(resNum);
		}
		catch(NumberFormatException e) {
			throw new RuntimeException("ERROR invalid residue number: "+fullResNum);
		}
	}
	
	
	public String[] toStringArray() {
		String[] ans = new String[2];
		ans[0] = lBound;
		ans[1] = uBound;
		return ans;
	}
	
}
