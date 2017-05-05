package edu.duke.cs.osprey.multistatekstar;

import java.io.Serializable;
import java.util.ArrayList;

/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 * 
 */

@SuppressWarnings("serial")
public class MSSearchSettings implements Serializable {
	
	public ArrayList<String> mutRes;//reduced flex res
	public ArrayList<ArrayList<String>> AATypeOptions;//reduced allowed AAs
	public double pruningWindow;
	public double stericThreshold;
	public boolean energyLBs;//use to compute either lower or upper bounds
	
	public MSSearchSettings() {}
	
	public String getFormattedSequence() {
		if(AATypeOptions.size()!=mutRes.size())
			throw new RuntimeException("ERROR: sequence length != number of mutable residues");
		
		StringBuilder sb = new StringBuilder();
		for(int i=0;i<mutRes.size();++i) {
			if(mutRes.get(i).equals("-1")) continue;
			StringBuilder sb0 = new StringBuilder();
			for(String aa : AATypeOptions.get(i)) sb0.append(aa+",");
			String aas = sb0.toString();
			aas = aas.substring(0, aas.length()-1);
			sb.append(aas+"-"+mutRes.get(i)+" ");
		}
		return sb.toString().trim();
	}
}
