package edu.duke.cs.osprey.tools;

import java.util.HashMap;
import java.util.Map;
import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.confspace.PositionConfSpace;

public class ResidueIndexMap {

	private Map<Integer, Integer> designIndexToPDBIndex = new HashMap<>();
	private Map<Integer, Integer> PDBIndexToDesignIndex = new HashMap<>();
	
	public static ResidueIndexMap createResidueIndexMap(ConfSpace conformationSpace)
	{
		return new ResidueIndexMap(conformationSpace);
	}
	
	public int designIndexToPDBIndex(int designIndex)
	{
		return designIndexToPDBIndex.get(designIndex);
	}
	
	public int PDBIndexToDesignIndex(int PDBIndex)
	{
		if(!PDBIndexToDesignIndex.containsKey(PDBIndex))
		{
			System.err.println("unrecognized PDBIndex.");
		}
		return PDBIndexToDesignIndex.get(PDBIndex);
	}
	
	private ResidueIndexMap (ConfSpace conformationSpace)
	{
		for(PositionConfSpace positionSpace : conformationSpace.posFlex)
		{
			int PDBIndex = positionSpace.res.getPDBIndex();
			int designIndex = positionSpace.designIndex;
			designIndexToPDBIndex.put(designIndex, PDBIndex);
			PDBIndexToDesignIndex.put(PDBIndex, designIndex);
		}
	}
}
