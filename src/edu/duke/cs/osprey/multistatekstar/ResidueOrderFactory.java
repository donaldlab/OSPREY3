package edu.duke.cs.osprey.multistatekstar;

import edu.duke.cs.osprey.control.ParamSet;

public class ResidueOrderFactory {
	
	public static ResidueOrder getResidueOrder(ParamSet msParams, MSSearchProblem[][] objFcnSearch) {
		String val = msParams.getValue("RESIDUEORDER");
		switch(val.toLowerCase()) {
		case "staticsequential":
			return new ResidueOrderStaticSequential();
		case "dynamicfscore":
			return new ResidueOrderDynamicScore(objFcnSearch);
		default:
			throw new UnsupportedOperationException("ERROR: unsupported residue order type: "+val);
		}
	}
	
}
