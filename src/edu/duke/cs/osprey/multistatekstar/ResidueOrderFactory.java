package edu.duke.cs.osprey.multistatekstar;

import edu.duke.cs.osprey.control.ParamSet;

public class ResidueOrderFactory {
	
	public static ResidueOrder getResidueOrder(ParamSet msParams, MSSearchProblem[][] objFcnSearch) {
		String val = msParams.getValue("RESIDUEORDER");
		switch(val.toLowerCase()) {
		case "staticsequential":
			return new ResidueOrderStaticSequential(objFcnSearch);
		case "staticmindom":
		case "staticobjFunchmean":
		case "dynamicobjfunchmean":
		default:
			throw new UnsupportedOperationException("ERROR: unsupported residue order type: "+val);
		}
	}
	
}
