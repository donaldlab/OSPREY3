package edu.duke.cs.osprey.multistatekstar;

import edu.duke.cs.osprey.control.ParamSet;

public class ResidueOrderFactory {
	
	public static ResidueOrder getResidueOrder(ParamSet msParams, MSSearchProblem[][] objFcnSearch) {
		String val = msParams.getValue("RESIDUEORDER");
		switch(val.toLowerCase()) {
		case "staticsequential":
			return new ResidueOrderStaticSequential();
		case "staticmindomp":
			return new ResidueOrderStaticMinDomain(objFcnSearch, true);
		case "staticmindoms":
			return new ResidueOrderStaticMinDomain(objFcnSearch, false);
		case "dynamicscoref":
			return new ResidueOrderDynamicScore(objFcnSearch, true);
		case "dynamicscoreh":
			return new ResidueOrderDynamicScore(objFcnSearch, false);
		case "dynamicscoremindom":
			return new ResidueOrderDynamicScoreMinDom(objFcnSearch);
		default:
			throw new UnsupportedOperationException("ERROR: unsupported residue order type: "+val);
		}
	}
	
}
