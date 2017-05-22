package edu.duke.cs.osprey.multistatekstar;

import java.util.StringTokenizer;

import edu.duke.cs.osprey.control.ParamSet;

public class ResidueOrderFactory {
	
	public static ResidueOrder getResidueOrder(ParamSet msParams, MSSearchProblem[][] objFcnSearch) {
		
		String val = msParams.getValue("RESIDUEORDER");
		StringTokenizer st = new StringTokenizer(val);
		int numTokens = st.countTokens();
		val = numTokens == 1 ? val : st.nextToken().trim();
		
		switch(val.toLowerCase()) {
		case "staticsequential":
			return new ResidueOrderStaticSequential();
		case "staticmindomp":
			return new ResidueOrderStaticMinDomain(objFcnSearch, true);
		case "staticmindoms":
			return new ResidueOrderStaticMinDomain(objFcnSearch, false);
		case "dynamicscoref":
			return new ResidueOrderDynamicScore(objFcnSearch, "fscore");
		case "dynamicscoreh":
			return new ResidueOrderDynamicScore(objFcnSearch, "hscore");
		case "dynamicscored":
			return new ResidueOrderDynamicScore(objFcnSearch, "discrepancy");
		case "dynamicscoremindom":
			if(numTokens != 3) throw new RuntimeException("ERROR: "
					+ "DynamicScoreMinDom needs coefficients for dynamic and mindom");
			double dCoeff = Double.valueOf(st.nextToken());
			double mCoeff = Double.valueOf(st.nextToken());
			return new ResidueOrderDynamicScoreMinDom(objFcnSearch, dCoeff, mCoeff);
		default:
			throw new UnsupportedOperationException("ERROR: unsupported residue order type: "+val);
		}
	}
	
}
