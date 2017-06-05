package edu.duke.cs.osprey.multistatekstar;

import java.util.StringTokenizer;

import edu.duke.cs.osprey.control.ParamSet;

public class ResidueOrderFactory {

	public static ResidueOrder getResidueOrder(ParamSet msParams, MSSearchProblem[][] objFcnSearch) {

		String val = msParams.getValue("RESIDUEORDER");
		StringTokenizer st = new StringTokenizer(val);
		int numTokens = st.countTokens();
		val = numTokens == 1 ? val : st.nextToken().trim();

		val = val.toLowerCase();	
		
		if(val.equals("staticsequential"))
			return new ResidueOrderStaticSequential();
		
		else if(val.equals("staticwtdist"))
			return new ResidueOrderWTDistance(objFcnSearch);
		
		else if(val.equals("staticmindomp"))
			return new ResidueOrderStaticMinDomain(objFcnSearch, "product");
		
		else if(val.equals("staticmindoms"))
			return new ResidueOrderStaticMinDomain(objFcnSearch, "sum");
		
		else if(val.equals("dynamicscoref"))
			return new ResidueOrderDynamicScore(objFcnSearch, "fscore");
		
		else if(val.equals("dynamicscoreh"))	
			return new ResidueOrderDynamicScore(objFcnSearch, "hscore");
		
		else if(val.equals("dynamicscored"))
			return new ResidueOrderDynamicScore(objFcnSearch, "discrepancy");
		
		else if(val.contains("dynamicscore") && val.contains("mindom")) {
			
			String errMsg = "ERROR: parameter mut be of type DynamicScore{F/H/D}MinDom{S,P} X Y";
			
			if(numTokens != 3) throw new RuntimeException(errMsg);
			double dCoeff = Double.valueOf(st.nextToken());
			double mCoeff = Double.valueOf(st.nextToken());
			
			ResidueOrderDynamicScore dynamic = null;
			ResidueOrderStaticMinDomain mindom = null;
			
			if(val.contains("dynamicscoref")) dynamic = new ResidueOrderDynamicScore(objFcnSearch, "fscore");
			else if(val.contains("dynamicscoreh")) dynamic = new ResidueOrderDynamicScore(objFcnSearch, "hscore");
			else if(val.contains("dynamicscored")) dynamic = new ResidueOrderDynamicScore(objFcnSearch, "discrepancy");
			else throw new RuntimeException(errMsg);
			
			if(val.contains("mindoms")) mindom = new ResidueOrderStaticMinDomain(objFcnSearch, "sum");
			else if(val.contains("mindomp")) mindom = new ResidueOrderStaticMinDomain(objFcnSearch, "product");
			else throw new RuntimeException(errMsg);
			
			return new ResidueOrderDynamicScoreMinDom(objFcnSearch, dynamic, dCoeff, mindom, mCoeff);
		}
		
		else
			throw new UnsupportedOperationException("ERROR: unsupported residue order type: "+val);
	}

}
