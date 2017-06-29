package edu.duke.cs.osprey.multistatekstar;

import java.util.StringTokenizer;

import edu.duke.cs.osprey.control.ParamSet;

public class ResidueOrderFactory {

	public static ResidueOrder getResidueOrder(ParamSet msParams, MSSearchProblem[][] objFcnSearch) {

		String val = msParams.getValue("RESIDUEORDER");
		StringTokenizer st = new StringTokenizer(val);
		int numTokens = st.countTokens();
		val = numTokens == 1 ? val : st.nextToken().trim();

		val = val.toLowerCase().trim();	

		if(val.equals("sequential"))
			return new ResidueOrderSequential();

		else if(val.startsWith("mindom") && !val.contains("gmec")) {
			String method = null;
			val = val.replace("mindom", "");
			switch(val) {
			case "p": method = "product"; break;
			case "s": method = "sum"; break;
			case "n": method = "numsplits"; break;
			default: 
				throw new RuntimeException("ERROR: parameter must be of type MinDom{N,P,S}");
			}
			return new ResidueOrderMinDomain(objFcnSearch, method);
		}

		else if(val.startsWith("gmec") && !val.contains("mindom")) {
			val = val.replace("gmec", "");
			if(val.length()==0) 
				return new ResidueOrderGMEC(objFcnSearch);
			else {
				String method = null;
				if(val.equalsIgnoreCase("p")) method = "lowest";
				else if(val.equalsIgnoreCase("ph")) method = "hmean";
				else throw new RuntimeException("ERROR: parameter must be of type GMEC{P,PH}");

				return new ResidueOrderGMECProxy(objFcnSearch, method);
			}
		}

		else if(val.startsWith("gmec") && val.contains("mindom")) {
			String errMsg = "ERROR: parameter must be of type GMEC{P}MinDom{N,P,S} X Y";
			if(numTokens != 3) throw new RuntimeException(errMsg);

			double gCoeff = Double.valueOf(st.nextToken());
			double mCoeff = Double.valueOf(st.nextToken());

			String method = null;
			ResidueOrderGMEC gmec = null;

			if(!val.contains("gmecp")) {
				gmec = new ResidueOrderGMEC(objFcnSearch);
			}
			else {
				if(val.contains("gmecph")) method = "hmean";
				else if(val.contains("gmecp")) method = "lowest";
				else throw new RuntimeException("ERROR: parameter must be of type GMEC{P,PH}");
				gmec = new ResidueOrderGMECProxy(objFcnSearch, method);
			}
			
			if(val.contains("mindoms")) method = "sum";
			else if(val.contains("mindomp")) method = "product";
			else if(val.contains("mindomn")) method = "numsplits";
			else throw new RuntimeException(errMsg);

			ResidueOrderMinDomain mindom = new ResidueOrderMinDomain(objFcnSearch, method);
			return new ResidueOrderMeta(gmec, gCoeff, mindom, mCoeff);
		}

		else throw new UnsupportedOperationException("ERROR: unsupported residue order type: "+val);
	}

}
