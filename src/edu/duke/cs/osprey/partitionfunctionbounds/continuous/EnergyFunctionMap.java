package edu.duke.cs.osprey.partitionfunctionbounds.continuous;

import java.util.HashMap;
import edu.duke.cs.osprey.confspace.SearchProblem;

/**
 * Represents a mapping between points in the CRMF domain to rotamers and energies
 * in the protein voxel space
 * @author Adegoke Ojewole (ao68@duke.edu)
 *
 */

public class EnergyFunctionMap {
	
	SearchProblem searchProblem;
	HashMap<double[],int[]> pointConfMap;
	
	public EnergyFunctionMap(SearchProblem otherSearchProblem, 
			HashMap<double[],int[]> otherPointConfMap) {
		this.searchProblem = otherSearchProblem;
		this.pointConfMap = otherPointConfMap;
	}
	
	double getEnergy(double[] point, boolean single, boolean pairWise) {
		return searchProblem.getEnergy(pointConfMap.get(point), point, single, pairWise);
	}
}
