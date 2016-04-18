package edu.duke.cs.osprey.control;

import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.minimization.CCDMinimizer;
import edu.duke.cs.osprey.minimization.MolecEObjFunction;

public class EnergyCalculator {

	public void run(ConfigFileParser cfp) {
		
		SearchProblem search = cfp.getSearchProblem();
		
		System.out.println();
		
		double energy;
		if (cfp.getParams().getBool("doMinimize")) {
			
			// run CCD over the continuous degrees of freedom
			// (on existing structure, not any RCs)
			// NOTE: you probably want cfg option standardizeConformation=false for this
			System.out.println("Building energy function...");
			MolecEObjFunction objFunc = new MolecEObjFunction(search.fullConfE, search.confSpace);
			System.out.println(String.format("Minimizing %d degrees of freedom...", objFunc.getNumDOFs()));
			DoubleMatrix1D optDOFVals = new CCDMinimizer(objFunc, false).minimize();
			
			System.out.println("Calculating energy...");
			energy = objFunc.getValue(optDOFVals);
			
		} else {
			
			// just evaluate the energy function
			System.out.println("Calculating energy...");
			energy = search.fullConfE.getEnergy();
		}
		
		System.out.println(String.format("Energy: %f kcal/mol\n", energy));
	}
}
