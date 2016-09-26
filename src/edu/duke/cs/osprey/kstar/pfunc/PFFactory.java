package edu.duke.cs.osprey.kstar.pfunc;

import java.util.ArrayList;

import edu.duke.cs.osprey.control.ConfEnergyCalculator;
import edu.duke.cs.osprey.control.ConfSearchFactory;
import edu.duke.cs.osprey.control.MinimizingEnergyCalculator;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.kstar.KSConfigFileParser;
import edu.duke.cs.osprey.kstar.KSSearchProblem;
import edu.duke.cs.osprey.kstar.pfunc.impl.PFAdapter;
import edu.duke.cs.osprey.kstar.pfunc.impl.PFParallel0;
import edu.duke.cs.osprey.kstar.pfunc.impl.PFParallel1;
import edu.duke.cs.osprey.kstar.pfunc.impl.PFParallel2;
import edu.duke.cs.osprey.kstar.pfunc.impl.PFTraditional;
import edu.duke.cs.osprey.kstar.pfunc.impl.PFUB;
import edu.duke.cs.osprey.pruning.PruningMatrix;


/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 * 
 * Chooses the type of partition function implementation, depending upon user
 * supplied value
 */
public class PFFactory {

	public static PFAbstract getPartitionFunction( String implementation, int strand, 
			ArrayList<String> sequence, ArrayList<Integer> absolutePos, String checkPointPath, 
			String searchProblemName, KSConfigFileParser cfp, KSSearchProblem sp ) {

		switch( implementation.toLowerCase() ) {

		case "traditional":
			return new PFTraditional( strand, sequence, absolutePos, checkPointPath, searchProblemName, cfp, sp );
		
		case "ub":
			return new PFUB( strand, sequence, absolutePos, checkPointPath, searchProblemName, cfp, sp );

		case "parallel0":
			return new PFParallel0( strand, sequence, absolutePos, checkPointPath, searchProblemName, cfp, sp );
			
		case "parallel1":
			return new PFParallel1( strand, sequence, absolutePos, checkPointPath, searchProblemName, cfp, sp );

		case "parallel2":
			return new PFParallel2( strand, sequence, absolutePos, checkPointPath, searchProblemName, cfp, sp );
			
		case "simple": {
			
			// make the adapter
			PFAdapter adapter = new PFAdapter("simple", strand, sequence, absolutePos, checkPointPath, searchProblemName, cfp, sp);
			
			// make the pfunc and attach it to the adapter
			KSSearchProblem search = adapter.getReducedSearchProblem();
			EnergyMatrix emat = search.emat;
			PruningMatrix pmat = search.reducedMat; // why not just replace pruneMat in the SearchProblem instance?
			ConfSearchFactory confSearchFactory = ConfSearchFactory.Tools.makeFromConfig(search, pmat, cfp);
			ConfEnergyCalculator.Async ecalc = MinimizingEnergyCalculator.makeFromConfig(search, cfp, 0);
			PartitionFunction pfunc = new SimplePartitionFunction(emat, pmat, confSearchFactory, ecalc);
			adapter.setPartitionFunction(pfunc);
			
			return adapter;
		}
		
		default:
			throw new RuntimeException("ERROR: specified value of parameter kStarPFuncMethod is invalid");

		}
	}

}
