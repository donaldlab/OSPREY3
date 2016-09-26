package edu.duke.cs.osprey.kstar.pfunc.impl;

import java.util.ArrayList;

import edu.duke.cs.osprey.control.ConfEnergyCalculator;
import edu.duke.cs.osprey.control.ConfSearchFactory;
import edu.duke.cs.osprey.control.MinimizingEnergyCalculator;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.kstar.KSConfigFileParser;
import edu.duke.cs.osprey.kstar.KSSearchProblem;
import edu.duke.cs.osprey.kstar.pfunc.PFAbstract;
import edu.duke.cs.osprey.kstar.pfunc.SimplePartitionFunction;
import edu.duke.cs.osprey.pruning.PruningMatrix;

public class PFAdapterSimple extends PFAdapter {
	
	private static final long serialVersionUID = 7010927447911304739L;

	public static final String PfImpl = "simple";
	
	private SimplePartitionFunction pfunc;
	private ConfEnergyCalculator.Async ecalc; 
	
	public PFAdapterSimple(int strand, ArrayList<String> sequence, ArrayList<Integer> absolutePos, String checkPointPath, String reducedSPName, KSConfigFileParser cfp, KSSearchProblem panSP) {
		super(PfImpl, strand, sequence, absolutePos, checkPointPath, reducedSPName, cfp, panSP);
		
		// get search things
		KSSearchProblem search = getReducedSearchProblem();
		EnergyMatrix emat = search.emat;
		PruningMatrix pmat = search.reducedMat; // why not just replace pruneMat in the SearchProblem instance?
		
		ConfSearchFactory confSearchFactory = ConfSearchFactory.Tools.makeFromConfig(search, pmat, cfp);
		ecalc = MinimizingEnergyCalculator.makeFromConfig(search, cfp, 0);
		pfunc = new SimplePartitionFunction(emat, pmat, confSearchFactory, ecalc);
		pfunc.setReportProgress(!PFAbstract.suppressOutput);
		setPartitionFunction(pfunc);
	}
	
	@Override
	public void cleanup() {
		super.cleanup();
		setPartitionFunction(null);
		ecalc.cleanup();
		ecalc = null;
	}
}
