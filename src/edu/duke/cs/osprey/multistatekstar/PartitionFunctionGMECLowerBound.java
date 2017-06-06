package edu.duke.cs.osprey.multistatekstar;

import edu.duke.cs.osprey.control.ConfEnergyCalculator.Async;
import edu.duke.cs.osprey.control.ConfSearchFactory;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.pruning.PruningMatrix;

public class PartitionFunctionGMECLowerBound extends PartitionFunctionGMECUpperBound {

	public PartitionFunctionGMECLowerBound(
			EnergyMatrix emat, 
			PruningMatrix pmat, 
			PruningMatrix invmat,
			ConfSearchFactory confSearchFactory, 
			Async ecalc
			) {
		super(emat, pmat, invmat, confSearchFactory, ecalc);
	}

	@Override
	public void init(double targetEpsilon) {
		super.init(targetEpsilon);
	}
	
	@Override
	public void compute(int maxNumConfs) {
		super.compute(maxNumConfs);
	}
}
