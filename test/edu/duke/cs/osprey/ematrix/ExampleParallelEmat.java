package edu.duke.cs.osprey.ematrix;

import java.util.List;

import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.structure.Residue;

public class ExampleParallelEmat {
	
	@SuppressWarnings("unused")
	private static EnergyMatrix calcEmat(int numThreads, ForcefieldParams ffparams, ConfSpace confSpace, List<Residue> shellResidues) {
		
		// make the energy matrix calculator
		SimpleEnergyMatrixCalculator.Cpu ematcalc = new SimpleEnergyMatrixCalculator.Cpu(numThreads, ffparams, confSpace, shellResidues);
		// NOTE: the super tiny energy matrix minimizations run slowly on the GPU, so just use the CPU for now
		
		// calculate the emat
		EnergyMatrix emat = ematcalc.calcEnergyMatrix();
		
		// don't forget to cleanup so we don't leave extra threads hanging around
		ematcalc.cleanup();
		
		return emat;
	}
}
