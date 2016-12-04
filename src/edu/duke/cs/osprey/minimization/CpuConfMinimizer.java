package edu.duke.cs.osprey.minimization;

import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.energy.EnergyFunctionGenerator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldInteractions;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.tools.Factory;

public class CpuConfMinimizer extends SpecializedConfMinimizer {
	
	public CpuConfMinimizer(int numThreads, ForcefieldParams ffparams, Factory<ForcefieldInteractions,Molecule> interactions, ConfSpace confSpace) {
		
		// make the minimizer
		EnergyFunctionGenerator egen = new EnergyFunctionGenerator(ffparams, Double.POSITIVE_INFINITY, false);
		Factory<? extends EnergyFunction,Molecule> efuncs = new Factory<EnergyFunction,Molecule>() {
			@Override
			public EnergyFunction make(Molecule mol) {
				return egen.interactionEnergy(interactions.make(mol));
			}
		};
		Factory<? extends Minimizer,MoleculeModifierAndScorer> minimizers = new Factory<SimpleCCDMinimizer,MoleculeModifierAndScorer>() {
			@Override
			public SimpleCCDMinimizer make(MoleculeModifierAndScorer mof) {
				SimpleCCDMinimizer minimizer = new SimpleCCDMinimizer();
				minimizer.init(mof);
				return minimizer;
			}
		};
		
		init(numThreads, efuncs, minimizers, confSpace);
	}
}
