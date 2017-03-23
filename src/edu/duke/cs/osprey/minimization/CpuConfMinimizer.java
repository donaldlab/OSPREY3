package edu.duke.cs.osprey.minimization;

import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.energy.EnergyFunctionGenerator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldInteractions;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.tools.Factory;

public class CpuConfMinimizer extends ConfMinimizer {
	
	public static class Builder {
		
		public final ForcefieldParams ffparams;
		public final Factory<ForcefieldInteractions,Molecule> interactions;
		public final ConfSpace confSpace;
		
		public int numThreads;
		Factory<Minimizer,MoleculeModifierAndScorer> minimizers;
		
		public Builder(ForcefieldParams ffparams, Factory<ForcefieldInteractions,Molecule> interactions, ConfSpace confSpace) {
			
			this.ffparams = ffparams;
			this.interactions = interactions;
			this.confSpace = confSpace;
			
			numThreads = 1;
			minimizers = (mof) -> new SimpleCCDMinimizer(mof);
		}
		
		public Builder setNumThreads(int val) {
			numThreads = val;
			return this;
		}
		
		public Builder setMinimizers(Factory<Minimizer,MoleculeModifierAndScorer> val) {
			minimizers = val;
			return this;
		}
		
		public CpuConfMinimizer build() {
			return new CpuConfMinimizer(numThreads, ffparams, interactions, confSpace, minimizers);
		}
	}
	
	public CpuConfMinimizer(int numThreads, ForcefieldParams ffparams, Factory<ForcefieldInteractions,Molecule> interactions, ConfSpace confSpace, Factory<? extends Minimizer,MoleculeModifierAndScorer> minimizers) {
		
		// make the energy function factory
		EnergyFunctionGenerator egen = new EnergyFunctionGenerator(ffparams);
		Factory<? extends EnergyFunction,Molecule> efuncs = new Factory<EnergyFunction,Molecule>() {
			@Override
			public EnergyFunction make(Molecule mol) {
				return egen.interactionEnergy(interactions.make(mol));
			}
		};
		
		init(numThreads, efuncs, minimizers, confSpace);
	}
}
