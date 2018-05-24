/*
** This file is part of OSPREY 3.0
** 
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
** 
** OSPREY relies on grants for its development, and since visibility
** in the scientific literature is essential for our success, we
** ask that users of OSPREY cite our papers. See the CITING_OSPREY
** document in this distribution for more information.
** 
** Contact Info:
**    Bruce Donald
**    Duke University
**    Department of Computer Science
**    Levine Science Research Center (LSRC)
**    Durham
**    NC 27708-0129
**    USA
**    e-mail: www.cs.duke.edu/brd/
** 
** <signature of Bruce Donald>, Mar 1, 2018
** Bruce Donald, Professor of Computer Science
*/

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
