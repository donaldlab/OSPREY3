package edu.duke.cs.osprey.ematrix;

import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.energy.ResInterGen;
import edu.duke.cs.osprey.energy.ResidueInteractions;

public enum EnergyPartition {
	
	// TODO: entropy energies
	
	/** intras and shell on singles, inters on pairs */
	Traditional {
		
		@Override
		public ResidueInteractions makeSingleInters(SimpleConfSpace confSpace, SimpleReferenceEnergies eref, int pos, int rc) {
			
			ResInterGen gen = ResInterGen.of(confSpace)
				.addIntra(pos, 1)
				.addShell(pos);
			
			// put reference energies all on singles
			if (eref != null) {
				gen.add(eref.getOffset(confSpace, pos, rc));
			}
			
			return gen.make();
		}
		
		@Override
		public ResidueInteractions makePairInters(SimpleConfSpace confSpace, SimpleReferenceEnergies eref, int pos1, int rc1, int pos2, int rc2) {
			return ResInterGen.of(confSpace)
				.addInter(pos1, pos2)
				.make();
		}
	},
	
	/** inters on pairs, intras and shell distributed evenly among pairs */
	AllOnPairs {
		
		@Override
		public ResidueInteractions makeSingleInters(SimpleConfSpace confSpace, SimpleReferenceEnergies eref, int pos, int rc) {
			// no energies on singles
			return null;
		}
		
		@Override
		public ResidueInteractions makePairInters(SimpleConfSpace confSpace, SimpleReferenceEnergies eref, int pos1, int rc1, int pos2, int rc2) {
			
			int n = confSpace.positions.size();
			double weight = 1.0/(n - 1);
			
			ResInterGen gen = ResInterGen.of(confSpace)
				.addInter(pos1, pos2)
				.addIntra(pos1, weight)
				.addIntra(pos2, weight)
				.addShell(pos1, weight)
				.addShell(pos2, weight);
			
			// distribute eref evenly among pairs
			if (eref != null) {
				gen.add(eref.getOffset(confSpace, pos1, rc1)*weight);
				gen.add(eref.getOffset(confSpace, pos2, rc2)*weight);
			}
			
			return gen.make();
		}
	};
	
	public abstract ResidueInteractions makeSingleInters(SimpleConfSpace confSpace, SimpleReferenceEnergies eref, int pos, int rc);
	public abstract ResidueInteractions makePairInters(SimpleConfSpace confSpace, SimpleReferenceEnergies eref, int pos1, int rc1, int pos2, int rc2);
}
