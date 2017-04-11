package edu.duke.cs.osprey.energy;

import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.SimpleReferenceEnergies;

public interface EnergyPartition {
	
	ResidueInteractions makeSingle(SimpleConfSpace confSpace, SimpleReferenceEnergies eref, int pos, int rc);
	ResidueInteractions makePair(SimpleConfSpace confSpace, SimpleReferenceEnergies eref, int pos1, int rc1, int pos2, int rc2);
	ResidueInteractions makeFragment(SimpleConfSpace confSpace, SimpleReferenceEnergies eref, RCTuple frag);
	
	// TODO: entropy energies
	
	/** intras and shell on singles, inters on pairs */
	public static class Traditional implements EnergyPartition {
		
		@Override
		public ResidueInteractions makeSingle(SimpleConfSpace confSpace, SimpleReferenceEnergies eref, int pos, int rc) {
			
			double offset = 0;
			if (eref != null) {
				offset += eref.getOffset(confSpace, pos, rc);
			}
			
			return ResInterGen.of(confSpace)
				.addIntra(pos, 1, offset)
				.addShell(pos)
				.make();
		}
		
		@Override
		public ResidueInteractions makePair(SimpleConfSpace confSpace, SimpleReferenceEnergies eref, int pos1, int rc1, int pos2, int rc2) {
			return ResInterGen.of(confSpace)
				.addInter(pos1, pos2)
				.make();
		}
		
		@Override
		public ResidueInteractions makeFragment(SimpleConfSpace confSpace, SimpleReferenceEnergies eref, RCTuple frag) {
			return ResInterGen.of(confSpace)
				.addIntras(frag, 1, (int pos, int rc) -> {
					double offset = 0;
					if (eref != null) {
						offset += eref.getOffset(confSpace, pos, rc);
					}
					return offset;
				})
				.addInters(frag)
				.addShell(frag)
				.make();
		}
	}
	
	/** inters on pairs, intras and shell distributed evenly among pairs */
	public static class AllOnPairs implements EnergyPartition {
		
		@Override
		public ResidueInteractions makeSingle(SimpleConfSpace confSpace, SimpleReferenceEnergies eref, int pos, int rc) {
			// no energies on singles
			return ResInterGen.of(confSpace)
				.make();
		}
		
		@Override
		public ResidueInteractions makePair(SimpleConfSpace confSpace, SimpleReferenceEnergies eref, int pos1, int rc1, int pos2, int rc2) {
			
			double weight = calcWeight(confSpace);
			
			double offset1 = 0;
			double offset2 = 0;
			if (eref != null) {
				offset1 += eref.getOffset(confSpace, pos1, rc1);
				offset2 += eref.getOffset(confSpace, pos2, rc2);
			}
			
			return ResInterGen.of(confSpace)
				.addInter(pos1, pos2)
				.addIntra(pos1, weight, offset1*weight)
				.addIntra(pos2, weight, offset2*weight)
				.addShell(pos1, weight, 0)
				.addShell(pos2, weight, 0)
				.make();
		}
		
		@Override
		public ResidueInteractions makeFragment(SimpleConfSpace confSpace, SimpleReferenceEnergies eref, RCTuple frag) {
			
			double weight = calcWeight(confSpace);
			
			return ResInterGen.of(confSpace)
				.addIntras(frag)
				.addInters(frag, 1, (int pos1, int rc1, int pos2, int rc2) -> {
					double offset = 0;
					if (eref != null) {
						offset += eref.getOffset(confSpace, pos1, rc1)*weight;
						offset += eref.getOffset(confSpace, pos2, rc2)*weight;
					}
					return offset;
				})
				.addShell(frag)
				.make();
		}
		
		private double calcWeight(SimpleConfSpace confSpace) {
			return 1.0/(confSpace.positions.size() - 1);
		}
	};
}
