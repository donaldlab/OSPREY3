package edu.duke.cs.osprey.energy;

import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.SimpleReferenceEnergies;
import edu.duke.cs.osprey.restypes.ResidueTemplate;

public enum EnergyPartition {
	
	/** intras and shell on singles, inters on pairs */
	Traditional {
		
		@Override
		public ResidueInteractions makeSingle(SimpleConfSpace confSpace, SimpleReferenceEnergies eref, boolean addResEntropy, int pos, int rc) {
			
			double offset = 0;
			if (eref != null) {
				offset += eref.getOffset(confSpace, pos, rc);
			}
			if (addResEntropy) {
				offset += getResEntropy(confSpace, pos, rc);
			}
			
			return ResInterGen.of(confSpace)
				.addIntra(pos, 1, offset)
				.addShell(pos)
				.make();
		}
		
		@Override
		public ResidueInteractions makePair(SimpleConfSpace confSpace, SimpleReferenceEnergies eref, boolean addResEntropy, int pos1, int rc1, int pos2, int rc2) {
			return ResInterGen.of(confSpace)
				.addInter(pos1, pos2)
				.make();
		}
	},
	
	/** inters on pairs, intras and shell distributed evenly among pairs */
	AllOnPairs {
		
		@Override
		public ResidueInteractions makeSingle(SimpleConfSpace confSpace, SimpleReferenceEnergies eref, boolean addResEntropy, int pos, int rc) {
			
			// only energies on singles if there's exactly one position in the conf space
			// (meaning, there's no pairs to put all the energies onto)
			if (confSpace.positions.size() == 1) {
				return Traditional.makeSingle(confSpace, eref, addResEntropy, pos, rc);
			}
			
			// otherwise, no energies on singles
			return ResInterGen.of(confSpace)
				.make();
		}
		
		@Override
		public ResidueInteractions makePair(SimpleConfSpace confSpace, SimpleReferenceEnergies eref, boolean addResEntropy, int pos1, int rc1, int pos2, int rc2) {
			
			double offset1 = 0;
			double offset2 = 0;
			if (eref != null) {
				offset1 += eref.getOffset(confSpace, pos1, rc1);
				offset2 += eref.getOffset(confSpace, pos2, rc2);
			}
			if (addResEntropy) {
				offset1 += getResEntropy(confSpace, pos1, rc1);
				offset2 += getResEntropy(confSpace, pos2, rc2);
			}
			
			double weight = calcWeight(confSpace);
			
			return ResInterGen.of(confSpace)
				.addInter(pos1, pos2)
				.addIntra(pos1, weight, offset1)
				.addIntra(pos2, weight, offset2)
				.addShell(pos1, weight, 0)
				.addShell(pos2, weight, 0)
				.make();
		}
		
		private double calcWeight(SimpleConfSpace confSpace) {
			return 1.0/(confSpace.positions.size() - 1);
		}
	};
	
	public abstract ResidueInteractions makeSingle(SimpleConfSpace confSpace, SimpleReferenceEnergies eref, boolean addResEntropy, int pos, int rc);
	public abstract ResidueInteractions makePair(SimpleConfSpace confSpace, SimpleReferenceEnergies eref, boolean addResEntropy, int pos1, int rc1, int pos2, int rc2);
	
	public static ResidueInteractions makeFragment(SimpleConfSpace confSpace, SimpleReferenceEnergies eref, boolean addResEntropy, RCTuple frag) {
		return ResInterGen.of(confSpace)
			.addIntras(frag, 1, (int pos, int rc) -> {
				double offset = 0;
				if (eref != null) {
					offset += eref.getOffset(confSpace, pos, rc);
				}
				if (addResEntropy) {
					offset += getResEntropy(confSpace, pos, rc);
				}
				return offset;
			})
			.addInters(frag)
			.addShell(frag)
			.make();
	}
	
	public static double getResEntropy(SimpleConfSpace confSpace, int pos, int rc) {
		ResidueTemplate template = confSpace.positions.get(pos).resConfs.get(rc).template;
		return confSpace.positions.get(pos).strand.templateLib.getResEntropy(template.name);
	}
}
