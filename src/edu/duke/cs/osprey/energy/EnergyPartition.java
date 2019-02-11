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

package edu.duke.cs.osprey.energy;

import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
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

		@Override
		public ResidueInteractions makeTuple(SimpleConfSpace confSpace, SimpleReferenceEnergies eref, boolean addResEntropy, RCTuple tuple) {
			return ResInterGen.of(confSpace)
				.addInters(tuple)
				.make();
		}

		@Override
		public ResidueInteractions makeTripleCorrection(SimpleConfSpace confSpace, SimpleReferenceEnergies eref, boolean addResEntropy, int pos1, int rc1, int pos2, int rc2, int pos3, int rc3) {
			double weight = 1.0/(confSpace.positions.size() - 2);
			double offset = 0.0;
			return ResInterGen.of(confSpace)
				.addInter(pos1, pos2, weight, offset)
				.addInter(pos1, pos3, weight, offset)
				.addInter(pos2, pos3, weight, offset)
				.make();
		}

		@Override
		public double offsetTripleEnergy(int pos1, int rc1, int pos2, int rc2, int pos3, int rc3, EnergyMatrix emat) {
			double weight = 1.0/(emat.getNumPos() - 2);
			return weight*(
				  emat.getPairwise(pos1, rc1, pos2, rc2)
				+ emat.getPairwise(pos1, rc1, pos3, rc3)
				+ emat.getPairwise(pos2, rc2, pos3, rc3)
			);
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
			
			double weight = 1.0/(confSpace.positions.size() - 1);
			
			return ResInterGen.of(confSpace)
				.addInter(pos1, pos2)
				.addIntra(pos1, weight, offset1)
				.addIntra(pos2, weight, offset2)
				.addShell(pos1, weight, 0)
				.addShell(pos2, weight, 0)
				.make();
		}

		@Override
		public ResidueInteractions makeTuple(SimpleConfSpace confSpace, SimpleReferenceEnergies eref, boolean addResEntropy, RCTuple tuple) {

			// short circuit
			if (tuple.size() == 1) {
				return makeSingle(confSpace, eref, addResEntropy, tuple.pos.get(0), tuple.RCs.get(0));
			}

			double weight = (double)(tuple.size() - 1)/(confSpace.positions.size() - 1);

			return ResInterGen.of(confSpace)
				.addIntras(tuple, weight, (pos, rc) -> {
					double offset = 0;
					if (eref != null) {
						offset += eref.getOffset(confSpace, pos, rc);
					}
					if (addResEntropy) {
						offset += getResEntropy(confSpace, pos, rc);
					}
					return offset;
				})
				.addInters(tuple)
				.addShell(tuple, weight, (pos, rc, shellResNum) -> 0.0)
				.make();
		}

		@Override
		public ResidueInteractions makeTripleCorrection(SimpleConfSpace confSpace, SimpleReferenceEnergies eref, boolean addResEntropy, int pos1, int rc1, int pos2, int rc2, int pos3, int rc3) {
			double offset1 = 0.0;
			double offset2 = 0.0;
			double offset3 = 0.0;
			if (eref != null) {
				offset1 += eref.getOffset(confSpace, pos1, rc1);
				offset2 += eref.getOffset(confSpace, pos2, rc2);
				offset3 += eref.getOffset(confSpace, pos3, rc3);
			}
			if (addResEntropy) {
				offset1 += getResEntropy(confSpace, pos1, rc1);
				offset2 += getResEntropy(confSpace, pos2, rc2);
				offset3 += getResEntropy(confSpace, pos3, rc3);
			}
			double singleWeight = 2.0/(confSpace.positions.size() - 2)/(confSpace.positions.size() - 1);
			double pairWeight = 1.0/(confSpace.positions.size() - 2);
			return ResInterGen.of(confSpace)
				.addIntra(pos1, singleWeight, offset1)
				.addIntra(pos2, singleWeight, offset2)
				.addIntra(pos3, singleWeight, offset3)
				.addInter(pos1, pos2, pairWeight, 0.0)
				.addInter(pos1, pos3, pairWeight, 0.0)
				.addInter(pos2, pos3, pairWeight, 0.0)
				.addShell(pos1, singleWeight, 0.0)
				.addShell(pos2, singleWeight, 0.0)
				.addShell(pos3, singleWeight, 0.0)
				.make();
		}

		@Override
		public double offsetTripleEnergy(int pos1, int rc1, int pos2, int rc2, int pos3, int rc3, EnergyMatrix emat) {
			// NOTE: no energy on singles, so don't need to add those here
			double weight = 1.0/(emat.getNumPos() - 2);
			return weight*(
				  emat.getPairwise(pos1, rc1, pos2, rc2)
				+ emat.getPairwise(pos1, rc1, pos3, rc3)
				+ emat.getPairwise(pos2, rc2, pos3, rc3)
			);
		}
	};

	public abstract ResidueInteractions makeSingle(SimpleConfSpace confSpace, SimpleReferenceEnergies eref, boolean addResEntropy, int pos, int rc);
	public abstract ResidueInteractions makePair(SimpleConfSpace confSpace, SimpleReferenceEnergies eref, boolean addResEntropy, int pos1, int rc1, int pos2, int rc2);
	public abstract ResidueInteractions makeTuple(SimpleConfSpace confSpace, SimpleReferenceEnergies eref, boolean addResEntropy, RCTuple frag);

	public abstract ResidueInteractions makeTripleCorrection(SimpleConfSpace confSpace, SimpleReferenceEnergies eref, boolean addResEntropy, int pos1, int rc1, int pos2, int rc2, int pos3, int rc3);
	public abstract double offsetTripleEnergy(int pos1, int rc1, int pos2, int rc2, int pos3, int rc3, EnergyMatrix emat);
	
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
