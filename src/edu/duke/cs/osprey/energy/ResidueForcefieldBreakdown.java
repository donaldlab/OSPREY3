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
import edu.duke.cs.osprey.energy.forcefield.ResidueForcefieldEnergy;

import java.util.*;
import java.util.stream.Collectors;

public class ResidueForcefieldBreakdown {

	public static enum Type {

		Electrostatics(true) {
			@Override
			public double getEnergy(ResidueForcefieldEnergy efunc) {
				return efunc.getElectrostaticsEnergy();
			}
		},
		VanDerWaals(true) {
			@Override
			public double getEnergy(ResidueForcefieldEnergy efunc) {
				return efunc.getVanDerWaalsEnergy();
			}
		},
		Solvation(true) {
			@Override
			public double getEnergy(ResidueForcefieldEnergy efunc) {
				return efunc.getSolvationEnergy();
			}
		},
		Offsets(true) {
			@Override
			public double getEnergy(ResidueForcefieldEnergy efunc) {
				return efunc.getOffsetsEnergy();
			}
		},
		All(false) {
			@Override
			public double getEnergy(ResidueForcefieldEnergy efunc) {
				return efunc.getElectrostaticsEnergy()
					+ efunc.getVanDerWaalsEnergy()
					+ efunc.getSolvationEnergy()
					+ efunc.getOffsetsEnergy();
			}
		};

		public final boolean isAtomic;

		private Type(boolean isAtomic) {
			this.isAtomic = isAtomic;
		}

		public abstract double getEnergy(ResidueForcefieldEnergy efunc);

		public static List<Type> atomics() {
			return Arrays.stream(values())
				.filter((val) -> val.isAtomic)
				.collect(Collectors.toList());
		}
	}

	public static class ByResidue {

		public final ResidueForcefieldEnergy efunc;

		public ByResidue(ResidueForcefieldEnergy efunc) {
			this.efunc = efunc;
		}

		public ByResidue(ConfEnergyCalculator confEcalc, int[] assignments) {
			this(confEcalc, confEcalc.calcEnergy(new RCTuple(assignments)));
		}

		public ByResidue(ConfEnergyCalculator confEcalc, EnergyCalculator.EnergiedParametricMolecule epmol) {
			this.efunc = (ResidueForcefieldEnergy)confEcalc.ecalc.makeEnergyFunction(epmol);
		}

		public EnergyMatrix makeEmat() {
			int[] numRCsAtPos = new int[efunc.residues.size()];
			Arrays.fill(numRCsAtPos, 1);
			return new EnergyMatrix(efunc.residues.size(), numRCsAtPos, 0.0);
		}

		public EnergyMatrix breakdownForcefield(Type type) {
			EnergyMatrix breakdown = makeEmat();
			for (ResidueInteractions.Pair pair : efunc.inters) {
				int i1 = efunc.residues.findIndexOrThrow(pair.resNum1);
				int i2 = efunc.residues.findIndexOrThrow(pair.resNum2);
				double energy = type.getEnergy(efunc.makeSubset(pair));
				if (i1 == i2) {
					breakdown.setOneBody(i1, 0, energy);
				} else {
					breakdown.setPairwise(i1, 0, i2, 0, energy);
				}
			}
			return breakdown;
		}
	}

	public static class ByPosition {

		public final ConfEnergyCalculator confEcalc;
		public final int[] assignments;
		public final EnergyCalculator.EnergiedParametricMolecule epmol;
		public final ResidueForcefieldEnergy efunc;

		public ByPosition(ConfEnergyCalculator confEcalc, int[] assignments) {
			this(confEcalc, assignments, confEcalc.calcEnergy(new RCTuple(assignments)));
		}

		public ByPosition(ConfEnergyCalculator confEcalc, int[] assignments, EnergyCalculator.EnergiedParametricMolecule epmol) {

			this.confEcalc = confEcalc;
			this.assignments = assignments;
			this.epmol = epmol;

			// get the forcefield
			this.efunc = (ResidueForcefieldEnergy)confEcalc.ecalc.makeEnergyFunction(epmol);
		}

		public EnergyMatrix makeEmat() {
			int[] numRCsAtPos = new int[confEcalc.confSpace.positions.size()];
			Arrays.fill(numRCsAtPos, 1);
			return new EnergyMatrix(confEcalc.confSpace.positions.size(), numRCsAtPos, 0.0);
		}

		public EnergyMatrix breakdownForcefield(Type type) {
			EnergyMatrix breakdown = makeEmat();
			for (SimpleConfSpace.Position pos1 : confEcalc.confSpace.positions) {
				int rc1 = assignments[pos1.index];

				{
					ResidueInteractions inters = confEcalc.makeSingleInters(pos1.index, rc1);
					double energy = type.getEnergy(efunc.makeSubset(inters));
					breakdown.setOneBody(pos1.index, 0, energy);
				}

				for (SimpleConfSpace.Position pos2 : confEcalc.confSpace.positions) {
					if (pos2.index < pos1.index) {
						int rc2 = assignments[pos2.index];

						ResidueInteractions inters = confEcalc.makePairInters(pos1.index, rc1, pos2.index, rc2);
						double energy = type.getEnergy(efunc.makeSubset(inters));
						breakdown.setPairwise(pos1.index, 0, pos2.index, 0, energy);
					}
				}
			}
			return breakdown;
		}

		public EnergyMatrix breakdownScore(EnergyMatrix emat) {
			EnergyMatrix breakdown = makeEmat();
			for (SimpleConfSpace.Position pos1 : confEcalc.confSpace.positions) {
				int rc1 = assignments[pos1.index];

				{
					double energy = emat.getOneBody(pos1.index, rc1);
					breakdown.setOneBody(pos1.index, 0, energy);
				}

				for (SimpleConfSpace.Position pos2 : confEcalc.confSpace.positions) {
					if (pos2.index < pos1.index) {
						int rc2 = assignments[pos2.index];

						double energy = emat.getPairwise(pos1.index, rc1, pos2.index, rc2);
						breakdown.setPairwise(pos1.index, 0, pos2.index, 0, energy);
					}
				}
			}
			return breakdown;
		}
	}
}
