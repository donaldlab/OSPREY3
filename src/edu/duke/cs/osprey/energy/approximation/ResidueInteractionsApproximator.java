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

package edu.duke.cs.osprey.energy.approximation;

import edu.duke.cs.osprey.dof.DofInfo;
import edu.duke.cs.osprey.energy.approximation.ApproximatedObjectiveFunction.Approximator;
import edu.duke.cs.osprey.energy.ResidueInteractions;


public class ResidueInteractionsApproximator {

	public static class Builder {

		public final DofInfo dofInfo;

		private Approximator.Addable approximator = null;
		private ResidueInteractions approxInters = new ResidueInteractions();
		private ResidueInteractions ffInters = new ResidueInteractions();

		public Builder(DofInfo dofInfo) {
			this.dofInfo = dofInfo;
		}

		public boolean isEmpty() {
			return approximator == null;
		}

		public double error() {
			if (approximator == null) {
				return 0;
			} else {
				return approximator.error();
			}
		}

		public void approximate(ResidueInteractions.Pair inter, Approximator.Addable approximator) {
			init(approximator);
			this.approximator.add(approximator, inter.weight, inter.offset);
			approxInters.add(inter);
		}

		public void dontApproximate(ResidueInteractions.Pair inter) {
			ffInters.add(inter);
		}

		private void init(Approximator.Addable approximator) {
			if (this.approximator == null) {
				this.approximator = approximator.makeIdentity(dofInfo.ids, dofInfo.counts);
			}
		}

		public ResidueInteractionsApproximator build() {

			if (approximator != null) {
				return new ResidueInteractionsApproximator(
					approximator,
					ffInters,
					approxInters
				);
			}

			// otherwise, make a no-op approximator
			return new ResidueInteractionsApproximator(
				new NOPApproximator(dofInfo.numDofs()),
				ffInters,
				new ResidueInteractions()
			);
		}
	}

	public final Approximator approximator;
	public final ResidueInteractions ffInters;
	public final ResidueInteractions approxInters;

	public ResidueInteractionsApproximator(Approximator approximator, ResidueInteractions ffInters, ResidueInteractions approxInters) {
		this.approximator = approximator;
		this.ffInters = ffInters;
		this.approxInters = approxInters;
	}
}
