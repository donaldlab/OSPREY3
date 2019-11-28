package edu.duke.cs.osprey.confspace.compiled;

import java.util.ArrayList;
import java.util.List;


/**
 * Defines how position interactions should be distributed
 * among conformation fragments.
 */
public enum PosInterDist {

	/**
	 * Uses the traditional distribution introduced in the original DEE paper:
	 *
	 * The dead-end elimination theorem and its use in protein side-chain positioning
	 * J Desmet, M De Maeyer, B Hazes, I Lasters
	 * Nature, 1992
	 * https://doi.org/10.1038/356539a0
	 *
	 * Namely, pos-static and pos interactions are placed on "single" conf fragments.
	 * And pos-pos interactions are placed on "pair" conf fragments.
	 */
	DesmetEtAl1992 {

		/**
		 * pos and pos-static interactions go on the single conf.
		 */
		@Override
		public List<PosInter> single(ConfSpace confSpace, int posi1) {
			List<PosInter> inters = new ArrayList<>();
			inters.add(new PosInter(posi1, posi1, 1.0));
			inters.add(new PosInter(posi1, PosInter.StaticPos, 1.0));
			return inters;
		}

		/**
		 * pos-pos interactions go on the pair conf.
		 */
		@Override
		public List<PosInter> pair(ConfSpace confSpace, int posi1, int posi2) {
			List<PosInter> inters = new ArrayList<>();
			inters.add(new PosInter(posi1, posi2, 1.0));
			return inters;
		}
	},

	/**
	 * A newer distribution that yields tighter lower bounds on conformation energy.
	 * No interactions are placed on "single" conf fragments.
	 * Instead, the pos-static and pos interactions are distributed evenly among the
	 * "pair" conf fragments involving that position, without double-counting.
	 * This re-distribution allows minimizing "pair" conf fragments in the presence of the
	 * fixed atoms, which generally results in less optimistic lower bounds for "pair" energies.
	 */
	TighterBounds {

		/**
		 * No interactions go on the single conf.
		 */
		@Override
		public List<PosInter> single(ConfSpace confSpace, int posi1) {
			return new ArrayList<>();
		}

		/**
		 * All interactions involving pos1 and pos2 go on the pair,
		 * but the pos and pos-static interactions are weighted to avoid double-counting.
		 */
		@Override
		public List<PosInter> pair(ConfSpace confSpace, int posi1, int posi2) {

			double weight = 1.0/(confSpace.positions.length - 1);

			List<PosInter> inters = new ArrayList<>();
			inters.add(new PosInter(posi1, posi1, weight));
			inters.add(new PosInter(posi2, posi2, weight));
			inters.add(new PosInter(posi1, PosInter.StaticPos, weight));
			inters.add(new PosInter(posi2, PosInter.StaticPos, weight));
			inters.add(new PosInter(posi1, posi2, 1.0));
			return inters;
		}
	};

	/**
	 * Include all interactions. Doesn't depend on the distribution.
	 */
	public static List<PosInter> all(ConfSpace confSpace) {
		List<PosInter> inters = new ArrayList<>();
		for (int posi1=0; posi1<confSpace.positions.length; posi1++) {

			// pos and pos-static interactions go on singles
			inters.add(new PosInter(posi1, posi1, 1.0));
			inters.add(new PosInter(posi1, PosInter.StaticPos, 1.0));

			// pos-pos interactions go on pairs
			for (int posi2=0; posi2<posi1; posi2++) {
				inters.add(new PosInter(posi1, posi2, 1.0));
			}
		}
		return inters;
	}

	public abstract List<PosInter> single(ConfSpace confSpace, int posi1);
	public abstract List<PosInter> pair(ConfSpace confSpace, int posi1, int posi2);
}
