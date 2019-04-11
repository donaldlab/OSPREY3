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

package edu.duke.cs.osprey.sofea;


import edu.duke.cs.osprey.confspace.*;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.tools.MathTools.BigDecimalBounds;
import edu.duke.cs.osprey.tools.MathTools.DoubleBounds;

import java.util.HashSet;
import java.util.Set;

import static edu.duke.cs.osprey.tools.Log.log;


/**
 * Finished SOFEA computation when the state free energies for a single sequence
 * reach a desired precision.
 */
public class SequenceLMFE implements Sofea.Criterion {

	/**
	 * The sequence.
	 */
	public final Sequence seq;

	/**
	 * The multi-state objective function.
	 */
	public final MultiStateConfSpace.LMFE lmfe;

	/**
	 * The best desired precision for a state free energy.
	 */
	public final double minFreeEnergyWidth;

	private final Set<Integer> finishedStates = new HashSet<>();

	public SequenceLMFE(Sequence seq, MultiStateConfSpace.LMFE lmfe, double minFreeEnergyWidth) {
		this.seq = seq;
		this.lmfe = lmfe;
		this.minFreeEnergyWidth = minFreeEnergyWidth;
	}

	@Override
	public Filter filterNode(MultiStateConfSpace.State state, int[] conf, BoltzmannCalculator bcalc) {

		// is this a state we care about?
		if (!lmfe.states().contains(state) || finishedStates.contains(state.index)) {
			return Filter.Requeue;
		}

		// is this a sequence we care about?
		Sequence nodeSeq = state.confSpace.seqSpace.makeSequence(state.confSpace, conf);
		if (!nodeSeq.equals(seq)) {
			return Filter.Requeue;
		}

		// the state and sequence matches, process this node
		return Filter.Process;
	}

	@Override
	public Satisfied isSatisfied(SeqDB seqdb, FringeDB fringedbLower, FringeDB fringedbUpper, long pass1Step, long pass2Step, BoltzmannCalculator bcalc) {

		// get the state free energies
		SeqDB.SeqInfo seqInfo = seqdb.getSequencedZSumBounds(seq);
		DoubleBounds[] stateFreeEnergies = lmfe.collectFreeEnergies(state -> {
			BigDecimalBounds zSumBounds;
			if (state.isSequenced) {
				zSumBounds = seqInfo.zSumBounds[state.sequencedIndex];
			} else {
				zSumBounds = seqdb.getUnsequencedZSumBounds(state);
			}
			return bcalc.freeEnergyPrecise(zSumBounds);
		});

		// are the state free energies precise enough?
		boolean isFinished = true;
		for (MultiStateConfSpace.State state : lmfe.states()) {
			if (stateFreeEnergies[state.index].size() <= minFreeEnergyWidth) {
				finishedStates.add(state.index);
			} else {
				isFinished = false;
			}
		}

		// report progress
		DoubleBounds lmfeBounds = lmfe.calc().addAll(stateFreeEnergies).bounds;
		log("sequence [%s]:", seq);
		log("%10s=%s w=%9.4f",
			"LMFE",
			lmfeBounds.toString(4, 9),
			lmfeBounds.size()
		);
		for (MultiStateConfSpace.State state : lmfe.states()) {
			DoubleBounds g = stateFreeEnergies[state.index];
			log("%10s=%s w=%9.4f  <  %9.4f   ? %s",
				state.name,
				g.toString(4, 9),
				g.size(),
				minFreeEnergyWidth,
				g.size() <= minFreeEnergyWidth ? "yes, finished" : "no, keep refining"
			);
		}

		return isFinished ? Satisfied.Terminate : Satisfied.KeepSweeping;
	}
}
