package edu.duke.cs.osprey.sofea;


import edu.duke.cs.osprey.confspace.*;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.tools.MathTools.BigDecimalBounds;
import edu.duke.cs.osprey.tools.MathTools.DoubleBounds;

import static edu.duke.cs.osprey.tools.Log.log;


public class SequenceLMFE implements Sofea.Criterion {

	public final Sequence seq;
	public final MultiStateConfSpace.LMFE lmfe;
	public final double maxLMFEWidth;

	public SequenceLMFE(Sequence seq, MultiStateConfSpace.LMFE lmfe, double maxLMFEWidth) {
		this.seq = seq;
		this.lmfe = lmfe;
		this.maxLMFEWidth = maxLMFEWidth;
	}

	@Override
	public Filter filterNode(MultiStateConfSpace.State state, int[] conf, BoltzmannCalculator bcalc) {

		for (SimpleConfSpace.Position confPos : state.confSpace.positions) {

			// skip non-mutable positions
			if (confPos.seqPos == null) {
				continue;
			}

			// skip unassigned positions
			int rc = conf[confPos.index];
			if (rc == Conf.Unassigned) {
				continue;
			}

			// if the RTs don't match, requeue this node
			SeqSpace.ResType seqResType = seq.get(confPos.seqPos.resNum);
			String confRT = confPos.resConfs.get(rc).template.name;
			if (!seqResType.name.equalsIgnoreCase(confRT)) {
				return Filter.Requeue;
			}
		}

		// the sequence matches, process this node
		return Filter.Process;
	}

	@Override
	public Satisfied isSatisfied(SeqDB seqdb, FringeDB fringedb, long sweepCount, BoltzmannCalculator bcalc) {

		SeqDB.SeqInfo seqInfo = seqdb.getSequencedZSumBounds(seq);

		// get the state free energies
		DoubleBounds[] stateFreeEnergies = lmfe.collectFreeEnergies(state -> {
			BigDecimalBounds zSumBounds;
			if (state.isSequenced) {
				zSumBounds = seqInfo.zSumBounds[state.sequencedIndex];
			} else {
				zSumBounds = seqdb.getUnsequencedZSumBounds(state);
			}
			return bcalc.freeEnergyPrecise(zSumBounds);
		});

		// compute the lmfe bounds
		DoubleBounds lmfeBounds = lmfe.calc().addAll(stateFreeEnergies).bounds;

		boolean isSatisfied = lmfeBounds.size() <= maxLMFEWidth;

		// report progress
		log("sequence [%s]:", seq);
		log("%10s=%s w=%9.4f  <  %9.4f   ? %s",
			"LMFE",
			lmfeBounds.toString(4, 9),
			lmfeBounds.size(),
			maxLMFEWidth,
			isSatisfied ? "yes, finished" : "no, keep refining"
		);
		for (MultiStateConfSpace.State state : seqdb.confSpace.states) {
			log("%10s=%s w=%9.4f",
				state.name,
				stateFreeEnergies[state.index].toString(4, 9),
				stateFreeEnergies[state.index].size()
			);
		}

		return isSatisfied ? Satisfied.Terminate : Satisfied.KeepIterating;
	}
}
