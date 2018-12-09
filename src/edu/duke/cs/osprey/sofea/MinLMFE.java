package edu.duke.cs.osprey.sofea;


import edu.duke.cs.osprey.confspace.MultiStateConfSpace;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.tools.MathTools;
import edu.duke.cs.osprey.tools.MathTools.DoubleBounds;

import java.math.BigInteger;
import java.math.MathContext;
import java.util.*;

import static edu.duke.cs.osprey.tools.Log.log;
import static edu.duke.cs.osprey.tools.Log.logf;


/**
 * Finishes SOFEA computation when we've found the lowest K sequences by LMFE
 */
public class MinLMFE implements Sofea.Criterion {

	public final MultiStateConfSpace.LMFE objective;
	public final int numSequences;
	public final MathContext mathContext;

	private final BoltzmannCalculator bcalc;

	public MinLMFE(MultiStateConfSpace.LMFE objective, int numSequences, MathContext mathContext) {

		// make sure we have enough sequences
		BigInteger maxNumSequences = objective.confSpace.seqSpace.getNumSequences();
		if (BigInteger.valueOf(numSequences).compareTo(maxNumSequences) > 0) {
			throw new IllegalArgumentException(String.format("conf space only has %d sequences, can't find %d",
				maxNumSequences, numSequences
			));
		}

		this.objective = objective;
		this.numSequences = numSequences;
		this.mathContext = mathContext;

		this.bcalc = new BoltzmannCalculator(mathContext);
	}

	public class TopSequences {

		public final TreeSet<Sofea.SeqResult> sequences = new TreeSet<>(Comparator.comparing(r -> r.lmfeFreeEnergy.lower));
		public double nextLowest = Double.POSITIVE_INFINITY;

		public boolean isTop(DoubleBounds bounds) {
			return bounds.lower < nextLowest;
		}

		public void add(Sofea.SeqResult result) {

			assert (isTop(result.lmfeFreeEnergy));

			sequences.add(result);

			// trim down to K sequences, and update the cutoff
			while (sequences.size() > numSequences) {
				nextLowest = sequences.pollLast().lmfeFreeEnergy.lower;
			}
		}
	}

	public TopSequences getTopSequences(SeqDB seqdb) {

		assert (objective.confSpace == seqdb.confSpace);
		MultiStateConfSpace confSpace = seqdb.confSpace;

		TopSequences topSequences = new TopSequences();

		// get the unsequenced z values
		MathTools.DoubleBounds[] unsequencedFreeEnergy = new MathTools.DoubleBounds[confSpace.unsequencedStates.size()];
		for (MultiStateConfSpace.State state : confSpace.unsequencedStates) {
			unsequencedFreeEnergy[state.unsequencedIndex] = bcalc.freeEnergyPrecise(seqdb.getUnsequenced(state.unsequencedIndex));
		}

		// for each sequence and partial sequence encountered so far...
		for (Map.Entry<Sequence,SeqDB.SeqInfo> entry : seqdb.getSequences()) {

			Sequence seq = entry.getKey();
			SeqDB.SeqInfo seqInfo = entry.getValue();

			DoubleBounds[] stateFreeEnergies = new DoubleBounds[confSpace.states.size()];
			Arrays.fill(stateFreeEnergies, null);

			// compute bounds on the objective function
			MathTools.DoubleBounds objectiveBounds = new MathTools.DoubleBounds(0.0, 0.0);
			for (MultiStateConfSpace.State state : objective.states()) {

				MathTools.DoubleBounds g;
				if (state.isSequenced) {
					g = bcalc.freeEnergyPrecise(seqInfo.bounds[state.sequencedIndex]);
				} else {
					g = unsequencedFreeEnergy[state.unsequencedIndex];
				}

				objectiveBounds.lower += g.lower* objective.getWeight(state);
				objectiveBounds.upper += g.upper* objective.getWeight(state);

				stateFreeEnergies[state.index] = g;
			}

			// update the top sequences
			if (topSequences.isTop(objectiveBounds)) {

				// yup, get the rest of the free energies if needed
				for (MultiStateConfSpace.State state : confSpace.states) {
					if (stateFreeEnergies[state.index] == null) {
						if (state.isSequenced) {
							stateFreeEnergies[state.index] = bcalc.freeEnergyPrecise(seqInfo.bounds[state.sequencedIndex]);
						} else {
							stateFreeEnergies[state.index] = unsequencedFreeEnergy[state.unsequencedIndex];
						}
					}
				}

				topSequences.add(new Sofea.SeqResult(
					seq,
					objectiveBounds,
					stateFreeEnergies
				));
			}
		}

		return topSequences;
	}

	@Override
	public boolean isFinished(SeqDB seqdb) {

		assert (objective.confSpace == seqdb.confSpace);
		MultiStateConfSpace confSpace = seqdb.confSpace;

		TopSequences topSequences = getTopSequences(seqdb);

		// TEMP
		log("lowest %d/%d sequences by LMFE:", topSequences.sequences.size(), numSequences);
		int i = 0;
		for (Sofea.SeqResult result : topSequences.sequences) {
			logf("\tseq %2d  obj=%s w=%9.4f",
				i++,
				result.lmfeFreeEnergy.toString(4, 9),
				result.lmfeFreeEnergy.size()
			);
			for (MultiStateConfSpace.State state : confSpace.sequencedStates) {
				DoubleBounds stateFreeEnergy = result.stateFreeEnergies[state.index];
				logf("    %s=%s w=%.4f",
					state.name,
					stateFreeEnergy.toString(4, 9),
					stateFreeEnergy.size()
				);
			}
			log("    [%s]", result.sequence);
		}
		log("\tnext lowest: %9.4f", topSequences.nextLowest);

		// do we have enough sequences yet?
		if (topSequences.sequences.size() < numSequences) {
			return false;
		}

		// all the top K sequences must be ...
		for (Sofea.SeqResult result : topSequences.sequences) {

			// fully assigned
			if (!result.sequence.isFullyAssigned()) {
				return false;
			}

			// have finite lower bounds
			if (!Double.isFinite(result.lmfeFreeEnergy.lower)) {
				return false;
			}

			// must not overlap the next lowest
			if (result.lmfeFreeEnergy.upper >= topSequences.nextLowest) {
				return false;
			}
		}

		// all is well!
		return true;
	}
}
