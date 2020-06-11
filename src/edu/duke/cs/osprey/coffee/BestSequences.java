package edu.duke.cs.osprey.coffee;

import edu.duke.cs.osprey.coffee.seqdb.SeqDB;
import edu.duke.cs.osprey.coffee.seqdb.SeqInfo;
import edu.duke.cs.osprey.confspace.MultiStateConfSpace;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.tools.MathTools;
import edu.duke.cs.osprey.tools.UnpossibleError;

import java.math.MathContext;
import java.math.RoundingMode;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;


public class BestSequences {

	// TODO: keep this in sync with SeqDB instead of reading the whole SeqDB every time?
	// TODO: is this really useful for affinity designs?
	//  LMFE isn't quite what we want...
	//  what we really want is a bi-criterion optimization of the complex and design states
	//  without necessarily choosing a weighting between the two criteria beforehand

	/**
	 * Get the best K sequences by LMFE upper bound from the sequence database
	 */
	public static BestSequences getK(SeqDB seqdb, int K, MultiStateConfSpace.LMFE objective) {

		// only computing ln, so only need enough precision to fill a double
		var mc = new MathContext(16, RoundingMode.HALF_UP);
		var bcalc = new BoltzmannCalculator(mc);
		var bestSeqs = new BestSequences(K);

		// get the unsequenced z values
		MathTools.DoubleBounds[] unsequencedFreeEnergy = new MathTools.DoubleBounds[seqdb.confSpace.unsequencedStates.size()];
		for (MultiStateConfSpace.State state : seqdb.confSpace.unsequencedStates) {
			MathTools.DoubleBounds freeEnergyBounds = bcalc.freeEnergyPrecise(seqdb.boundsUnsequenced(state));
			unsequencedFreeEnergy[state.unsequencedIndex] = freeEnergyBounds;
		}

		// for each sequence and partial sequence encountered so far...
		for (var entry : seqdb.boundsSequenced()) {
			Sequence seq = entry.getKey();
			SeqInfo seqInfo = entry.getValue();

			MathTools.DoubleBounds[] stateFreeEnergies = objective.collectFreeEnergies(state -> {
				if (state.isSequenced) {
					return bcalc.freeEnergyPrecise(seqInfo.zSumBounds[state.sequencedIndex]);
				} else {
					return unsequencedFreeEnergy[state.unsequencedIndex];
				}
			});

			// compute bounds on the objective function
			MathTools.DoubleBounds objectiveBounds = objective.calc().addAll(stateFreeEnergies).bounds;

			// update the best sequences
			if (bestSeqs.isBest(objectiveBounds)) {

				// yup, get the free energies not in the objective, if needed
				for (MultiStateConfSpace.State state : seqdb.confSpace.states) {
					if (stateFreeEnergies[state.index] == null) {
						if (state.isSequenced) {
							stateFreeEnergies[state.index] = bcalc.freeEnergyPrecise(seqInfo.zSumBounds[state.sequencedIndex]);
						} else {
							stateFreeEnergies[state.index] = unsequencedFreeEnergy[state.unsequencedIndex];
						}
					}
				}

				bestSeqs.add(new ScoredSeq(
					seq,
					objectiveBounds,
					stateFreeEnergies
				));
			}
		}

		return bestSeqs;
	}


	public final int maxNumSequences;
	public final List<ScoredSeq> sequences = new ArrayList<>();

	public ScoredSeq nextLowest = null;

	public BestSequences(int maxNumSequences) {
		this.maxNumSequences = maxNumSequences;
	}

	public boolean isBest(MathTools.DoubleBounds bounds) {
		return nextLowest == null || bounds.upper < nextLowest.lmfe.upper;
	}

	public void add(ScoredSeq scoredSeq) {

		assert (isBest(scoredSeq.lmfe));

		sequences.add(scoredSeq);

		// sort the sequences by LMFE upper bound
		sequences.sort(Comparator.comparing(s -> s.lmfe.upper));

		// trim down to the best K sequences (allowing ties, but kick out +inf), and update the cutoff
		int numTopK = 0;
		Double lastVal = null;
		for (int i=0; i<sequences.size(); i++) {

			ScoredSeq r = sequences.get(i);
			double val = r.lmfe.upper;

			// immediately trim +inf out of best K
			if (val == Double.POSITIVE_INFINITY) {
				trimAt(i);
				break;
			}

			if (lastVal == null) {

				// no best K yet results, accept the first result
				numTopK++;
				lastVal = val;

			} else if (val == lastVal) {

				// same as last result in best K, allow ties
				numTopK++;

			} else if (val > lastVal) {

				// worse than last result in best K
				if (numTopK < maxNumSequences) {

					// have room for more in best K, add this result
					numTopK++;
					lastVal = val;

				} else {

					// no more room in best k, trim the rest of the list
					trimAt(i);
					break;
				}
			} else {
				throw new UnpossibleError();
			}
		}
	}

	private void trimAt(int i) {
		nextLowest = sequences.get(i);
		// but don't keep +inf sequences even at the next lowest
		if (nextLowest.lmfe.upper == Double.POSITIVE_INFINITY) {
			nextLowest = null;
		}
		sequences.subList(i, sequences.size()).clear();
	}
}
