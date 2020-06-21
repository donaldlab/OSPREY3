package edu.duke.cs.osprey.coffee.directors;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.coffee.Coffee;
import edu.duke.cs.osprey.coffee.FreeEnergyCalculator;
import edu.duke.cs.osprey.coffee.NodeProcessor;
import edu.duke.cs.osprey.coffee.directions.Directions;
import edu.duke.cs.osprey.coffee.seqdb.StateZ;
import edu.duke.cs.osprey.confspace.MultiStateConfSpace;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.tools.MathTools.DoubleBounds;
import edu.duke.cs.osprey.tools.Stopwatch;

import java.util.concurrent.TimeUnit;


public class SequenceDirector implements Coffee.Director {

	public final MultiStateConfSpace confSpace;
	public final Sequence seq;
	public final double gWidthMax;

	private final FreeEnergyCalculator gcalc = new FreeEnergyCalculator();

	private DoubleBounds[] freeEnergies;

	public SequenceDirector(MultiStateConfSpace confSpace, Sequence seq, double gWidthMax) {
		this.confSpace = confSpace;
		this.seq = seq;
		this.gWidthMax = gWidthMax;
	}

	@Override
	public void init(Directions directions, NodeProcessor processor) {

		// set the node trees for each state to just the specified sequence
		directions.setTrees(confSpace.states.stream()
			.map(state -> seq.makeRCs(state.confSpace))
			.toArray(RCs[]::new)
		);
	}

	@Override
	public void direct(Directions directions, NodeProcessor processor) {

		// process all the states to the desired precision
		for (var state : confSpace.states) {

			directions.member.log("Processing state %s", state.name);
			Stopwatch stopwatch = new Stopwatch().start();

			// tell the cluster to focus on this state
			directions.focus(state.index);

			// process the state until the free energy is precise enough
			while (true) {

				// process nodes for a little bit
				// within the first 5 minutes of computation, use smaller batches, but gradually grow to 1 min batches
				long nodesTimeMs = Math.max(500L, (long)(Math.min(stopwatch.getTimeM()/5.0, 1.0)*60_000));
				var foundNodes = processor.processFor(directions, nodesTimeMs, TimeUnit.MILLISECONDS);
				if (!foundNodes) {
					break;
				}

				// TEMP
				synchronized (processor.seqdb) {

					// check progress
					StateZ statez = getStateZ(processor, state);
					DoubleBounds g = gcalc.calc(statez.zSumBounds);
					double gWidth = g.size();

					// what's the best precision we could ever get given the nodes we've already dropped?
					double gWidthMin = Math.abs(gcalc.calc(
						processor.seqdb.bigMath()
							.set(statez.zSumBounds.upper)
							.sub(statez.zSumDropped)
							.recip()
							.mult(statez.zSumBounds.upper)
							.get()
					));

					/* PROOF:
						At step i, zSum_i is bounded by [min_i, max_i].
						Let U_i = max_i - min_i, be the uncertainty at step i.
						Let D_i be the part of the uncertainty due to dropped nodes at step i.

						When the computations finishes at step n, the same still holds:
							zSum_n is bounded by [min_n, max_n].
						Assuming we didn't drop any more nodes since step i, then:
							U_n = D_i
							min_n >= min_i
							max_n <= max_i.

						By definition:
							min_n = max_n - U_n.
						Then after some substitutions:
							min_n = max_n - D_i
								 <= max_i - D_i.

						So min_n is bounded by [min_i, max_i - D_i].

						At step n, zSum_n is bounded by [min_n, min_n + U_n] by definition
							and therefore by [min_n, min_n + D_i], snce U_n = D_i.

						After converting bounds on zSum_n to free energy, the width of the free energy bounds depends
						not only on D_i, but also on the value of min_n.
						Since we're intersted in the minimim free energy width possible after dropping nodes,
						we must choose the min_n in [min_i, max_i - D_i] that minimizes the free energy width.
						The optimial value for min_n happens to be max_i - D_i, so that's convenient. =P
						(proof of the last step omitted, it's not the hard step in this proof)
					*/

					// what fraction of the uncertainty is not dropped, and hence still reducible?
					// (this should asymptotically approach zero)
					var reducibleRatio = processor.seqdb.bigMath()
						.set(statez.zSumBounds.upper)
						.sub(statez.zSumBounds.lower)
						.div(statez.zSumDropped)
						.get()
						.doubleValue()
						- 1.0;

					// report progress
					directions.member.log("\tfree energy %s   width %.6f of %.6f  nodedb usage %5.1f%%   r %.6f   time %s",
						g.toString(6), gWidth, gWidthMin,
						processor.nodedb.usage()*100f,
						reducibleRatio,
						stopwatch.getTime(2)
					);
					// TODO: show node processing speeds?

					// are we there yet?
					if (g.size() <= gWidthMax) {
						break;
					}
				}
			}
		}

		// all done, stop the computation
		directions.stop();

		// collect the results
		freeEnergies = confSpace.states.stream()
			.map(state -> gcalc.calc(getStateZ(processor, state).zSumBounds))
			.toArray(DoubleBounds[]::new);
	}

	private StateZ getStateZ(NodeProcessor processor, MultiStateConfSpace.State state) {
		if (state.isSequenced) {
			return processor.seqdb.get(seq).statezs[state.sequencedIndex];
		} else {
			return processor.seqdb.getUnsequenced(state);
		}
	}

	public DoubleBounds getFreeEnergy(MultiStateConfSpace.State state) {
		return freeEnergies[state.index];
	}
}
