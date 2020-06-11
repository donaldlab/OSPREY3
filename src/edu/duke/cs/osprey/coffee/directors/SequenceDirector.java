package edu.duke.cs.osprey.coffee.directors;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.coffee.Coffee;
import edu.duke.cs.osprey.coffee.FreeEnergyCalculator;
import edu.duke.cs.osprey.coffee.NodeProcessor;
import edu.duke.cs.osprey.coffee.directions.Directions;
import edu.duke.cs.osprey.confspace.MultiStateConfSpace;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.tools.MathTools.BigDecimalBounds;
import edu.duke.cs.osprey.tools.MathTools.DoubleBounds;


public class SequenceDirector implements Coffee.Director {

	public final MultiStateConfSpace confSpace;
	public final Sequence seq;
	public final double gWidthMax;

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

			// tell the cluster to focus on this state
			directions.focus(state.index);

			// process the state until the free energy is precise enough
			while (getFreeEnergy(processor, state).size() > gWidthMax) {

				// TEMP: just process nodes locally for now
				var result = processor.process(directions);
				if (!result.processedNode) {
					directions.member.log("ran out of nodes early somehow!");
					break;
				}
			}
		}

		// all done, stop the computation
		directions.stop();

		// TEMP
		directions.member.log("seqdb: %s", processor.seqdb.dump());

		// show the results
		for (var state : confSpace.states) {
			var g = getFreeEnergy(processor, state);
			directions.member.log("%20s   G = %s  w=%.2f", state.name, g, g.size());
		}
	}

	private DoubleBounds getFreeEnergy(NodeProcessor processor, MultiStateConfSpace.State state) {

		// get the zSum bounds
		BigDecimalBounds zSum;
		if (state.isSequenced) {
			zSum = processor.seqdb.bounds(seq).zSumBounds[state.sequencedIndex];
		} else {
			zSum = processor.seqdb.boundsUnsequenced(state);
		}

		return new FreeEnergyCalculator().calc(zSum);
	}
}
