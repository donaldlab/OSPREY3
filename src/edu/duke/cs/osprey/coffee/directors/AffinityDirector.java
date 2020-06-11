package edu.duke.cs.osprey.coffee.directors;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.coffee.BestSequences;
import edu.duke.cs.osprey.coffee.Coffee;
import edu.duke.cs.osprey.coffee.NodeProcessor;
import edu.duke.cs.osprey.coffee.directions.Directions;
import edu.duke.cs.osprey.confspace.MultiStateConfSpace;

import java.util.concurrent.TimeUnit;

import static edu.duke.cs.osprey.tools.Log.log;


/**
 * Directs COFFEE to find the K best sequences by binding affinity
 *
 * The complex and design should be mutable, the target should not.
 */
public class AffinityDirector implements Coffee.Director {

	public final MultiStateConfSpace confSpace;
	public final MultiStateConfSpace.State complex;
	public final MultiStateConfSpace.State design;
	public final MultiStateConfSpace.State target;
	public final int K;

	/**
	 * The multi-state objective function to be minimized.
	 */
	public final MultiStateConfSpace.LMFE objective;


	public AffinityDirector(MultiStateConfSpace confSpace, String complexName, String designName, String targetName, int K) {
		this(
			confSpace,
			confSpace.getState(complexName),
			confSpace.getState(designName),
			confSpace.getState(targetName),
			K
		);
	}

	public AffinityDirector(MultiStateConfSpace confSpace, MultiStateConfSpace.State complex, MultiStateConfSpace.State design, MultiStateConfSpace.State target, int K) {

		// complex and design states should be sequenced, target state should not

		this.confSpace = confSpace;
		this.complex = complex;
		this.design = design;
		this.target = target;
		this.K = K;

		objective = confSpace.lmfe()
			.addPositive(complex)
			.addNegative(design)
			.addNegative(target)
			.build();
	}

	public void init(Directions directions, NodeProcessor processor) {

		// set the node trees to the whole space
		directions.setTrees(confSpace.states.stream()
			.map(state -> new RCs(state.confSpace))
			.toArray(RCs[]::new)
		);
	}

	public void direct(Directions directions, NodeProcessor processor) {

		// where we at?
		// TODO: how to choose the "best" K sequences? LMFE isn't necessarily the right answer here
		var bestSequences = BestSequences.getK(processor.seqdb, K, objective);
		log("Best Sequences: %d/%d", bestSequences.sequences.size(), K);
		for (int i=0; i<bestSequences.sequences.size(); i++) {
			var scoredSeq = bestSequences.sequences.get(i);
			log("\t%3d: %s, %s", i, scoredSeq.sequence, scoredSeq.lmfe.toString(4));
		}

		// TEMP: focus on the complex state
		directions.focus(complex.index);

		// TEMP: wait 5 seconds, then stop
		directions.member.sleep(5, TimeUnit.SECONDS);
		directions.stop();

		// TEMP: output results
	}

	private void findBestSequences() {
	}
}
