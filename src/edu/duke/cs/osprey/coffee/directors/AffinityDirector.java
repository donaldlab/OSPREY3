package edu.duke.cs.osprey.coffee.directors;

import edu.duke.cs.osprey.coffee.Coffee;
import edu.duke.cs.osprey.coffee.NodeProcessor;
import edu.duke.cs.osprey.coffee.commands.Commands;
import edu.duke.cs.osprey.confspace.MultiStateConfSpace;

import java.util.concurrent.TimeUnit;


/**
 * Directs COFFEE to find the K best sequences with good binding affinity
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

	public void direct(Commands commands, NodeProcessor processor) {

		// TEMP: focus on the complex state
		commands.focus(complex.index);

		// TEMP: wait 5 seconds, then stop
		commands.member.sleep(5, TimeUnit.SECONDS);
		commands.stop();

		// TEMP: output results
	}
}
