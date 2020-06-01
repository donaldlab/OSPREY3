package edu.duke.cs.osprey.coffee.drivers;

import edu.duke.cs.osprey.coffee.Coffee;
import edu.duke.cs.osprey.confspace.MultiStateConfSpace;


/**
 * Drives COFFEE to find the K best sequences with good binding affinity
 */
public class AffinityDriver implements Coffee.Driver {

	public final MultiStateConfSpace confSpace;
	public final MultiStateConfSpace.State complex;
	public final MultiStateConfSpace.State design;
	public final MultiStateConfSpace.State target;
	public final int K;

	/**
	 * The multi-state objective function to be minimized.
	 */
	public final MultiStateConfSpace.LMFE objective;


	public AffinityDriver(MultiStateConfSpace confSpace, String complexName, String designName, String targetName, int K) {
		this(
			confSpace,
			confSpace.getState(complexName),
			confSpace.getState(designName),
			confSpace.getState(targetName),
			K
		);
	}

	public AffinityDriver(MultiStateConfSpace confSpace, MultiStateConfSpace.State complex, MultiStateConfSpace.State design, MultiStateConfSpace.State target, int K) {

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
}
