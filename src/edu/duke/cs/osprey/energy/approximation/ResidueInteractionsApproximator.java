package edu.duke.cs.osprey.energy.approximation;

import edu.duke.cs.osprey.energy.ResidueInteractions;


public class ResidueInteractionsApproximator {

	public final ResidueInteractions inters;
	public final ApproximatedObjectiveFunction.Approximator approximator;

	public ResidueInteractionsApproximator(ResidueInteractions inters, ApproximatedObjectiveFunction.Approximator approximator) {
		this.inters = inters;
		this.approximator = approximator;
	}
}
