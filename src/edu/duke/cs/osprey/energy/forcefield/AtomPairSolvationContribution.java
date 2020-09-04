package edu.duke.cs.osprey.energy.forcefield;

public class AtomPairSolvationContribution extends AtomPairEnergyContribution {

	public AtomPairSolvationContribution(ResPairCache.ResPair resPair, int atom1, int atom2, double energy) {
		super(resPair, atom1, atom2, energy);
	}
}
