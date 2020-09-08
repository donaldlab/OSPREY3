package edu.duke.cs.osprey.energy.forcefield;

import edu.duke.cs.osprey.structure.Atom;

public class AtomPairEnergyContribution {

	public static final AtomPairEnergyContribution BrokenConformation = new AtomPairEnergyContribution(null, 0, 0, Double.POSITIVE_INFINITY);

	private final int atom1;
	private final int atom2;
	private final double energy;
	private final ResPairCache.ResPair resPair;

	public AtomPairEnergyContribution(ResPairCache.ResPair resPair, int atom1, int atom2, double energy) {
		this.resPair = resPair;
		this.atom1 = atom1;
		this.atom2 = atom2;
		this.energy = energy;
	}

	public double getEnergy() {
		return energy;
	}

	public Atom getAtom2() {
		return resPair.res2.atoms.get(atom2);
	}

	public Atom getAtom1() {
		return resPair.res1.atoms.get(atom1);
	}
}
