package edu.duke.cs.osprey.energy.forcefield;

public class AtomPairVanDerWaalsContribution extends AtomPairEnergyContribution {

	private final double aij;
	private final double bij;
	private final double r2;

	public AtomPairVanDerWaalsContribution(ResPairCache.ResPair resPair, int atom1, int atom2, double energy, double aij, double bij, double r2) {
		super(resPair, atom1, atom2, energy);

		this.aij = aij;
		this.bij = bij;
		this.r2 = r2;
	}

	public double getR2() {
		return r2;
	}

	public double getBij() {
		return bij;
	}

	public double getAij() {
		return aij;
	}
}
