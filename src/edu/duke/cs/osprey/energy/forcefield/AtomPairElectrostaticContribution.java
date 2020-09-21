package edu.duke.cs.osprey.energy.forcefield;

public class AtomPairElectrostaticContribution extends AtomPairEnergyContribution {

	private final double radius;
	private final double charge;
	private final boolean is14Bonded;
	private final boolean isHeavyPair;
	private final boolean usesDistanceDependentDialectic;

	public AtomPairElectrostaticContribution(ResPairCache.ResPair resPair, int atom1, int atom2, double energy,
											 double radius, double charge, boolean is14Bonded,
											 boolean isHeavyPair, boolean usesDistanceDependentDialectic) {
		super(resPair, atom1, atom2, energy);
		this.radius = radius;
		this.charge = charge;
		this.is14Bonded = is14Bonded;
		this.isHeavyPair = isHeavyPair;
		this.usesDistanceDependentDialectic = usesDistanceDependentDialectic;
	}

	public double getRadius() {
		return radius;
	}

	public double getCharge() {
		return charge;
	}

	public boolean isIs14Bonded() {
		return is14Bonded;
	}

	public boolean isHeavyPair() {
		return isHeavyPair;
	}

	public boolean isUsesDistanceDependentDialectic() {
		return usesDistanceDependentDialectic;
	}
}
