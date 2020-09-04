package edu.duke.cs.osprey.energy.forcefield;

import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

public class ResPairEnergyContribution {

	public static final ResPairEnergyContribution BrokenConformation = new ResPairEnergyContribution(null, null);
	private final ResPairCache.ResPair resPair;
	private final List<AtomPairEnergyContribution> atomPairs;

	public ResPairEnergyContribution(ResPairCache.ResPair resPair, List<AtomPairEnergyContribution> atomPairs) {
		this.atomPairs = atomPairs;
		this.resPair = resPair;
	}

	public double getEnergy() {
		var energy = atomPairs.stream().map(AtomPairEnergyContribution::getEnergy).reduce(Double::sum).orElseThrow();
		return (energy + resPair.offset + resPair.solvEnergy) * resPair.weight;
	}

	public ResPairCache.ResPair getResPair() {
		return resPair;
	}

	public List<AtomPairEnergyContribution> getTopEnergyContributions(int n) {
		return atomPairs.stream()
				.sorted(Comparator.comparingDouble(AtomPairEnergyContribution::getEnergy).reversed())
				.limit(n)
				.collect(Collectors.toList());
	}

	public List<AtomPairEnergyContribution> getAtomPairEnergyContributions() {
		return atomPairs;
	}
}
