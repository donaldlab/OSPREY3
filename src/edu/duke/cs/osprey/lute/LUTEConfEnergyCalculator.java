package edu.duke.cs.osprey.lute;

import edu.duke.cs.osprey.confspace.Conf;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.ResidueInteractions;


public class LUTEConfEnergyCalculator extends ConfEnergyCalculator {

	public final LUTE lute;

	public LUTEConfEnergyCalculator(LUTE lute, EnergyCalculator ecalc) {
		super(lute.confSpace, ecalc, null, null, false);

		this.lute = lute;
	}

	@Override
	public ResidueInteractions makeSingleInters(int pos, int rc) {
		throw new UnsupportedOperationException("LUTE can only be used to compute full-conformation energies");
	}

	@Override
	public ResidueInteractions makePairInters(int pos1, int rc1, int pos2, int rc2) {
		throw new UnsupportedOperationException("LUTE can only be used to compute full-conformation energies");
	}

	@Override
	public ConfSearch.EnergiedConf calcEnergy(ConfSearch.ScoredConf conf) {
		return new ConfSearch.EnergiedConf(conf, calcEnergy(conf.getAssignments()));
	}

	@Override
	public ConfSearch.EnergiedConf calcEnergy(ConfSearch.ScoredConf conf, ResidueInteractions inters) {
		throw new UnsupportedOperationException("Not implemented yet... don't think anyone uses this anyway");
	}

	public double calcEnergy(int[] conf) {

		numCalculations.incrementAndGet();

		LUTE.LinearSystem system = lute.getTrainingSystem();

		// does this conf have a pair we didn't train on?
		for (RCTuple pair : Conf.getPairs(conf)) {
			if (!lute.getTuples().contains(pair)) {

				// yup, we can't give an energy for this conf
				return Double.NaN;
			}
		}

		// TODO: make a lookup table so we can do better than time linear in the number of tuples?
		// triples are sparse though, so hard to do lookups efficiently
		double energy = 0.0;
		for (int i=0; i<system.tuples.size(); i++) {
			RCTuple tuple = system.tuples.get(i);
			if (Conf.containsTuple(conf, tuple)) {
				energy += system.x.getEntry(i);
			}
		}

		return energy;
	}
}
