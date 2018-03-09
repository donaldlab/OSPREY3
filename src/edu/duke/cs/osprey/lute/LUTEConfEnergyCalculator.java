package edu.duke.cs.osprey.lute;

import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.ResidueInteractions;


public class LUTEConfEnergyCalculator extends ConfEnergyCalculator {

	public final EnergyMatrix lutemat;

	public LUTEConfEnergyCalculator(SimpleConfSpace confSpace, EnergyCalculator ecalc, EnergyMatrix lutemat) {
		super(confSpace, ecalc, null, null, false);

		this.lutemat = lutemat;
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
		// TODO: add energy contribution of triples
		double energy = lutemat.getInternalEnergy(new RCTuple(conf.getAssignments()));
		return new ConfSearch.EnergiedConf(conf, energy);
	}

	@Override
	public ConfSearch.EnergiedConf calcEnergy(ConfSearch.ScoredConf conf, ResidueInteractions inters) {
		throw new UnsupportedOperationException("Not implemented yet... don't think anyone uses this anyway");
	}
}
