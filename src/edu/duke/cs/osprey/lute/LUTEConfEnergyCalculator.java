package edu.duke.cs.osprey.lute;

import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.TuplesIndex;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.ResidueInteractions;
import edu.duke.cs.osprey.minimization.MoleculeObjectiveFunction;
import edu.duke.cs.osprey.parallelism.TaskExecutor;


public class LUTEConfEnergyCalculator extends ConfEnergyCalculator {

	public final LUTEState state;
	public final TuplesIndex tuples;

	public LUTEConfEnergyCalculator(SimpleConfSpace confSpace, LUTEState state) {
		super(confSpace, new TaskExecutor()); // TODO: parallelism?

		this.state = state;
		this.tuples = new TuplesIndex(confSpace, state.tuples);
	}

	private static class NotSupportedByLUTEException extends RuntimeException {
		public NotSupportedByLUTEException() {
			super("LUTE can only be used to compute full-conformation energies");
		}
	}

	@Override
	public ResidueInteractions makeSingleInters(int pos, int rc) {
		throw new NotSupportedByLUTEException();
	}

	@Override
	public ResidueInteractions makePairInters(int pos1, int rc1, int pos2, int rc2) {
		throw new NotSupportedByLUTEException();
	}

	@Override
	public ConfSearch.EnergiedConf calcEnergy(ConfSearch.ScoredConf conf) {
		return new ConfSearch.EnergiedConf(conf, calcEnergy(conf.getAssignments()));
	}

	@Override
	public ConfSearch.EnergiedConf calcEnergy(ConfSearch.ScoredConf conf, ResidueInteractions inters) {
		throw new UnsupportedOperationException("Not implemented yet... don't think anyone uses this anyway");
	}

	@Override
	public EnergyCalculator.EnergiedParametricMolecule calcEnergy(RCTuple frag, ResidueInteractions inters) {
		throw new NotSupportedByLUTEException();
	}

	@Override
	public MoleculeObjectiveFunction makeIntraShellObjFcn(int pos, int rc) {
		throw new NotSupportedByLUTEException();
	}

	@Override
	public MoleculeObjectiveFunction makePairwiseObjFcn(int pos1, int rc1, int pos2, int rc2) {
		throw new NotSupportedByLUTEException();
	}

	public double calcEnergy(int[] conf) {

		numCalculations.incrementAndGet();

		final boolean throwIfMissingSingle = false; // we're not fitting singles
		final boolean throwIfMissingPair = true; // we always fit to dense pairs, confs shouldn't be using pruned pairs

		// silly Java... If this were Kotlin, we wouldn't have to use an array to make the energy modifiable
		// hopefully JVM escape analysis will stack-allocate this?
		final double[] energy = new double[] { 0.0 };
		tuples.forEachIn(conf, throwIfMissingSingle, throwIfMissingPair, (t) -> {
			energy[0] += state.tupleEnergies[t];
		});
		return energy[0] + state.tupleEnergyOffset;
	}

	public boolean hasTuple(int pos, int rc) {
		return tuples.getIndex(pos, rc) != null;
	}

	public boolean hasTuple(int pos1, int rc1, int pos2, int rc2) {
		return tuples.getIndex(pos1, rc1, pos2, rc2) != null;
	}

	public boolean hasTuple(int pos1, int rc1, int pos2, int rc2, int pos3, int rc3) {
		return tuples.getIndex(pos1, rc1, pos2, rc2, pos3, rc3) != null;
	}

	public double getEnergy(int pos, int rc) {
		return getEnergy(tuples.getIndex(pos, rc));
	}

	public double getEnergy(int pos1, int rc1, int pos2, int rc2) {
		return getEnergy(tuples.getIndex(pos1, rc1, pos2, rc2));
	}

	public double getEnergy(int pos1, int rc1, int pos2, int rc2, int pos3, int rc3) {
		return getEnergy(tuples.getIndex(pos1, rc1, pos2, rc2, pos3, rc3));
	}

	private double getEnergy(Integer index) {
		if (index == null) {
			return 0.0;
		}
		return state.tupleEnergies[index];
	}
}
