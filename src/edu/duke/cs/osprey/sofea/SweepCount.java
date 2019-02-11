package edu.duke.cs.osprey.sofea;

import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;


public class SweepCount implements Sofea.Criterion {

	public final int sweepCount;

	public SweepCount(int sweepCount) {
		this.sweepCount = sweepCount;
	}

	@Override
	public Satisfied isSatisfied(SeqDB seqdb, FringeDB fringedb, long sweepCount, BoltzmannCalculator bcalc) {
		return sweepCount >= this.sweepCount ? Satisfied.Terminate : Satisfied.KeepIterating;
	}
}
