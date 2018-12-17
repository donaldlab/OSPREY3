package edu.duke.cs.osprey.sofea;


public class SweepCount implements Sofea.Criterion {

	public final int sweepCount;

	public SweepCount(int sweepCount) {
		this.sweepCount = sweepCount;
	}

	@Override
	public boolean isFinished(SeqDB seqdb, FringeDB fringedb, long sweepCount) {
		return sweepCount >= this.sweepCount;
	}
}
