package edu.duke.cs.osprey.sofea;

import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;


public class StepCount implements Sofea.Criterion {

	public final int maxPass1Steps;
	public final int maxPass2Steps;

	public StepCount(int maxPass1Steps, int maxPass2Steps) {
		this.maxPass1Steps = maxPass1Steps;
		this.maxPass2Steps = maxPass2Steps;
	}

	@Override
	public Satisfied isSatisfied(SeqDB seqdb, FringeDB fringedbLower, FringeDB fringedbUpper, long pass1Step, long pass2Step, BoltzmannCalculator bcalc) {
		return pass1Step <= maxPass1Steps && pass2Step <= maxPass2Steps ? Satisfied.Terminate : Satisfied.KeepSweeping;
	}
}
