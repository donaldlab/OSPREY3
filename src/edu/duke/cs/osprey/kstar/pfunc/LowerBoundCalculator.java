package edu.duke.cs.osprey.kstar.pfunc;

import edu.duke.cs.osprey.confspace.ConfDB;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.parallelism.TaskExecutor;
import edu.duke.cs.osprey.tools.MathTools;
import edu.duke.cs.osprey.tools.TimeTools;

import java.math.BigDecimal;


public class LowerBoundCalculator {

	public static enum Status {
		HasLowEnergies,
		OutOfConfs,
		OutOfLowEnergies
	}

	public final ConfSearch tree;
	public final ConfEnergyCalculator ecalc;

	public BigDecimal weightedScoreSum = BigDecimal.ZERO;
	public BigDecimal weightedEnergySum = BigDecimal.ZERO;
	public int numConfsScored = 0;
	public int numConfsEnergied = 0;
	public ConfDB confDB = null;
	public ConfDB.ConfTable confTable = null;

	private BoltzmannCalculator boltzmann = new BoltzmannCalculator(PartitionFunction.decimalPrecision);

	public LowerBoundCalculator(ConfSearch tree, ConfEnergyCalculator ecalc) {
		this.tree = tree;
		this.ecalc = ecalc;
	}

	public Status run(int numConfs) {
		return run(numConfs, null);
	}

	public Status run(int numConfs, TaskExecutor.TaskListener<ConfSearch.EnergiedConf> confListener) {
		Status status = Status.HasLowEnergies;
		for (int i=0; i<numConfs && status == Status.HasLowEnergies; i++) {
			status = energyNextConfAsync(confListener);
		}
		waitForFinish();
		return status;
	}

	public Status energyNextConfAsync() {
		return energyNextConfAsync(null);
	}

	public Status energyNextConfAsync(TaskExecutor.TaskListener<ConfSearch.EnergiedConf> confListener) {

		ConfSearch.ScoredConf conf = tree.nextConf();
		if (conf == null) {
			return Status.OutOfConfs;
		}
		numConfsScored++;
		if (Double.isInfinite(conf.getScore()) || MathTools.isZero(boltzmann.calc(conf.getScore()))) {
			return Status.OutOfLowEnergies;
		}

		// what to do when we get a conf energy?
		TaskExecutor.TaskListener<ConfSearch.EnergiedConf> onEconf = (ConfSearch.EnergiedConf econf) -> {

			// we might be on the listener thread, so sync to keep from racing the main thread
			synchronized (LowerBoundCalculator.this) {

				// energy calculation done, update pfunc values
				if (!Double.isInfinite(econf.getEnergy())) {
					weightedEnergySum = weightedEnergySum.add(boltzmann.calc(econf.getEnergy()));
				}
				weightedScoreSum = weightedScoreSum.add(boltzmann.calc(econf.getScore()));
				numConfsEnergied++;
			}

			if (confListener != null) {
				confListener.onFinished(econf);
			}
		};

		// do we have a db table for this conf?
		ConfDB.ConfTable confTable;
		if (this.confTable != null) {
			confTable = this.confTable;
		} else if (this.confDB != null) {
			Sequence sequence = confDB.confSpace.makeSequenceFromAssignments(conf.getAssignments());
			confTable = confDB.getSequence(sequence);
		} else {
			confTable = null;
		}

		// is this conf in the db?
		if (confTable != null) {
			ConfDB.Conf dbconf = confTable.get(conf.getAssignments());
			if (dbconf != null) {

				// yup, don't calculate the energy again
				// (but make sure the energy gets reported on the listener thread
				// so we behave exactly the same as if we did the energy calculation)
				ecalc.tasks.submit(
					() -> new ConfSearch.EnergiedConf(conf, dbconf.upper.energy),
					(econf) -> onEconf.onFinished(econf)
				);

				return Status.HasLowEnergies;
			}
		}

		// nope, do the energy calculation asynchronously
		ecalc.calcEnergyAsync(conf, (econf) -> {

			// save the conf in the db if needed
			if (confTable != null) {
				confTable.setBounds(
					conf.getAssignments(),
					conf.getScore(),
					econf.getEnergy(),
					TimeTools.getTimestampNs()
				);
				confTable.flush();
			}

			onEconf.onFinished(econf);
		});

		return Status.HasLowEnergies;
	}

	public void waitForFinish() {
		ecalc.tasks.waitForFinish();
	}

	@Override
	public String toString() {
		return String.format("LowerBoundCalculator   scored: %6d   energied: %6d   energies: %e   scores: %e",
			numConfsScored,
			numConfsEnergied,
			weightedEnergySum.doubleValue(),
			weightedScoreSum.doubleValue()
		);
	}
}
