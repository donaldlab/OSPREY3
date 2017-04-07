package edu.duke.cs.osprey.energy;

import edu.duke.cs.osprey.confspace.ConfSearch.EnergiedConf;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.parallelism.TaskExecutor;
import edu.duke.cs.osprey.parallelism.TaskExecutor.TaskListener;

public interface ConfEnergyCalculator {
	
	EnergiedConf calcEnergy(ScoredConf conf);

	// use asynchronous techniques so we can parallelize conformation evaluation
	public static interface Async extends ConfEnergyCalculator {

		void calcEnergyAsync(ScoredConf conf, Listener listener);
		TaskExecutor getTasks();
		void cleanup();

		public static interface Listener extends TaskListener<EnergiedConf> {
			// nothing else to do
		}

		public static class Adapter implements Async {

			private ConfEnergyCalculator calc;
			private TaskExecutor tasks;

			public Adapter(ConfEnergyCalculator calc) {
				this.calc = calc;
				this.tasks = new TaskExecutor();
			}

			@Override
			public EnergiedConf calcEnergy(ScoredConf conf) {
				return calc.calcEnergy(conf);
			}

			@Override
			public void calcEnergyAsync(ScoredConf conf, Listener listener) {
				listener.onFinished(calc.calcEnergy(conf));
			}
			
			@Override
			public TaskExecutor getTasks() {
				return tasks;
			}

			@Override
			public void cleanup() {
				// nothing to do
			}
		}
	}
}
