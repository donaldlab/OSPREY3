package edu.duke.cs.osprey.gmec;

import edu.duke.cs.osprey.confspace.ConfSearch.EnergiedConf;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.parallelism.TaskExecutor;
import edu.duke.cs.osprey.parallelism.TaskExecutor.TaskListener;

@Deprecated
public interface GMECConfEnergyCalculator {
	
	EnergiedConf calcEnergy(ScoredConf conf);

	// use asynchronous techniques so we can parallelize conformation evaluation
	public static interface Async extends GMECConfEnergyCalculator {

		void calcEnergyAsync(ScoredConf conf, Listener listener);
		TaskExecutor getTasks();

		public static interface Listener extends TaskListener<EnergiedConf> {
			// nothing else to do
		}

		public static class Adapter implements Async {

			private GMECConfEnergyCalculator calc;
			private TaskExecutor tasks;

			public Adapter(GMECConfEnergyCalculator calc) {
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
		}
	}
}
