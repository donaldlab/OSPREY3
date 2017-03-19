package edu.duke.cs.osprey.control;

import edu.duke.cs.osprey.confspace.ConfSearch.EnergiedConf;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;

public interface ConfEnergyCalculator {
	
	EnergiedConf calcEnergy(ScoredConf conf);

	// use asynchronous techniques so we can parallelize conformation evaluation
	public static interface Async extends ConfEnergyCalculator {

		void calcEnergyAsync(ScoredConf conf, Listener listener);
		int getParallelism();
		void waitForFinish();
		void cleanup();

		public static interface Listener {
			void onEnergy(EnergiedConf conf);
		}

		public static class Adapter implements Async {

			private ConfEnergyCalculator calc;

			public Adapter(ConfEnergyCalculator calc) {
				this.calc = calc;
			}

			@Override
			public int getParallelism() {
				return 1;
			}
			
			@Override
			public EnergiedConf calcEnergy(ScoredConf conf) {
				return calc.calcEnergy(conf);
			}

			@Override
			public void calcEnergyAsync(ScoredConf conf, Listener listener) {
				listener.onEnergy(calc.calcEnergy(conf));
			}
			
			@Override
			public void waitForFinish() {
				// nothing to do
			}

			@Override
			public void cleanup() {
				// nothing to do
			}
		}
	}
}
