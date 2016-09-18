package edu.duke.cs.osprey.control;

import edu.duke.cs.osprey.confspace.ConfSearch.EnergiedConf;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;

public interface ConfEnergyCalculator {
	
	EnergiedConf calcEnergy(ScoredConf conf);

	// use asynchronous techniques so we can parallelize conformation evaluation
	public static interface Async extends ConfEnergyCalculator {

		void setListener(Listener listener);
		void calcEnergyAsync(ScoredConf conf);
		void waitForFinish();
		void cleanup();

		public static interface Listener {
			void onEnergy(EnergiedConf conf);
		}

		public static class Adapter implements Async {

			private ConfEnergyCalculator calc;
			private Listener listener;

			public Adapter(ConfEnergyCalculator calc) {
				this.calc = calc;
			}

			@Override
			public EnergiedConf calcEnergy(ScoredConf conf) {
				return calc.calcEnergy(conf);
			}

			@Override
			public void setListener(Listener listener) {
				this.listener = listener;
			}

			@Override
			public void calcEnergyAsync(ScoredConf conf) {
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
