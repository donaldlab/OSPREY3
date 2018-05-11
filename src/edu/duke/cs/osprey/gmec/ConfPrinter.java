package edu.duke.cs.osprey.gmec;

import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;

public interface ConfPrinter {
	
	void print(ConfSearch.EnergiedConf conf, SimpleConfSpace confSpace, EnergyRange range);
	
	default void print(ConfSearch.EnergiedConf conf, SimpleConfSpace confSpace) {
		print(conf, confSpace, null);
	}

	default void print(ConfSearch.ScoredConf conf, SimpleConfSpace confSpace) {
		print(new ConfSearch.EnergiedConf(conf, Double.NaN), confSpace, null);
	}
	
	default void print(ConfSearch.EnergiedConf conf) {
		print(conf, null, null);
	}

	default void print(ConfSearch.ScoredConf conf) {
		print(new ConfSearch.EnergiedConf(conf, Double.NaN), null, null);
	}

	default void cleanup() {
		// nothing to do
	}
	
	public static class Nop implements ConfPrinter {

		@Override
		public void print(ConfSearch.EnergiedConf conf, SimpleConfSpace confSpace, EnergyRange range) {
			// do nothing, ie, a no-op
		}
	}
}
