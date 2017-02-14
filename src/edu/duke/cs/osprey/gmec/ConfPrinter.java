package edu.duke.cs.osprey.gmec;

import edu.duke.cs.osprey.confspace.ConfSearch.EnergiedConf;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;

public interface ConfPrinter {
	
	void print(EnergiedConf conf, SimpleConfSpace confSpace, EnergyRange range);
	
	default void print(EnergiedConf conf, SimpleConfSpace confSpace) {
		print(conf, confSpace, null);
	}
	
	default void print(EnergiedConf conf) {
		print(conf, null, null);
	}
	
	default void cleanup() {
		// nothing to do
	}
	
	public static class Nop implements ConfPrinter {

		@Override
		public void print(EnergiedConf conf, SimpleConfSpace confSpace, EnergyRange range) {
			// do nothing, ie, a no-op
		}
	}
}
