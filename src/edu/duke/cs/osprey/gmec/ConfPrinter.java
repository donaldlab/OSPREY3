package edu.duke.cs.osprey.gmec;

import edu.duke.cs.osprey.confspace.ConfSearch.EnergiedConf;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;

public interface ConfPrinter {
	
	void print(SimpleConfSpace confSpace, EnergiedConf conf, EnergyWindow window);
	
	default void cleanup() {
		// nothing to do
	}
	
	default void print(SimpleConfSpace confSpace, EnergiedConf conf) {
		print(confSpace, conf, null);
	}
	
	public static class Nop implements ConfPrinter {

		@Override
		public void print(SimpleConfSpace confSpace, EnergiedConf conf, EnergyWindow window) {
			// do nothing, ie, a no-op
		}
	}
}
