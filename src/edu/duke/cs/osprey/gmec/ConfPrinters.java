package edu.duke.cs.osprey.gmec;

import java.util.Arrays;
import java.util.List;

import edu.duke.cs.osprey.confspace.ConfSearch.EnergiedConf;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;

public class ConfPrinters implements ConfPrinter {
	
	public final List<ConfPrinter> printers;
	
	public ConfPrinters(ConfPrinter ... printers) {
		this.printers = Arrays.asList(printers);
	}

	@Override
	public void print(EnergiedConf conf, SimpleConfSpace confSpace, EnergyRange window) {
		for (ConfPrinter printer : printers) {
			printer.print(conf, confSpace, window);
		}
	}
}
