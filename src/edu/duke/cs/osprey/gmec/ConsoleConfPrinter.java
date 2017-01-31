package edu.duke.cs.osprey.gmec;

import edu.duke.cs.osprey.confspace.ConfSearch.EnergiedConf;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.SimpleConfSpace.Position;
import edu.duke.cs.osprey.confspace.SimpleConfSpace.ResidueConf;

public class ConsoleConfPrinter implements ConfPrinter {
	
	private static final int LabelSize = 30;
	private static final String LabelFormat = "\t%-" + LabelSize + "s";
	
	public static String makeReport(SimpleConfSpace confSpace, EnergiedConf conf) {
		return makeReport(confSpace, conf, null);
	}

	public static String makeReport(SimpleConfSpace confSpace, EnergiedConf conf, EnergyWindow window) {
		
		StringBuilder buf = new StringBuilder();
		
		buf.append(makeConfReport(confSpace, conf.getAssignments()));
		
		buf.append(String.format(LabelFormat + " %.6f", "Energy", conf.getEnergy()));
		if (window != null) {
			buf.append(String.format(" (best so far: %.6f)", window.getMin()));
		}
		buf.append("\n");
		
		buf.append(String.format(LabelFormat + " %.6f (gap: %.6f", "Score", conf.getScore(), Math.abs(conf.getScore() - conf.getEnergy())));
		if (window != null) {
			buf.append(String.format(", remaining: %.6f", window.getMax() - conf.getScore()));
		}
		buf.append(")\n");
		
		return buf.toString();
	}
	
	public static String makeConfReport(SimpleConfSpace confSpace, int[] conf) {
		
		StringBuilder buf = new StringBuilder();
	
		buf.append(String.format(LabelFormat, "RCs (residue-based numbers)"));
		for (int rc : conf) {
			buf.append(String.format(" %3d", rc));
		}
		buf.append("\n");

		buf.append(String.format(LabelFormat, "Residue types"));
		for (Position pos : confSpace.positions) {
			ResidueConf resConf = pos.resConfs.get(conf[pos.index]);
			buf.append(String.format(" %3s", resConf.template.name));
		}
		buf.append("\n");

		buf.append(String.format(LabelFormat, "Rotamer numbers"));
		for (Position pos : confSpace.positions) {
			ResidueConf resConf = pos.resConfs.get(conf[pos.index]);
			buf.append(String.format(" %3s", resConf.getRotamerCode()));
		}
		buf.append("\n");
		
		return buf.toString();
	}

	@Override
	public void print(SimpleConfSpace confSpace, EnergiedConf conf, EnergyWindow window) {
		System.out.print(makeReport(confSpace, conf, window));
	}
}
