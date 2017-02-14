package edu.duke.cs.osprey.gmec;

import edu.duke.cs.osprey.confspace.ConfSearch.EnergiedConf;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.SimpleConfSpace.Position;
import edu.duke.cs.osprey.confspace.SimpleConfSpace.ResidueConf;

public class ConsoleConfPrinter implements ConfPrinter {
	
	private static final int LabelSize = 20;
	private static final String LabelFormat = "\t%-" + LabelSize + "s";
	
	@Override
	public void print(EnergiedConf conf, SimpleConfSpace confSpace, EnergyRange range) {
		System.out.print(makeReport(conf, confSpace, range));
	}
	
	public static String makeReport(EnergiedConf conf, SimpleConfSpace confSpace, EnergyRange range) {
		
		StringBuilder buf = new StringBuilder();
		
		buf.append(String.format(LabelFormat, "Residue Conf Ids"));
		for (int rc : conf.getAssignments()) {
			buf.append(String.format(" %3d", rc));
		}
		buf.append("\n");
		
		if (confSpace != null) {

			buf.append(String.format(LabelFormat, "Residue types"));
			for (Position pos : confSpace.positions) {
				ResidueConf resConf = pos.resConfs.get(conf.getAssignments()[pos.index]);
				buf.append(String.format(" %3s", resConf.template.name));
			}
			buf.append("\n");

			buf.append(String.format(LabelFormat, "Rotamer numbers"));
			for (Position pos : confSpace.positions) {
				ResidueConf resConf = pos.resConfs.get(conf.getAssignments()[pos.index]);
				buf.append(String.format(" %3s", resConf.getRotamerCode()));
			}
			buf.append("\n");
		}
		
		buf.append(String.format(LabelFormat + " %.6f", "Energy", conf.getEnergy()));
		if (range != null) {
			buf.append(String.format(" (best so far: %.6f)", range.getMin()));
		}
		buf.append("\n");
		
		buf.append(String.format(LabelFormat + " %.6f (gap: %.6f", "Score", conf.getScore(), Math.abs(conf.getScore() - conf.getEnergy())));
		if (range != null) {
			buf.append(String.format(", remaining: %.6f", range.getMax() - conf.getScore()));
		}
		buf.append(")\n");
		
		return buf.toString();
	}
}
