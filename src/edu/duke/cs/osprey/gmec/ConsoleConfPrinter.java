package edu.duke.cs.osprey.gmec;

import edu.duke.cs.osprey.confspace.ConfSearch.EnergiedConf;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.SimpleConfSpace.Position;
import edu.duke.cs.osprey.confspace.SimpleConfSpace.ResidueConf;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

public class ConsoleConfPrinter implements ConfPrinter {
	
	private static final int LabelSize = 20;
	private static final String LabelFormat = "\t%-" + LabelSize + "s";
	
	@Override
	public void print(EnergiedConf conf, SimpleConfSpace confSpace, EnergyRange range) {
		System.out.print(makeReport(conf, confSpace, range));
	}
	
	public static String makeReport(EnergiedConf conf, SimpleConfSpace confSpace, EnergyRange range) {
		
		StringBuilder buf = new StringBuilder();
		
		if (confSpace != null) {

			buf.append(String.format(LabelFormat, "Residue numbers"));
			buf.append(" ");
			buf.append(confSpace.formatResidueNumbers());
			buf.append("\n");

			buf.append(String.format(LabelFormat, "Residue types"));
			buf.append(" ");
			buf.append(confSpace.formatConfSequence(conf));
			buf.append("\n");

			buf.append(String.format(LabelFormat, "Rotamer numbers"));
			buf.append(" ");
			buf.append(confSpace.formatConfRotamers(conf));
			buf.append("\n");
		}

		buf.append(String.format(LabelFormat, "Residue Conf Ids"));
		buf.append(" ");
		buf.append(SimpleConfSpace.formatConfRCs(conf));
		buf.append("\n");

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
	
	
	public static HashMap<String,List<String>> makeReportMap(EnergiedConf conf, SimpleConfSpace confSpace, EnergyRange range) {
		
                HashMap<String,List<String>> map = new HashMap<>();
                
                ArrayList<String> confRCs = new ArrayList<>();
                for(int rc : conf.getAssignments())
                    confRCs.add(String.valueOf(rc));
                map.put("CONF", confRCs);
                
		
		if (confSpace != null) {
                        
                        ArrayList<String> seq = new ArrayList<>();
			for (Position pos : confSpace.positions) {
				ResidueConf resConf = pos.resConfs.get(conf.getAssignments()[pos.index]);
				seq.add(resConf.template.name);
			}
			map.put("SEQ", seq);
                        
                        ArrayList<String> rots = new ArrayList<>();
			for (Position pos : confSpace.positions) {
				ResidueConf resConf = pos.resConfs.get(conf.getAssignments()[pos.index]);
				rots.add(resConf.getRotamerCode());
			}
                        map.put("ROTS", rots);
                }
                
                map.put("ENERGY", Arrays.asList(String.format("%.2f",conf.getEnergy())));//DEBUG!!
                map.put("SCORE", Arrays.asList(String.valueOf(conf.getScore())));
                if(range!=null){
                    map.put("BESTENERGY", Arrays.asList(String.valueOf(range.getMin())));
                    map.put("BESTSCORE", Arrays.asList(String.valueOf(range.getMax() - conf.getScore())));
                }
                
		return map;
	}

}
