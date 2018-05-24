/*
** This file is part of OSPREY 3.0
** 
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
** 
** OSPREY relies on grants for its development, and since visibility
** in the scientific literature is essential for our success, we
** ask that users of OSPREY cite our papers. See the CITING_OSPREY
** document in this distribution for more information.
** 
** Contact Info:
**    Bruce Donald
**    Duke University
**    Department of Computer Science
**    Levine Science Research Center (LSRC)
**    Durham
**    NC 27708-0129
**    USA
**    e-mail: www.cs.duke.edu/brd/
** 
** <signature of Bruce Donald>, Mar 1, 2018
** Bruce Donald, Professor of Computer Science
*/

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
