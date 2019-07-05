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

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import org.apache.commons.io.IOUtils;

import edu.duke.cs.osprey.confspace.ConfSearch.EnergiedConf;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.SimpleConfSpace.Position;
import edu.duke.cs.osprey.confspace.SimpleConfSpace.ResidueConf;

public class LoggingConfPrinter implements ConfPrinter {
	
	private File file;
	private FileWriter fout;
	private int numConfs;
	private double minEnergy;
	
	public LoggingConfPrinter(File file) {
		
		// open the log file
		this.file = file;
		try {
			fout = new FileWriter(file);
		} catch (IOException ex) {
			throw new RuntimeException("can't open log file at: " + file.getAbsolutePath(), ex);
		}
		
		numConfs = 0;
		minEnergy = Double.POSITIVE_INFINITY;
	}
	
	@Override
	public void print(EnergiedConf conf, SimpleConfSpace confSpace, EnergyRange range) {
		
		// track min energies
		minEnergy = Math.min(minEnergy, conf.getEnergy());
		
		// log the conformation
		try {
			
			// write the conf
			fout.write(Integer.toString(numConfs));
			fout.write(" CONF:");
			for (int rc : conf.getAssignments()) {
				fout.write(" ");
				fout.write(Integer.toString(rc));
			}
			
			if (confSpace != null) {
			
				// write the residue types
				fout.write(" RESTYPES:");
				for (Position pos : confSpace.positions) {
					ResidueConf resConf = pos.resConfs.get(conf.getAssignments()[pos.index]);
					fout.write(" ");
					fout.write(resConf.template.name);
				}
				
				// write the rotamers
				fout.write(" ROTS:");
				for (Position pos : confSpace.positions) {
					ResidueConf resConf = pos.resConfs.get(conf.getAssignments()[pos.index]);
					fout.write(" ");
					fout.write(resConf.getRotamerCode());
				}
			}
			
			// write the energies
			fout.write(" Score: ");
			fout.write(Double.toString(conf.getScore()));
			fout.write(" Energy: ");
			fout.write(Double.toString(conf.getEnergy()));
			fout.write(" Best so far: ");
			fout.write(Double.toString(minEnergy));
			
			fout.write("\n");
			
			// flush writes to disk immediately so we save as much info as possible after a failure
			fout.flush();
			
		} catch (IOException ex) {
			throw new RuntimeException("can't write conf to log file at: " + file.getAbsolutePath(), ex);
		}
		
		numConfs++;
	}
	
	@Override
	public void cleanup() {
		IOUtils.closeQuietly(fout);
	}
}
