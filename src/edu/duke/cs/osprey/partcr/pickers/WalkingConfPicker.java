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

package edu.duke.cs.osprey.partcr.pickers;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;

public class WalkingConfPicker implements ConfPicker {

	private int numItersPerConf;
	
	// NOTE: don't use int[] as the key value, primitive arrays doesn't have a good hashCode() method
	private Map<List<Integer>,Integer> pickedConfs;
	
	public WalkingConfPicker() {
		this(1);
	}
	
	public WalkingConfPicker(int numItersPerConf) {
		this.numItersPerConf = numItersPerConf;
		pickedConfs = new HashMap<>(); // yeah, ok to match on instance
	}

	@Override
	public ScoredConf pick(List<ScoredConf> confs) {
		
		for (ScoredConf conf : confs) {
			
			// use Java 8 magic to convert the array to a list
			// sadly, Arrays.toList() won't work here =(
			List<Integer> confList = Arrays.stream(conf.getAssignments()).boxed().collect(Collectors.toList());
			
			// have we picked this conf yet?
			Integer numIters = pickedConfs.get(confList);
			if (numIters == null) {
				
				// nope, let's start
				pickedConfs.put(confList, 1);
				return conf;
				
			} else if (numIters < numItersPerConf) {
				
				// yup, but not enough
				pickedConfs.put(confList, numIters + 1);
				return conf;
			}
			
			// yeah, but picked it too much already, so look to the next conf
		}
		
		// we ran out of confs, start over
		pickedConfs.clear();
		return pick(confs);
	}
}
