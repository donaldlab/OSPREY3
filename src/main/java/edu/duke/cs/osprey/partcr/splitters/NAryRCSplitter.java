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

package edu.duke.cs.osprey.partcr.splitters;

import java.util.ArrayList;
import java.util.List;

import edu.duke.cs.osprey.confspace.RC;

public class NAryRCSplitter extends RCSplitter {
	
	private static enum Side {
		
		Low,
		High;
		
		public static Side get(int rcIndex, int dofIndex) {
			// whats the nth bit set to?
			int ord = (rcIndex >> dofIndex) & 0x1;
			return values()[ord];
		}
	}
	
	@Override
	public List<RC> split(int pos, RC rc) {
		
		List<RC> rcs = new ArrayList<>();
		
		// get the splittable DoFs
		List<Integer> indices = new ArrayList<>();
		for (int i=0; i<rc.DOFs.size(); i++) {
			if (rc.DOFmin.get(i) != rc.DOFmax.get(i)) {
				indices.add(i);
			}
		}
		int numDofs = indices.size();
		
		// split all the dimensions at once
		int numRcs = 1 << numDofs;
		for (int rcIndex=0; rcIndex<numRcs; rcIndex++) {
			
			// clone the RC
			ArrayList<Double> mins = new ArrayList<>(rc.DOFmin);
			ArrayList<Double> maxs = new ArrayList<>(rc.DOFmax);
			RC subRc = new RC(rc.AAType, rc.template, rc.rotNum, rc.DOFs, mins, maxs, rc.RCIndex);
			
			// change the dof bounds
			for (int dofIndex : indices) {
				
				// split the dof down the middle
				double min = rc.DOFmin.get(dofIndex);
				double max = rc.DOFmax.get(dofIndex);
				double mid = (max + min)/2;
				
				switch (Side.get(rcIndex, dofIndex)) {
					case Low:
						mins.set(dofIndex, min);
						maxs.set(dofIndex, mid);
					break;
					case High:
						mins.set(dofIndex, mid);
						maxs.set(dofIndex, max);
					break;
				}
			}
			
			rcs.add(subRc);
		}
		
		// TEMP: check the split RCs
		double volume = calcVolume(rc, indices);
		for (RC subRc : rcs) {
			volume -= calcVolume(subRc, indices);
		}
		assert (Math.abs(volume) < 1e-10);
		
		return rcs;
	}
	
	private double calcVolume(RC rc, List<Integer> indices) {
		double volume = 1;
		for (int dofIndex : indices) {
			double min = rc.DOFmin.get(dofIndex);
			double max = rc.DOFmax.get(dofIndex);
			volume *= max - min;
		}
		return volume;
	}
}
