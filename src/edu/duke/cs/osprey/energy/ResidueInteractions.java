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

package edu.duke.cs.osprey.energy;

import java.util.*;

import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.structure.Residues;
import edu.duke.cs.osprey.tools.HashCalculator;

/**
 * Represents interactions between residues (or within single residues) based on residue
 * numbers.
 * 
 * Does not rely on specific molecules, so residue interactions are portable
 * between different molecules, or molecule copies.
 * 
 * The {@link ResInterGen} helper class can assist with creating residue interactions.
 */
public class ResidueInteractions implements Iterable<ResidueInteractions.Pair> {

	private static String makeId(String resNum1, String resNum2) {

		// sort res numbers so we always have a stable order
		if (resNum1.compareTo(resNum2) > 0) {
			String swap = resNum1;
			resNum1 = resNum2;
			resNum2 = swap;
		}

		return resNum1 + " " + resNum2;
	}
	
	public static class Pair {
		
		public static final double IdentityWeight = 1;
		public static final double IdentityOffset = 0;
		
		public final String resNum1;
		public final String resNum2;
		public final double weight;
		public final double offset;

		public final String id;
		
		public Pair(String resNum1, String resNum2, double weight, double offset) {
			
			// sort res numbers so we always have a stable order
			if (resNum1.compareTo(resNum2) > 0) {
				String swap = resNum1;
				resNum1 = resNum2;
				resNum2 = swap;
			}
			
			this.resNum1 = resNum1;
			this.resNum2 = resNum2;
			this.weight = weight;
			this.offset = offset;

			this.id = makeId(resNum1, resNum2);
		}

		@Override
		public int hashCode() {
			return HashCalculator.combineHashes(
				resNum1.hashCode(),
				resNum2.hashCode()
			);
		}
		
		@Override
		public boolean equals(Object other) {
			if (other instanceof Pair) {
				return equals((Pair)other);
			}
			return false;
		}
		
		public boolean equals(Pair other) {
			return resNum1.equals(other.resNum1)
				&& resNum2.equals(other.resNum2);
		}

		public String getOtherResNum(String resNum) {
			if (resNum.equals(resNum1)) {
				return resNum2;
			} else if (resNum.equals(resNum2)) {
				return resNum1;
			} else {
				return null;
			}
		}
	}
	
	private Set<String> resNums;
	private Map<String,Pair> pairs;
	
	public ResidueInteractions() {
		resNums = new LinkedHashSet<>();
		pairs = new LinkedHashMap<>();
	}

	public ResidueInteractions(Pair ... pairs) {
		this();
		for (Pair pair : pairs) {
			addPair(pair.resNum1, pair.resNum2, pair.weight, pair.offset);
		}
	}

	public boolean contains(Pair pair) {
		return pairs.containsKey(pair.id);
	}

	public Pair get(String resNum1, String resNum2) {
		return pairs.get(makeId(resNum1, resNum2));
	}

	public void add(Pair pair) {
		resNums.add(pair.resNum1);
		resNums.add(pair.resNum2);
		pairs.put(pair.id, pair);
	}
	
	public void addSingle(String resNum) {
		addSingle(resNum, Pair.IdentityWeight, Pair.IdentityOffset);
	}
	
	public void addSingle(String resNum, double weight, double offset) {
		resNums.add(resNum);
		add(new Pair(resNum, resNum, weight, offset));
	}
	
	public void addPair(String resNum1, String resNum2) {
		addPair(resNum1, resNum2, Pair.IdentityWeight, Pair.IdentityOffset);
	}
	
	public void addPair(String resNum1, String resNum2, double weight, double offset) {
		resNums.add(resNum1);
		resNums.add(resNum2);
		add(new Pair(resNum1, resNum2, weight, offset));
	}

	public void addComplete(Residues residues) {
		addComplete(residues, Pair.IdentityWeight, Pair.IdentityOffset);
	}

	public void addComplete(Residues residues, double weight, double offset) {
		for (int i=0; i<residues.size(); i++) {
			Residue res1 = residues.get(i);
			addSingle(res1.getPDBResNumber(), weight, offset);
			for (int j=0; j<i; j++) {
				Residue res2 = residues.get(j);
				addPair(res1.getPDBResNumber(), res2.getPDBResNumber(), weight, offset);
			}
		}
	}
	
	public Set<String> getResidueNumbers() {
		return resNums;
	}
	
	@Override
	public Iterator<Pair> iterator() {
		return pairs.values().iterator();
	}
	
	public int size() {
		return pairs.size();
	}
	
	public Residues filter(Residues residues) {
		Residues filtered = new Residues();
		for (String resNum : resNums) {
			filtered.add(residues.getOrThrow(resNum));
		}
		return filtered;
	}

	public static ResidueInteractions subtract(ResidueInteractions a, ResidueInteractions b) {
		ResidueInteractions out = new ResidueInteractions();
		for (Pair pair : a) {
			if (!b.contains(pair)) {
				out.add(pair);
			}
		}
		return out;
	}
}
