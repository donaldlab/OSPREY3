package edu.duke.cs.osprey.energy;

import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.Set;

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
	
	public static class Pair {
		
		public static final double IdentityWeight = 1;
		public static final double IdentityOffset = 0;
		
		public final String resNum1;
		public final String resNum2;
		public final double weight;
		public final double offset;
		
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
	}
	
	private Set<String> resNums;
	private Set<Pair> pairs;
	
	public ResidueInteractions() {
		resNums = new LinkedHashSet<>();
		pairs = new LinkedHashSet<>();
	}

	public ResidueInteractions(Pair ... pairs) {
		this();
		for (Pair pair : pairs) {
			addPair(pair.resNum1, pair.resNum2, pair.weight, pair.offset);
		}
	}
	
	public void addSingle(String resNum) {
		addSingle(resNum, Pair.IdentityWeight, Pair.IdentityOffset);
	}
	
	public void addSingle(String resNum, double weight, double offset) {
		resNums.add(resNum);
		pairs.add(new Pair(resNum, resNum, weight, offset));
	}
	
	public void addPair(String resNum1, String resNum2) {
		addPair(resNum1, resNum2, Pair.IdentityWeight, Pair.IdentityOffset);
	}
	
	public void addPair(String resNum1, String resNum2, double weight, double offset) {
		resNums.add(resNum1);
		resNums.add(resNum2);
		pairs.add(new Pair(resNum1, resNum2, weight, offset));
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
		return pairs.iterator();
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
}
