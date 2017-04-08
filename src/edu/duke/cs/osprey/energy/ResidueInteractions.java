package edu.duke.cs.osprey.energy;

import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

import edu.duke.cs.osprey.tools.HashCalculator;

public class ResidueInteractions implements Iterable<ResidueInteractions.Pair> {
	
	public static class Pair {
		
		public static final double IdentityWeight = 1;
		
		public final String resNum1;
		public final String resNum2;
		public final double weight;
		
		public Pair(String resNum1, String resNum2, double weight) {
			
			// sort res numbers so we always have a stable order
			if (resNum1.compareTo(resNum2) > 0) {
				String swap = resNum1;
				resNum1 = resNum2;
				resNum2 = swap;
			}
			
			this.resNum1 = resNum1;
			this.resNum2 = resNum2;
			this.weight = weight;
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
	
	private Set<Pair> pairs;
	private double offset;
	
	public ResidueInteractions() {
		pairs = new HashSet<>();
		offset = 0;
	}
	
	public void addSingle(String resNum) {
		addSingle(resNum, Pair.IdentityWeight);
	}
	
	public void addSingle(String resNum, double weight) {
		pairs.add(new Pair(resNum, resNum, weight));
	}
	
	public void addPair(String resNum1, String resNum2) {
		addPair(resNum1, resNum2, Pair.IdentityWeight);
	}
	
	public void addPair(String resNum1, String resNum2, double weight) {
		pairs.add(new Pair(resNum1, resNum2, weight));
	}
	
	public double getOffset() {
		return offset;
	}
	
	public void addOffset(double val) {
		offset += val;
	}
	
	@Override
	public Iterator<Pair> iterator() {
		return pairs.iterator();
	}
}
