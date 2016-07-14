package edu.duke.cs.osprey.tools;

public class HashCalculator {
	
	public static int combineHashes(long ... nums) {
		long hashCode = 1;
		for (long i : nums) {
			hashCode = hashCode * 31 + i;
		}
		return (int)((hashCode >> 32) ^ hashCode);
	}
	
	public static int combineHashes(int ... nums) {
		int hashCode = 1;
		for (int i : nums) {
			hashCode = hashCode * 31 + i;
		}
		return hashCode;
	}
	
	public static int combineHashes(short ... nums) {
		int hashCode = 1;
		for (int i : nums) {
			hashCode = hashCode * 31 + i;
		}
		return hashCode;
	}
	
	public static int combineHashes(byte ... nums) {
		int hashCode = 1;
		for (int i : nums) {
			hashCode = hashCode * 31 + i;
		}
		return hashCode;
	}
	
	public static int combineHashesCommutative(int ... nums) {
		int hashCode = 1;
		for (int i : nums) {
			hashCode += i;
		}
		return hashCode;
	}
	
	public static int hashIds(int ... nums) {
		int hashCode = 1;
		for (int i : nums) {
			hashCode = hashCode * 37 ^ (i + 1);
		}
		return hashCode;
	}
}
