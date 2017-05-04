package edu.duke.cs.osprey.tools;

public class MathTools {
	
	public static int divUp(int num, int denom) {
		return (num + denom - 1)/denom;
	}
	
	public static int roundUpToMultiple(int val, int base) {
		int mod = val % base;
		if (mod == 0) {
			return val;
		}
		return val + base - mod;
	}
}
