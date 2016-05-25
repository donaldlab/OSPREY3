package edu.duke.cs.osprey.tools;

import java.util.ArrayList;

public class Combinations {

	public static ArrayList<int[]> CombinationsIterative ( int[] in, int k ) {

		ArrayList<int[]> subsets = new ArrayList<>();

		int[] s = new int[k]; // here we'll keep indices pointing to elements in input array

		if (k <= in.length) {
			// first index sequence: 0, 1, 2, ...

			for (int i = 0; (s[i] = i) < k - 1; i++);  

			subsets.add(getSubset(in, s));

			for(;;) {

				int i;
				// find position of item that can be incremented
				for (i = k - 1; i >= 0 && s[i] == in.length - k + i; i--); 

				if (i < 0) {
					break;
				} 

				else {
					s[i]++; // increment this item

					for (++i; i < k; i++) { // fill up remaining items
						s[i] = s[i - 1] + 1; 
					}

					subsets.add(getSubset(in, s));
				}
			}
		}

		return subsets;
	}


	// generate actual subset by index sequence
	private static int[] getSubset(int[] input, int[] subset) {

		int[] result = new int[subset.length]; 

		for (int i = 0; i < subset.length; i++) 
			result[i] = input[subset[i]];

		return result;
	}

}
