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

package edu.duke.cs.osprey.structure.analysis;


import java.util.Arrays;
import java.util.Collection;

public class Histogram {

	public final int numData;
	public final int[] counts;
	public final double[] boundaries;

	public Histogram(Collection<Double> data, int numBuckets) {

		if (data.size() <= 0) {
			throw new IllegalArgumentException("no data");
		}
		if (numBuckets <= 0) {
			throw new IllegalArgumentException("wrong number of buckets: " + numBuckets);
		}

		numData = data.size();

		// get the data min,max
		double min = Double.POSITIVE_INFINITY;
		double max = Double.NEGATIVE_INFINITY;
		for (double d : data) {
			min = Math.min(min, d);
			max = Math.max(max, d);
		}

		// make the buckets
		counts = new int[numBuckets];
		boundaries = new double[numBuckets + 1];
		boundaries[0] = min;
		boundaries[numBuckets] = max;
		for (int i=1; i<numBuckets; i++) {
			boundaries[i] = min + (max - min)*i/numBuckets;
		}

		// count the data
		Arrays.fill(counts, 0);
		for (double d : data) {
			counts[findBucket(d)]++;
		}
	}

	private int findBucket(double d) {

		for (int i=1; i<boundaries.length; i++) {
			if (d < boundaries[i]) {
				return i - 1;
			}
		}

		return counts.length - 1;
	}

	@Override
	public String toString() {
		return toString(4, 1, 20);
	}

	public String toString(int size, int precision, int barLength) {

		int maxCount = 0;
		int countLength = 0;
		for (int count : counts) {
			maxCount = Math.max(maxCount, count);
			countLength = Math.max(countLength, Integer.toString(count).length());
		}

		String dataFormat = "%" + size + "." + precision + "f";
		String bucketFormat = "Bucket [" + dataFormat + "," + dataFormat + "] %" + countLength + "d  ";

		StringBuilder buf = new StringBuilder();
		for (int i=0; i<counts.length; i++) {

			if (i > 0) {
				buf.append('\n');
			}

			buf.append(String.format(bucketFormat, boundaries[i], boundaries[i+1], counts[i]));

			int numSegments = counts[i]*barLength/maxCount;
			for (int j=0; j<barLength; j++) {
				if (j < numSegments) {
					buf.append('-');
				} else {
					buf.append(' ');
				}
			}
		}
		return buf.toString();
	}
}
