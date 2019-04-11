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

import edu.duke.cs.osprey.tools.MathTools;

import java.util.Arrays;
import java.util.TreeSet;
import java.util.function.Function;


/** a cluster of data points in R^1 */
public class ClusterR1 {

	private final TreeSet<Double> values = new TreeSet<>();

	public void add(double a) {
		values.add(a);
	}

	public int size() {
		return values.size();
	}

	public boolean isEmpty() {
		return values.isEmpty();
	}


	public class Stats {

		public final ClusterR1 cluster = ClusterR1.this;

		public final MathTools.DoubleBounds bounds = bounds();
		public final double mean = mean();

		public final double mode;

		public Stats(double modeIncludesPercent, int modeBuckets) {
			mode = mode(getInterval(mean, modeIncludesPercent), modeBuckets);
		}

		public MathTools.DoubleBounds getInterval(double center, double includePercent) {

			// count how many points are on each side of the center
			int numNeg = 0;
			int numPos = 0;
			for (double a : values) {
				if (a < center) {
					numNeg++;
				} else if (a > center) {
					numPos++;
				}
			}

			// figure out how many samples to skip for each side
			double skipRatio = 1.0 - includePercent/100.0;
			int numNegToSkip = (int)(numNeg*skipRatio);
			int numPosToSkip = (int)(numPos*skipRatio);

			MathTools.DoubleBounds interval = new MathTools.DoubleBounds(center, center);

			// do the negative side
			int count = 0;
			for (double a : values) {
				if (count++ > numNegToSkip) {
					interval.lower = a;
					break;
				}
			}

			// do the positive side
			count = 0;
			for (double a : values.descendingSet()) {
				if (count++ > numPosToSkip) {
					interval.upper = a;
					break;
				}
			}

			return interval;
		}
	}

	public MathTools.DoubleBounds bounds() {
		MathTools.DoubleBounds bounds = new MathTools.DoubleBounds(values.first(), values.first());
		for (double a : values) {
			bounds.expand(a);
		}
		return bounds;
	}

	public double mean() {
		double sum = 0;
		for (double a : values) {
			 sum += a;
		}
		return sum/values.size();
	}

	/** use a histogram with a certain number of buckets to find the mode */
	public double mode(MathTools.DoubleBounds bounds, int numBuckets) {

		double min = bounds.lower;
		double width = bounds.size();

		Function<Double,Integer> valueToBucket = v -> (int)((v - min)*numBuckets/width);
		Function<Integer,Double> bucketToValue = i -> width*i/numBuckets + min;

		// build the histogram
		int[] counts = new int[numBuckets];
		Arrays.fill(counts, 0);
		for (double v : values) {
			int i = valueToBucket.apply(v);
			if (i >= 0 && i < numBuckets) {
				counts[i]++;
			}
		}

		// find the biggest bucket
		int bucket = 0;
		int maxCount = 0;
		for (int i=0; i<numBuckets; i++) {
			if (counts[i] > maxCount) {
				maxCount = counts[i];
				bucket = i;
			}
		}

		return bucketToValue.apply(bucket);
	}
}
