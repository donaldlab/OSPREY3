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

import edu.duke.cs.osprey.tools.Protractor;

import java.util.Arrays;
import java.util.TreeSet;
import java.util.function.Function;


/** a cluster of data points in S^1 */
public class ClusterS1 {

	private final TreeSet<Double> angles = new TreeSet<>();

	public void add(double a) {
		angles.add(a);
	}


	public class Stats {

		public final ClusterS1 cluster = ClusterS1.this;

		public final SmallAngleVoxel.Interval bounds = bounds();
		public final double mean = mean(bounds);

		public final double mode;

		public Stats(double modeIncludesPercent, int modeBuckets) {
			mode = mode(getInterval(mean, modeIncludesPercent), modeBuckets);
		}

		public SmallAngleVoxel.Interval getInterval(double center, double includePercent) {

			// count how many points are on each side of the center
			int numNeg = 0;
			int numPos = 0;
			for (double a : angles) {
				double delta = Protractor.getDeltaDegrees(center, a);
				if (delta < 0) {
					numNeg++;
				} else if (delta > 0) {
					numPos++;
				}
			}

			// figure out how many samples to skip for each side
			double skipRatio = 1.0 - includePercent/100.0;
			int numNegToSkip = (int)(numNeg*skipRatio);
			int numPosToSkip = (int)(numPos*skipRatio);

			SmallAngleVoxel.Interval interval = new SmallAngleVoxel.Interval(center);

			// do the negative side
			int count = 0;
			for (double a : angles) {
				if (count++ > numNegToSkip) {
					interval.less = Protractor.getDeltaDegrees(center, a);
					break;
				}
			}

			// do the positive side
			count = 0;
			for (double a : angles.descendingSet()) {
				if (count++ > numPosToSkip) {
					interval.more = Protractor.getDeltaDegrees(center, a);
					break;
				}
			}

			return interval;
		}
	}

	public SmallAngleVoxel.Interval bounds() {
		SmallAngleVoxel.Interval bounds = new SmallAngleVoxel.Interval(angles.first());
		for (double a : angles) {
			bounds.expand(a);
		}
		return bounds;
	}

	public double mean(SmallAngleVoxel.Interval bounds) {

		// use the interval to break the curvature of the space to we can compute a mean

		double sum = 0;
		for (double a : angles) {
			 sum += bounds.center + Protractor.getDeltaDegrees(bounds.center, a);
		}
		return Protractor.normalizeDegrees(sum/angles.size());
	}

	/** use a histogram with a certain number of buckets to find the mode */
	public double mode(SmallAngleVoxel.Interval bounds, int numBuckets) {

		double min = bounds.min();
		double width = bounds.size();

		Function<Double,Integer> angleToBucket = a -> {
			double deltaMin = Protractor.getDeltaDegrees(bounds.center, a) - bounds.less;
			return (int)(deltaMin*numBuckets/width);
		};
		Function<Integer,Double> bucketToAngle = i -> Protractor.normalizeDegrees(width*i/numBuckets + min);

		// build the histogram
		int[] counts = new int[numBuckets];
		Arrays.fill(counts, 0);
		for (double a : angles) {
			int i = angleToBucket.apply(a);
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

		return bucketToAngle.apply(bucket);
	}

	public double medoid() {

		// calculating means in curved spaces is hard, let's try medoids instead
		// this is a navie O(n2) algorithm, but it's fast enough for this amount of data
		// (barely!)

		double minScore = Double.POSITIVE_INFINITY;
		Double bestAngle = null;

		// copy angles into a primitive array for speed
		// iterating the tree is pretty slow!
		double[] fastAngles = new double[angles.size()];
		int i = 0;
		for (double a : angles) {
			fastAngles[i++] = a;
		}

		for (double a1 : fastAngles) {

			double score = 0.0;
			for (double a2 : fastAngles) {
				double d = Protractor.getDistDegrees(a1, a2);
				score += d*d;
			}

			if (score < minScore) {
				minScore = score;
				bestAngle = a1;
			}
		}

		assert (bestAngle != null);
		return bestAngle;
	}
}
