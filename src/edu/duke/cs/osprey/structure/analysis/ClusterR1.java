package edu.duke.cs.osprey.structure.analysis;

import edu.duke.cs.osprey.tools.MathTools;

import java.util.Arrays;
import java.util.TreeSet;
import java.util.function.Function;

import static edu.duke.cs.osprey.tools.Log.log;


/** a cluster of data points in R^1 */
public class ClusterR1 {

	private final TreeSet<Double> values = new TreeSet<>();

	public void add(double a) {
		values.add(a);
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
