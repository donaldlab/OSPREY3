package edu.duke.cs.osprey.structure.analysis;

import edu.duke.cs.osprey.tools.Protractor;

import java.util.TreeSet;
import java.util.function.BiConsumer;


public class SmallAngleCluster {

	private final TreeSet<Double> angles = new TreeSet<>();

	public void add(double a) {
		angles.add(a);
	}


	public class Stats {

		public final SmallAngleCluster cluster = SmallAngleCluster.this;

		public final SmallAngleVoxel.Interval bounds = bounds();
		public final double medoid = medoid();

		public final int numNeg;
		public final int numPos;

		public Stats() {

			int neg = 0;
			int pos = 0;
			for (double a : angles) {
				double delta = Protractor.getDeltaDegrees(medoid, a);
				if (delta < 0) {
					neg++;
				} else if (delta > 0) {
					pos++;
				}
			}

			numNeg = neg;
			numPos = pos;
		}

		public SmallAngleVoxel.Interval getInterval(double includePercent) {

			double skipRatio = 1.0 - includePercent/100.0;
			double numNegToSkip = numNeg*skipRatio;
			double numPosToSkip = numPos*skipRatio;

			SmallAngleVoxel.Interval interval = new SmallAngleVoxel.Interval(medoid);

			// do the negative side
			int count = 0;
			for (double a : angles) {
				if (count++ > numNegToSkip) {
					interval.less = Protractor.getDeltaDegrees(medoid, a);
					break;
				}
			}

			// do the positive side
			count = 0;
			for (double a : angles.descendingSet()) {
				if (count++ > numPosToSkip) {
					interval.more = Protractor.getDeltaDegrees(medoid, a);
					break;
				}
			}

			return interval;
		}

		public String toString(double ... includePercents) {
			StringBuilder buf = new StringBuilder();

			buf.append("Small angle cluster:");

			BiConsumer<String,String> addLine = (label, value) ->
				buf.append(String.format("\n%20s: %s", label, value));

			addLine.accept("angles", String.format("%d", angles.size()));
			addLine.accept("bounds", bounds.toString());
			for (double includePercent : includePercents) {
				addLine.accept(
					String.format("%4.1f%% interval", includePercent),
					getInterval(includePercent).toString()
				);
			}

			return buf.toString();
		}

		@Override
		public String toString() {
			return toString(99, 95, 90, 75, 50);
		}
	}

	private SmallAngleVoxel.Interval bounds() {
		SmallAngleVoxel.Interval bounds = new SmallAngleVoxel.Interval(angles.first());
		for (double a : angles) {
			bounds.expand(a);
		}
		return bounds;
	}

	private double medoid() {

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
