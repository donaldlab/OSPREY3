package edu.duke.cs.osprey.coffee;


import java.math.BigDecimal;
import java.util.Arrays;

/**
 * Essentially a histogram of gaps in the energy bounds.
 */
public class EnergyBoundStats {

	private final double[] boundaries = {
		-10.0,
		-5.0,
		-2.0,
		-1.0,
		-0.1,
		0.0,
		0.1,
		1.0,
		2.0,
		5.0,
		10.0,
		20.0,
		50.0,
		100.0
	};
	private final long[] buckets = new long[boundaries.length + 1];
	private long count = 0;
	private BigDecimal sum = BigDecimal.ZERO;

	public void add(double lowerBound, double energy) {

		// compute the signed energy gap
		double gap = energy - lowerBound;

		// find the bucket index
		int index = 0;
		for (; index<buckets.length; index++) {
			if (gap < boundaries[index]) {
				break;
			}
		}

		buckets[index]++;
		count++;
		sum = sum.add(BigDecimal.valueOf(gap));
	}

	@Override
	public String toString() {

		if (count <= 0) {
			return "Energy bound gap statistics: none";
		}

		// compute some summary statistics
		double[] percents = Arrays.stream(buckets)
			.mapToDouble(bucket -> 100.0*bucket/count)
			.toArray();
		double mean = sum.doubleValue()/count;

		var buf = new StringBuilder();
		buf.append(String.format("Energy bound gap statistics:   count %d   mean %.4f\n", count, mean));
		for (int i=0; i<buckets.length; i++) {
			if (i < boundaries.length) {
				buf.append(String.format(" < %5.1f: %9d  %5.1f%%", boundaries[i], buckets[i], percents[i]));
			} else {
				buf.append(String.format(">= %5.1f: %9d  %5.1f%%", boundaries[boundaries.length - 1], buckets[i], percents[i]));
			}
			buf.append("  |");
			buf.append("=".repeat((int)(percents[i]*20/100.0)));
			if (i < boundaries.length) {
				buf.append('\n');
			}
		}
		return buf.toString();
	}
}
