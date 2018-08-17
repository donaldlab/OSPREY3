package edu.duke.cs.osprey.structure.analysis;

import edu.duke.cs.osprey.tools.MathTools;

import java.util.*;


public class DegreesHistogram {

	public final int numAngles;

	// use sparse storage, so we can handle high dimensions
	public final Map<Long,Integer> buckets;

	public DegreesHistogram(int numAngles) {

		this.numAngles = numAngles;

		buckets = new TreeMap<>();
	}

	public long getIndex(double[] angles) {

		long i = 0;
		for (int d=0; d<numAngles; d++) {
			if (d > 0) {
				i <<= 16;
			}
			// put a bucket between every degree
			i |= (int)Math.floor(angles[d]) + 180;
		}
		return i;
	}

	public double[] makeDihedrals(long key) {
		double[] dihedrals = new double[numAngles];
		for (int d=numAngles-1; d>=0; d--) {
			long i = key & 0xffff;
			key >>= 16;
			dihedrals[d] = (double)(i - 180);
		}
		return dihedrals;
	}

	public boolean hasBucket(double[] angles) {
		return buckets.containsKey(getIndex(angles));
	}

	public void add(double[] angles) {
		// increment the right bucket
		buckets.merge(getIndex(angles), 1, (a, b) -> a + b);
	}

	public long getIndexDelta(int d) {
		long delta = 1;
		for (int i=numAngles-1; i>d; i--) {
			delta <<= 16;
		}
		return delta;
	}

	/** return the sum of all the buckets */
	public int count() {
		int count = 0;
		for (int c : buckets.values()) {
			count += c;
		}
		return count;
	}

	/** remove buckets whose local neighborhood contains fewer than the threshold count */
	public void filterDensityWindow(int radiusDegrees, int thresholdCount) {

		final int widthDegrees = radiusDegrees*2 + 1;

		final int[] widths = new int[numAngles];
		Arrays.fill(widths, widthDegrees);

		List<Long> bucketsToRemove = new ArrayList<>();
		for (Map.Entry<Long,Integer> bucket : buckets.entrySet()) {

			// sum over the window
			int windowCount = 0;
			for (int[] indices : new MathTools.GridIterable(widths)) {
				long index = bucket.getKey();
				for (int d=0; d<numAngles; d++) {
					index += getIndexDelta(d)*(indices[d] - radiusDegrees);
				}
				windowCount += buckets.getOrDefault(index, 0);
			}

			// what ratio of points are here?
			if (windowCount < thresholdCount) {
				bucketsToRemove.add(bucket.getKey());
			}
		}

		// remove the buckets
		for (long key : bucketsToRemove) {
			buckets.remove(key);
		}
	}
}
