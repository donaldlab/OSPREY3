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

import java.util.*;


public class DegreesHistogram {

	public final int numAngles;

	// use sparse storage, so we can handle high dimensions
	public final Map<Long,Integer> buckets = new TreeMap<>();

	public DegreesHistogram(int numAngles) {

		// need at least 9 bits per angle
		if (numAngles > 7) {
			throw new IllegalArgumentException("only up to 7 angles supported");
		}

		this.numAngles = numAngles;
	}

	private static int getNumAngles(Collection<double[]> angles) {

		if (angles.isEmpty()) {
			return 0;
		}

		return angles.iterator().next().length;
	}

	public DegreesHistogram(Collection<double[]> angles) {
		this(getNumAngles(angles));
		addAll(angles);
	}

	public void add(double[] angles) {
		// increment the right bucket
		buckets.merge(getIndex(angles), 1, (a, b) -> a + b);
	}

	public void addAll(Collection<double[]> angles) {
		for (double[] a : angles) {
			add(a);
		}
	}

	public long getIndex(double[] angles) {

		long i = 0;
		for (int d=0; d<numAngles; d++) {
			if (d > 0) {
				i <<= 9;
			}
			// put a bucket between every degree
			i |= (int)Math.floor(angles[d]) + 180;
		}
		return i;
	}

	public long getIndexDelta(int d) {
		long delta = 1;
		for (int i=numAngles-1; i>d; i--) {
			delta <<= 9;
		}
		return delta;
	}

	private int project(long index, int d) {
		return (int)(index >> 9*(numAngles - d - 1)) & 0x1ff;
	}

	/** return the sum of all the buckets */
	public int count() {
		int count = 0;
		for (int c : buckets.values()) {
			count += c;
		}
		return count;
	}

	public double[] makeDihedrals(long key) {
		double[] dihedrals = new double[numAngles];
		for (int d=numAngles-1; d>=0; d--) {
			long i = key & 0x1ff;
			key >>= 9;
			dihedrals[d] = (double)(i - 180);
		}
		return dihedrals;
	}

	public boolean hasBucket(double[] angles) {
		return buckets.containsKey(getIndex(angles));
	}

	/** get all buckets along angle d in [-180,179] order */
	public int[] getCounts(int d) {
		int[] out = new int[360];
		Arrays.fill(out, 0);
		for (Map.Entry<Long,Integer> bucket : buckets.entrySet()) {
			out[project(bucket.getKey(), d)] += bucket.getValue();
		}
		return out;
	}

	public static final String[] AngleLadder = {
		"-1        -1        -1        -1        -1        -1        -1        -1        -1                                                                                                                                                                                                       1         1         1         1         1         1         1         1        1 ",
		" 8         7         6         5         4         3         2         1         0        -9        -8        -7        -6        -5        -4        -3        -2        -1                   1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7        7 ",
		"[098765432109876543210987654321098765432109876543210987654321098765432109876543210987654321098765432109876543210987654321098765432109876543210987654321098765432109876543210987654321012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789]"
	};

	public String dump() {
		StringBuilder buf = new StringBuilder();

		buf.append("angles: ");
		buf.append(count());

		buf.append("   buckets filled: ");
		buf.append(buckets.size());

		buf.append("\n    ");
		buf.append(AngleLadder[0]);
		buf.append("\n    ");
		buf.append(AngleLadder[1]);
		buf.append("\n    ");
		buf.append(AngleLadder[2]);

		for (int d=0; d<numAngles; d++) {
			buf.append("\n ");
			buf.append(d);
			buf.append("  ");
			buf.append(dump(d));
		}

		return buf.toString();
	}

	public String dump(int d) {

		int[] counts = getCounts(d);

		// find the max
		int max = 0;
		for (int count : counts) {
			max = Math.max(max, count);
		}

		// render to a string
		StringBuilder buf = new StringBuilder();
		buf.append('[');
		for (int count : counts) {

			// no data at all?
			if (max <= 0) {
				buf.append(" ");
				continue;
			}

			// normalize the counts
			double n = 1.0*count/max;

			// quantize with the symbols " _.o8|"
			if (n <= 0.0) {
				buf.append(' ');
			} else if (n < 0.1) {
				buf.append('_');
			} else if (n < 0.5) {
				buf.append('.');
			} else if (n < 0.75) {
				buf.append('o');
			} else if (n < 1.0) {
				buf.append('8');
			} else if (n >= 1.0) {
				buf.append('|');
			} else if (Double.isNaN(n)) {
				buf.append('N');
			} else {
				buf.append('?');
			}
		}
		buf.append("]");
		return buf.toString();
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
