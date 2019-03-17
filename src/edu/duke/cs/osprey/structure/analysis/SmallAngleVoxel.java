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

import java.util.Collection;
import java.util.List;
import java.util.stream.Collectors;


/** a cartesian product of intervals in S^1 x S^1 x S^1 ... localized around a point */
public class SmallAngleVoxel {

	/** an interval on S^1, localized around a point */
	public static class Interval {
		
		public double center;
		public double less;
		public double more;

		public Interval(double center) {
			this(0, center, 0);
		}

		public Interval(double less, double center, double more) {
			this.center = center;
			this.less = less;
			this.more = more;
		}

		public void expand(double p) {
			double delta = Protractor.getDeltaDegrees(center, p);
			if (delta < less) {
				less = delta;
			}
			if (delta > more) {
				more = delta;
			}
		}

		public double size() {
			return more - less;
		}

		public boolean contains(double p) {
			double delta = Protractor.getDeltaDegrees(center, p);
			return delta >= less && delta <= more;
		}

		public double min() {
			return Protractor.normalizeDegrees(center + less);
		}

		public double max() {
			return Protractor.normalizeDegrees(center + more);
		}

		@Override
		public String toString() {
			return String.format("[%6.1f, %6.1f, %6.1f] -> [%6.1f,%6.1f]  width=%6.1f", less, center, more, min(), max(), more - less);
		}
	}

	public Interval[] intervals;


	public SmallAngleVoxel(Interval ... intervals) {
		this.intervals = intervals;
	}

	public SmallAngleVoxel(double[] p) {
		intervals = new Interval[p.length];
		for (int d=0; d<p.length; d++) {
			intervals[d] = new Interval(p[d]);
		}
	}

	public static SmallAngleVoxel makeFromBounds(double[] bounds) {

		if (bounds.length % 3 != 0) {
			throw new IllegalArgumentException("bounds could come in multiples of 3");
		}
		int n = bounds.length/3;

		SmallAngleVoxel.Interval[] voxels = new SmallAngleVoxel.Interval[n];
		for (int i=0; i<n; i++) {
			int i3 = i*3;
			voxels[i] = new SmallAngleVoxel.Interval(bounds[i3], bounds[i3 + 1], bounds[i3 + 2]);
		}

		return new SmallAngleVoxel(voxels);
	}

	public void expand(double[] p) {
		for (int d=0; d<p.length; d++) {
			intervals[d].expand(p[d]);
		}
	}

	public double volume() {
		if (intervals.length == 0) {
			return 0.0;
		}
		double v = 1.0;
		for (Interval interval : intervals) {
			v *= interval.size();
		}
		return v;
	}

	public boolean contains(double[] p) {
		for (int d=0; d<p.length; d++) {
			if (!intervals[d].contains(p[d])) {
				return false;
			}
		}
		return true;
	}

	@Override
	public String toString() {
		StringBuilder buf = new StringBuilder();
		for (int d=0; d<intervals.length; d++) {
			if (d > 0) {
				buf.append("\n");
			}
			buf.append(intervals[d].toString());
		}
		return buf.toString();
	}

	public List<double[]> filter(Collection<double[]> points) {
		return points.stream()
			.filter(p -> contains(p))
			.collect(Collectors.toList());
	}
}
