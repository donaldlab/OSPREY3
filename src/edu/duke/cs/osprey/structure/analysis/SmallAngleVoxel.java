package edu.duke.cs.osprey.structure.analysis;

import edu.duke.cs.osprey.tools.Protractor;

import java.util.Collection;
import java.util.List;
import java.util.stream.Collectors;


/** a cartesian product of intervals in S^1 x S^1 x S^1 ... localized around a point */
public class SmallAngleVoxel {

	/** an interval on S^1, localized around a point */
	public class Interval {
		
		public double center;
		public double less;
		public double more;

		public Interval(double center) {
			this.center = center;
			this.less = 0.0;
			this.more = 0.0;
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
			return String.format("[%6.1f,%6.1f,%6.1f] -> [%6.1f,%6.1f]", less, center, more, min(), max());
		}
	}

	public Interval[] intervals;
	

	public SmallAngleVoxel(double[] p) {
		intervals = new Interval[p.length];
		for (int d=0; d<p.length; d++) {
			intervals[d] = new Interval(p[d]);
		}
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
