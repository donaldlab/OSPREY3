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

package edu.duke.cs.osprey.tools.resultdoc;


import edu.duke.cs.osprey.tools.Streams;

import java.io.IOException;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.function.Supplier;


public class Plot {

	public static class Color {

		public final String  hex;

		public Color(String hex) {
			this.hex = hex;
		}

		public String render() {
			return "rgbcolor \"#" + hex + "\"";
		}
	}

	public static List<Color> defaultColors = Arrays.asList(
		new Color("004586"), // blue
		new Color("ff420e"), // red
		//new Color("ffd320"), // yellow (doesn't show up well against the white background)
		new Color("579d1c"), // green
		new Color("7e0021"), // magenta
		new Color("83caff"), // cyan
		new Color("314004"), // dark green
		new Color("aecf00"), // lime
		new Color("4b1f6f")  // purple
	);

	/** get the i-th item from the list, using modular arithmetic */
	public static <T> T modGet(List<T> items, int i) {
		return items.get(i % items.size());
	}

	public static enum PointType {

		Dot(0),
		Plus(1),
		X(2),
		Star(3),
		Square(4),
		FilledSquare(5),
		Circle(6),
		CircleFilled(7),
		TriangleUp(8),
		TriangleUpFilled(9),
		TriangleDown(10),
		TriangleDownFilled(11),
		Diamond(12),
		DiamondFilled(13);

		public final int id;

		PointType(int id) {
			this.id = id;
		}
	}

	public final List<Series> series = new ArrayList<>();

	private enum Axis {

		X,
		Y;

		public static <T> Map<Axis,T> makeMap() {
			return new EnumMap<>(Axis.class);
		}
	}

	private static abstract class AxisInfo {

		double minFinite = Double.NaN;
		double maxFinite = Double.NaN;
		boolean hasNeginf = false;
		boolean hasPosinf = false;

		/** update min,max and has inf flags */
		abstract void update();

		String neginf = null;
		String posinf = null;


		static class DoubleList extends AxisInfo {

			private Supplier<List<Double>> dataSupplier;

			DoubleList(Supplier<List<Double>> dataSupplier) {
				this.dataSupplier = dataSupplier;
			}

			@Override
			void update() {
				List<Double> data = dataSupplier.get();
				minFinite = data.stream()
					.filter(Objects::nonNull)
					.filter(Double::isFinite)
					.min(Double::compare)
					.orElse(0.0);
				maxFinite = data.stream()
					.filter(Objects::nonNull)
					.filter(Double::isFinite)
					.max(Double::compare)
					.orElse(0.0);
				hasNeginf = data.stream().anyMatch(d -> d == Double.NEGATIVE_INFINITY);
				hasPosinf = data.stream().anyMatch(d -> d == Double.POSITIVE_INFINITY);
			}
		}
	}

	public abstract class Series {

		public String name = null;
		public Color color = null;

		protected Map<Axis,AxisInfo> axisInfos = Axis.makeMap();

		public abstract String renderData(Map<Axis,Double> neginfVals, Map<Axis,Double> posinfVals);
		public abstract String renderCols();
		public abstract String renderStyle(int seriesIndex, List<Color> colors);

		public Series() {
			series.add(this);
		}

		public Color pickColor(int seriesIndex, List<Color> colors) {

			// use the explicitly set color if possible
			if (color != null) {
				return color;
			}

			// otherwise, pick from the default colors
			return modGet(colors, seriesIndex);
		}
	}

	public int width = 600;
	public int height = 400;
	public String xlabel = null;
	public String ylabel = null;
	public List<Color> colors = defaultColors;
	public List<String> xlabels = null;
	public Double xlabelrotate = null;
	public List<String> ylabels = null;
	public String key = null;

	public byte[] renderPng() {
		return GNUPlot.plotPng(width, height, c -> render(c));
	}

	public String renderSvg() {
		return GNUPlot.plotSvg(width, height, c -> render(c));
	}

	private Map<Number,String> makeIndexedTics(List<String> labels) {
		Map<Number,String> out = new TreeMap<>();
		for (int i=0; i<labels.size(); i++) {
			out.put(i, labels.get(i));
		}
		return out;
	}

	private void addInfLabels(Axis axis, GNUPlot.Cmd c, GNUPlot.Cmd.Tics tics, Map<Axis,Double> neginfVals, Map<Axis,Double> posinfVals)
	throws IOException {

		boolean hasNeginf = series.stream().anyMatch(s -> {
			AxisInfo axisInfo = s.axisInfos.get(axis);
			return axisInfo != null && axisInfo.hasNeginf;
		});
		if (hasNeginf) {

			// add the tic
			double neginf = neginfVals.get(axis);
			tics.add(neginf, "-Inf");

			// update the axis range so the inf tick shows up on the end
			// but pad just a little so we don't screw up the other axis auto ranges
			c.command("set %srange [%f:]", axis.name().toLowerCase(), neginf - Math.abs(neginf)*0.01);
		}

		boolean hasPosinf = series.stream().anyMatch(s -> {
			AxisInfo axisInfo = s.axisInfos.get(axis);
			return axisInfo != null && axisInfo.hasPosinf;
		});
		if (hasPosinf) {

			// add the tic
			double posinf = posinfVals.get(axis);
			tics.add(posinf, "+Inf");

			// update the axis range so the inf tick shows up on the end
			// but pad just a little so we don't screw up the other axis auto ranges
			c.command("set %srange [:%f]", axis.name().toLowerCase(), posinf + Math.abs(posinf)*0.01);
		}
	}

	private static String renderDouble(Double val, Double neginf, Double posinf) {
		if (val == null) {
			return "";
		} else if (val == Double.NEGATIVE_INFINITY) {
			return neginf.toString();
		} else if (val == Double.POSITIVE_INFINITY) {
			return posinf.toString();
		} else {
			return val.toString();
		}
	}

	private void render(GNUPlot.Cmd c)
	throws IOException {

		// update each axis info
		for (Series s : series) {
			for (AxisInfo axisInfo : s.axisInfos.values()) {
				axisInfo.update();
				assert (Double.isFinite(axisInfo.minFinite));
				assert (Double.isFinite(axisInfo.maxFinite));
			}
		}

		// map inf to finite values outside the min,max range
		Map<Axis,Double> neginfVals = Axis.makeMap();
		Map<Axis,Double> posinfVals = Axis.makeMap();
		for (Axis axis : Axis.values()) {

			// get min,max for this axis over all series
			Double minFinite = series.stream()
				.map(s -> s.axisInfos.get(axis))
				.filter(Objects::nonNull)
				.map(axisInfo -> axisInfo.minFinite)
				.min(Double::compare)
				.orElse(null);
			if (minFinite == null) {
				continue;
			}
			Double maxFinite = series.stream()
				.map(s -> s.axisInfos.get(axis))
				.filter(Objects::nonNull)
				.map(axisInfo -> axisInfo.maxFinite)
				.max(Double::compare)
				.orElse(null);
			if (maxFinite == null) {
				continue;
			}

			double padFactor = 0.1;
			double pad = (maxFinite - minFinite)*padFactor;
			neginfVals.put(axis, minFinite - pad);
			posinfVals.put(axis, maxFinite + pad);
		}

		// render each series as a datablock
		for (int i=0; i<series.size(); i++) {
			c.datablock("series" + i, series.get(i).renderData(neginfVals, posinfVals));
		}

		// global formatting
		c.key.set(key);

		// make the axes look prettier by default
		c.command("set border 3 back");
		c.command("set tics nomirror");
		c.command("set style line 102 linecolor rgb \"#cccccc\" linetype 0 linewidth 1");
		c.command("set grid back linestyle 102");

		// x axis formatting
		c.xlabel.set(xlabel);
		boolean showX = series.stream().anyMatch(s -> s.axisInfos.containsKey(Axis.X))
			|| xlabels != null
			|| xlabelrotate != null;
		if (showX) {
			if (xlabels != null) {
				c.xtics.set(makeIndexedTics(xlabels));
			}
			addInfLabels(Axis.X, c, c.xtics, neginfVals, posinfVals);
			// negate the rotation here so positive rotations don't clash with the plot
			c.xtics.rotate.set(-xlabelrotate);
		} else {
			c.xtics.hide();
		}

		// y axis formatting
		c.ylabel.set(ylabel);
		boolean showY = series.stream().anyMatch(s -> s.axisInfos.containsKey(Axis.Y))
			|| ylabels != null;
		if (showY) {
			if (ylabels != null) {
				c.ytics.set(makeIndexedTics(ylabels));
			}
			addInfLabels(Axis.Y, c, c.ytics, neginfVals, posinfVals);
		} else {
			c.ytics.hide();
		}

		// render all series
		List<String> seriesCmds = new ArrayList<>();
		for (int i=0; i<series.size(); i++) {
			Series s = series.get(i);
			seriesCmds.add(String.format("$series%d using %s title \"%s\" %s",
				i,
				s.renderCols(),
				s.name != null ? s.name : "Series " + i,
				s.renderStyle(i, colors)
			));
		}
		c.command("plot " + String.join(", ", seriesCmds));
	}


	public abstract class Points extends Series {

		public PointType type = null;
		public Double size = null;

		protected String renderData1D(List<Double> data, Double neginf, Double posinf) {
			return Streams.joinToString(data, "\n", d -> renderDouble(d, neginf, posinf));
		}

		@Override
		public String renderStyle(int seriesIndex, List<Color> colors) {
			return String.format("with points linecolor %s %s %s",
				pickColor(seriesIndex, colors).render(),
				type == null ? "" : "pointtype " + type.id,
				size == null ? "" : "pointsize " + size
			);
		}
	}

	public class PointsX extends Points {

		public List<Double> data = null;

		public PointsX() {
			axisInfos.put(Axis.X, new AxisInfo.DoubleList(() -> data));
		}

		@Override
		public String renderData(Map<Axis,Double> neginfVals, Map<Axis,Double> posinfVals) {
			return renderData1D(data, neginfVals.get(Axis.X), posinfVals.get(Axis.X));
		}

		@Override
		public String renderCols() {
			return "1:0";
		}
	}

	public class PointsY extends Points {

		public List<Double> data = null;

		public PointsY() {
			axisInfos.put(Axis.Y, new AxisInfo.DoubleList(() -> data));
		}

		@Override
		public String renderData(Map<Axis,Double> neginfVals, Map<Axis,Double> posinfVals) {
			return renderData1D(data, neginfVals.get(Axis.Y), posinfVals.get(Axis.Y));
		}

		@Override
		public String renderCols() {
			return "0:1";
		}
	}

	// TODO: PointsXY class?

	public abstract class Lines extends Series {

		protected String renderData1D(List<Double> data, Double neginf, Double posinf) {
			return Streams.joinToString(data, "\n", d -> renderDouble(d, neginf, posinf));
		}

		@Override
		public String renderStyle(int seriesIndex, List<Color> colors) {
			return String.format("with lines linecolor %s", pickColor(seriesIndex, colors).render());
		}
	}

	public class LinesX extends Lines {

		public List<Double> data = null;

		public LinesX() {
			axisInfos.put(Axis.X, new AxisInfo.DoubleList(() -> data));
		}

		@Override
		public String renderData(Map<Axis,Double> neginfVals, Map<Axis,Double> posinfVals) {
			return renderData1D(data, neginfVals.get(Axis.X), posinfVals.get(Axis.X));
		}

		@Override
		public String renderCols() {
			return "1:0";
		}
	}

	public class LinesY extends Lines {

		public List<Double> data = null;

		public LinesY() {
			axisInfos.put(Axis.Y, new AxisInfo.DoubleList(() -> data));
		}

		@Override
		public String renderData(Map<Axis,Double> neginfVals, Map<Axis,Double> posinfVals) {
			return renderData1D(data, neginfVals.get(Axis.Y), posinfVals.get(Axis.Y));
		}

		@Override
		public String renderCols() {
			return "0:1";
		}
	}

	// TODO: LinesXY?

	public static class Interval {

		public double min;
		public double max;

		public Interval(double min, double max) {
			this.min = min;
			this.max = max;
		}

		/** return the smallest finite endpoint, or null */
		public Double finiteMin() {
			if (Double.isFinite(min)) {
				return min;
			} else if (Double.isFinite(max)) {
				return max;
			} else {
				return null;
			}
		}

		/** return the largest finite endpoint, or null */
		public Double finiteMax() {
			if (Double.isFinite(max)) {
				return max;
			} else if (Double.isFinite(min)) {
				return min;
			} else {
				return null;
			}
		}
	}

	public abstract class Intervals extends Series {

		public List<Interval> data = null;

		protected Intervals(Axis axis) {
			axisInfos.put(axis, new AxisInfo() {

				@Override
				void update() {
					minFinite = data.stream()
						.filter(Objects::nonNull)
						.map(Interval::finiteMin)
						.filter(Objects::nonNull)
						.min(Double::compare)
						.orElse(0.0);
					maxFinite = data.stream()
						.filter(Objects::nonNull)
						.map(Interval::finiteMax)
						.filter(Objects::nonNull)
						.max(Double::compare)
						.orElse(0.0);
					hasNeginf = data.stream()
						.filter(Objects::nonNull)
						.anyMatch(i -> i.min == Double.NEGATIVE_INFINITY || i.max == Double.NEGATIVE_INFINITY);
					hasPosinf = data.stream()
						.filter(Objects::nonNull)
						.anyMatch(i -> i.min == Double.POSITIVE_INFINITY || i.max == Double.POSITIVE_INFINITY);
				}
			});
		}

		protected double getIntervalMid(Interval interval) {
			if (Double.isFinite(interval.min) && Double.isFinite(interval.max)) {
				return (interval.min + interval.max)/2;
			} else if (Double.isFinite(interval.min)) {
				return interval.min;
			} else if (Double.isFinite(interval.max)) {
				return interval.max;
			} else {
				return 0.0;
			}
		}

		@Override
		public String renderStyle(int seriesIndex, List<Color> colors) {
			return String.format("with boxxy linecolor %s", pickColor(seriesIndex, colors).render());
		}

		@Override
		public String renderCols() {
			return "1:2:3:4:5:6";
		}
	}

	public class IntervalsX extends Intervals {

		public IntervalsX() {
			super(Axis.X);
		}

		public double barheight = 0.8;

		@Override
		public String renderData(Map<Axis,Double> neginfVals, Map<Axis,Double> posinfVals) {
			Double neginf = neginfVals.get(Axis.X);
			Double posinf = posinfVals.get(Axis.X);
			double barHalfHeight = barheight/2;
			AtomicInteger i = new AtomicInteger(0);
			return Streams.joinToString(data, "\n", interval -> {
				if (interval == null) {
					i.getAndIncrement();
					return "";
				} else {
					double x = getIntervalMid(interval);
					double xlo = interval.min;
					double xhi = interval.max;
					double y = i.getAndIncrement();
					double ylo = y - barHalfHeight;
					double yhi = y + barHalfHeight;
					return String.format("%s %f %s %s %f %f",
						renderDouble(x, neginf, posinf),
						y,
						renderDouble(xlo, neginf, posinf),
						renderDouble(xhi, neginf, posinf),
						ylo, yhi
					);
				}
			});
		}
	}

	public class IntervalsY extends Intervals {

		public IntervalsY() {
			super(Axis.Y);
		}

		public double barwidth = 0.8;

		@Override
		public String renderData(Map<Axis,Double> neginfVals, Map<Axis,Double> posinfVals) {
			Double neginf = neginfVals.get(Axis.X);
			Double posinf = posinfVals.get(Axis.X);
			double barHalfWidth = barwidth/2;
			AtomicInteger i = new AtomicInteger(0);
			return Streams.joinToString(data, "\n", interval -> {
				if (interval == null) {
					i.getAndIncrement();
					return "";
				} else {
					double x = i.getAndIncrement();
					double xlo = x - barHalfWidth;
					double xhi = x + barHalfWidth;
					double y = getIntervalMid(interval);
					double ylo = interval.min;
					double yhi = interval.max;
					return String.format("%f %s %f %f %s %s",
						x,
						renderDouble(y, neginf, posinf),
						xlo, xhi,
						renderDouble(ylo, neginf, posinf),
						renderDouble(yhi, neginf, posinf)
					);
				}
			});
		}
	}
}
