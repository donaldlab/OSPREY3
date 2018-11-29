package edu.duke.cs.osprey.tools.resultdoc;


import edu.duke.cs.osprey.tools.Streams;

import java.io.IOException;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;


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
		new Color("ffd320"), // yellow
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

	public abstract class Series {

		public String name = null;
		public Color color = null;

		public boolean hasXData = false;
		public boolean hasYData = false;

		public abstract String renderData();
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

	private void render(GNUPlot.Cmd c)
	throws IOException {

		// render each series as a datablock
		for (int i=0; i<series.size(); i++) {
			c.datablock("series" + i, series.get(i).renderData());
		}

		// x axis formatting
		c.xlabel.set(xlabel);
		boolean showX = series.stream().anyMatch(series -> series.hasXData)
			|| xlabels != null
			|| xlabelrotate != null;
		if (showX) {
			if (xlabels != null) {
				c.xtics.set(makeIndexedTics(xlabels));
			}
			c.xtics.rotate.set(xlabelrotate);
		} else {
			c.xtics.hide();
		}

		// y axis formatting
		c.ylabel.set(ylabel);
		boolean showY = series.stream().anyMatch(series -> series.hasYData)
			|| ylabels != null;
		if (showY) {
			if (ylabels != null) {
				c.ytics.set(makeIndexedTics(ylabels));
			}
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
			hasXData = true;
		}

		@Override
		public String renderData() {
			return Streams.joinToString(data, "\n", d -> d.toString());
		}

		@Override
		public String renderCols() {
			return "1:0";
		}
	}

	public class PointsY extends Points {

		public List<Double> data = null;

		public PointsY() {
			hasYData = true;
		}

		@Override
		public String renderData() {
			return Streams.joinToString(data, "\n", d -> d.toString());
		}

		@Override
		public String renderCols() {
			return "0:1";
		}
	}

	// TODO: PointsXY class?

	public abstract class Lines extends Series {

		@Override
		public String renderStyle(int seriesIndex, List<Color> colors) {
			return String.format("with lines linecolor %s", pickColor(seriesIndex, colors).render());
		}
	}

	public class LinesX extends Lines {

		public List<Double> data = null;

		public LinesX() {
			hasXData = true;
		}

		@Override
		public String renderData() {
			return Streams.joinToString(data, "\n", d -> d.toString());
		}

		@Override
		public String renderCols() {
			return "1:0";
		}
	}

	public class LinesY extends Lines {

		public List<Double> data = null;

		public LinesY() {
			hasYData = true;
		}

		@Override
		public String renderData() {
			return Streams.joinToString(data, "\n", d -> d.toString());
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
	}

	public abstract class Intervals extends Series {

		public List<Interval> data = null;

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
			hasXData = true;
		}

		public double barheight = 0.8;

		@Override
		public String renderData() {
			double barHalfHeight = barheight/2;
			AtomicInteger i = new AtomicInteger(0);
			return Streams.joinToString(data, "\n", interval -> {
				double x = (interval.max + interval.min)/2;
				double xlo = interval.min;
				double xhi = interval.max;
				double y = i.getAndIncrement();
				double ylo = y - barHalfHeight;
				double yhi = y + barHalfHeight;
				return String.format("%f %f %f %f %f %f", x, y, xlo, xhi, ylo, yhi);
			});
		}
	}

	public class IntervalsY extends Intervals {

		public IntervalsY() {
			hasYData = true;
		}

		public double barwidth = 0.8;

		@Override
		public String renderData() {
			double barHalfWidth = barwidth/2;
			AtomicInteger i = new AtomicInteger(0);
			return Streams.joinToString(data, "\n", interval -> {
				double x = i.getAndIncrement();
				double xlo = x - barHalfWidth;
				double xhi = x + barHalfWidth;
				double y = (interval.max + interval.min)/2;
				double ylo = interval.min;
				double yhi = interval.max;
				return String.format("%f %f %f %f %f %f", x, y, xlo, xhi, ylo, yhi);
			});
		}
	}
}
