package edu.duke.cs.osprey.tools;

import java.util.ArrayList;
import java.util.List;

public abstract class SVGPlot {

	public static class LeftVerticalAxis {

		public double min = 0.0;
		public double max = 100.0;

		public double x = 0.0;
		public double y = 0.0;

		public List<Double> ticks = new ArrayList<>();

		public double tickLength = 4.0;
		public String tickFormat = "%.1f";
		public double tickTextMargin = 4.0;
		public double tickTextDy = -2;

		public SVG.StyleClass lineStyle = new SVG.StyleClass("plot-axis-line");
		public SVG.StyleClass tickTextStyle = new SVG.StyleClass("plot-axis-text");

		public LeftVerticalAxis() {

			// config default styles

			lineStyle.setStrokeColor(0x000000);
			lineStyle.setStrokeWidth(0.5);

			tickTextStyle.setFillColor(0x000000);
			tickTextStyle.setFontSize(5.0, SVG.LengthUnit.px);
			tickTextStyle.setNoStroke();
			tickTextStyle.setTextAnchor(SVG.StyleClass.TextAnchor.End);
		}

		public void addTicksOn(double multiple) {
			for (double val = Math.floor(min/multiple); val < max; val += multiple) {
				ticks.add(val);
			}
		}

		public void draw(SVG svg) {

			svg.putStyleClasses(lineStyle, tickTextStyle);

			// axis line
			svg.makeLine(
					x, y + min,
					x, y + max
				)
				.setStyleClasses(lineStyle)
				.draw();

			// draw the ticks
			drawTick(svg, min);
			for (double y : ticks) {
				drawTick(svg, y);
			}
			drawTick(svg, max);
		}

		private void drawTick(SVG svg, double y) {

			// draw the line
			svg.makeLine(
				x, y,
				x - tickLength, y
			)
				.setStyleClasses(lineStyle)
				.draw();

			// add the text
			svg.makeText(String.format(tickFormat, y))
				.setPos(x - tickLength - tickTextMargin, y)
				.setDY(tickTextDy, SVG.LengthUnit.px)
				.setStyleClasses(tickTextStyle)
				.draw();
		}
	}

	public static class Intervals {

		public static class Interval {

			public double min;
			public double max;
			public String id;
			public SVG.StyleClass extraStyle;

			public Interval(double min, double max) {
				this.min = min;
				this.max = max;
				this.id = null;
				this.extraStyle = null;
			}
		}

		public double intervalWidth = 8.0;
		public double intervalSpacing = 2.0;
		public double axisTicksOn = 10.0;
		public double minRectHeight = 0.2;

		public SVG.StyleClass intervalStyle = new SVG.StyleClass("plot-intervals-bar");

		public double xmin = 0.0;
		public double xmax = 0.0;
		public double ymin = 0.0;
		public double ymax = 0.0;

		public LeftVerticalAxis axis = null;

		private List<Interval> intervals = new ArrayList<>();

		public Intervals() {

			// config styles
			intervalStyle.setStrokeWidth(0.2);
			intervalStyle.setStrokeColor(0x666666);
			intervalStyle.setFillColor(0xcccccc);
		}

		public Interval addInterval(double min, double max) {

			// make the interval
			Interval interval = new Interval(min, max);
			intervals.add(interval);

			// update bounds
			xmax = intervals.size()*(intervalSpacing + intervalWidth);
			ymax = Math.max(ymax, max);

			return interval;
		}

		public LeftVerticalAxis makeAxis() {
			SVGPlot.LeftVerticalAxis axis = new SVGPlot.LeftVerticalAxis();
			axis.min = 0.0;
			axis.max = Math.ceil(ymax/axisTicksOn)*axisTicksOn;
			axis.addTicksOn(axisTicksOn);
			return axis;
		}

		public void draw(SVG svg) {

			svg.putStyleClasses(intervalStyle);

			int n = intervals.size();
			for (int i=0; i<n; i++) {
				Interval interval = intervals.get(i);

				double x = (i+1)* intervalSpacing + i* intervalWidth;
				double y1 = interval.min;
				double y2 = interval.max;

				if (Double.isInfinite(y1) || Double.isInfinite(y2)) {
					System.err.println(String.format("WARNING: interval [%f,%f] is infinite and will not be drawn", y1, y2));
					continue;
				}

				// if the interval is big enough, draw a rect
				SVG.Drawable d;
				if (y2 - y1 >= minRectHeight) {
					d = svg.makeRect(x, x + intervalWidth, y1, y2);
				} else {
					d = svg.makeLine(
						x, y1,
						x + intervalWidth, y1
					);
				}
				d.setStyleClasses(intervalStyle, interval.extraStyle)
					.setId(interval.id)
					.draw();
			}

			// make an axis if needed
			if (axis == null) {
				axis = makeAxis();
			}
			axis.draw(svg);
		}

		public void setBounds(SVG svg, double margin, double tickMargin) {
			svg.setBounds(
				xmin - margin - tickMargin, xmax + margin,
				ymin - margin, ymax + margin
			);
		}
	}

	public static class Boxes {

		public static class Box {

			public double xmin;
			public double xmax;
			public double ymin;
			public double ymax;
			public String id;
			public SVG.StyleClass extraStyle;

			public Box(double xmin, double xmax, double ymin, double ymax) {
				this.xmin = xmin;
				this.xmax = xmax;
				this.ymin = ymin;
				this.ymax = ymax;
				this.id = null;
				this.extraStyle = null;
			}

			public double getDX() {
				return xmax - xmin;
			}

			public double getDY() {
				return ymax - ymin;
			}
		}

		public double axisTicksOn = 10.0;
		public double minBoxSize = 0.2;

		public SVG.StyleClass boxStyle = new SVG.StyleClass("plot-boxes-box");

		public double xmin = 0.0;
		public double xmax = 0.0;
		public double ymin = 0.0;
		public double ymax = 0.0;

		public LeftVerticalAxis yaxis = null;
		// TODO
		// public BottomHorizontalAxis xaxis = null;

		private List<Box> boxes = new ArrayList<>();

		public Boxes() {

			// config styles
			boxStyle.setStrokeWidth(0.2);
			boxStyle.setStrokeColor(0x666666);
			boxStyle.setFillColor(0xcccccc);
		}

		public Box addBox(double xmin, double xmax, double ymin, double ymax) {

			// make the interval
			Box box = new Box(xmin, xmax, ymin, ymax);
			boxes.add(box);

			// update bounds
			this.xmax = Math.max(this.xmax, xmax);
			this.ymax = Math.max(this.ymax, ymax);

			return box;
		}

		public LeftVerticalAxis makeYAxis() {
			SVGPlot.LeftVerticalAxis axis = new SVGPlot.LeftVerticalAxis();
			axis.min = 0.0;
			axis.max = Math.ceil(ymax/axisTicksOn)*axisTicksOn;
			axis.addTicksOn(axisTicksOn);
			return axis;
		}

		/* TODO
		public LeftVerticalAxis makeXAxis() {
			SVGPlot.LeftVerticalAxis axis = new SVGPlot.LeftVerticalAxis();
			axis.min = 0.0;
			axis.max = Math.ceil(ymax/axisTicksOn)*axisTicksOn;
			axis.addTicksOn(axisTicksOn);
			return axis;
		}
		*/

		public void draw(SVG svg) {

			svg.putStyleClasses(boxStyle);

			for (Box box : boxes) {

				// just in case...
				if (Double.isInfinite(box.xmin) || Double.isInfinite(box.xmax) || Double.isInfinite(box.ymin) || Double.isInfinite(box.ymax)) {
					System.err.println(String.format("WARNING: box [%f,%f]x[%f,%f] is infinite and will not be drawn",
						box.xmin, box.xmax, box.ymin, box.ymax
					));
					continue;
				}

				// if the interval is big enough, draw a rect
				SVG.Drawable d;
				if (box.getDX() >= minBoxSize && box.getDY() >= minBoxSize) {
					// draw a real rect
					d = svg.makeRect(
						box.xmin, box.xmax,
						box.ymin, box.ymax
					);
				} else if (box.getDX() < minBoxSize && box.getDY() < minBoxSize) {
					// TODO: draw a point?
					throw new UnsupportedOperationException("draw a point?");
				} else if (box.getDX() < minBoxSize) {
					// draw a vertical line
					d = svg.makeLine(
						box.xmin, box.ymin,
						box.xmin, box.ymax
					);
				} else { // if (box.getDY() < minBoxSize) {
					// draw a horizontal line
					d = svg.makeLine(
						box.xmin, box.ymin,
						box.xmax, box.ymin
					);
				}
				d.setStyleClasses(boxStyle, box.extraStyle)
					.setId(box.id)
					.draw();
			}

			// make an axis if needed
			if (yaxis == null) {
				yaxis = makeYAxis();
			}
			yaxis.draw(svg);
		}

		public void setBounds(SVG svg, double margin, double tickMargin) {
			svg.setBounds(
				xmin - margin - tickMargin, xmax + margin,
				ymin - margin - tickMargin, ymax + margin
			);
		}
	}
}
