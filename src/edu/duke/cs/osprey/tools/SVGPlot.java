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

package edu.duke.cs.osprey.tools;

import java.util.ArrayList;
import java.util.List;

public abstract class SVGPlot {

	public static abstract class Axis {

		public double x = 0.0;
		public double y = 0.0;

		public double min = 0.0;
		public double max = 100.0;

		public List<Double> ticks = new ArrayList<>();

		public double tickLength = 4.0;
		public String tickFormat = "%.1f";
		public double tickTextMargin = 2.0;

		public SVG.StyleClass lineStyle = new SVG.StyleClass("plot-axis-line");
		public SVG.StyleClass tickTextStyle = new SVG.StyleClass("plot-axis-text");

		public Axis() {

			// config default styles

			lineStyle.setStrokeColor(0x000000);
			lineStyle.setStrokeWidth(0.5);

			tickTextStyle.setFillColor(0x000000);
			tickTextStyle.setFontSize(5.0, SVG.LengthUnit.px);
			tickTextStyle.setNoStroke();
		}

		public void addTicksOn(double multiple) {
			for (double val = Math.floor(min/multiple)*multiple + multiple; val < max; val += multiple) {
				ticks.add(val);
			}
		}

		protected void preDraw(SVG svg) {
			svg.putStyleClasses(lineStyle, tickTextStyle);
		}
	}

	public static class LeftAxis extends Axis {

		public double tickTextDy = -2;

		public LeftAxis() {

			// set text alignment
			tickTextStyle.setTextAnchor(SVG.StyleClass.TextAnchor.End);
		}

		public void draw(SVG svg) {
			preDraw(svg);

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

	public static class BottomAxis extends Axis {

		public int tickTextHeight = 5;

		public BottomAxis() {

			// set text alignment
			tickTextStyle.setTextAnchor(SVG.StyleClass.TextAnchor.Middle);
		}

		public void draw(SVG svg) {
			super.preDraw(svg);

			// axis line
			svg.makeLine(
					x + min, y,
					x + max, y
				)
				.setStyleClasses(lineStyle)
				.draw();

			// draw the ticks
			drawTick(svg, min);
			for (double x : ticks) {
				drawTick(svg, x);
			}
			drawTick(svg, max);
		}

		private void drawTick(SVG svg, double x) {

			// draw the line
			svg.makeLine(
					x, y,
					x, y - tickLength
				)
				.setStyleClasses(lineStyle)
				.draw();

			// add the text
			svg.makeText(String.format(tickFormat, x))
				.setPos(x, y - tickLength - tickTextMargin - tickTextHeight)
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

		public LeftAxis axis = null;

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

		public LeftAxis makeAxis() {
			SVGPlot.LeftAxis axis = new SVGPlot.LeftAxis();
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

				svg.makeRect(x, x + intervalWidth, y1, y2)
					.setStyleClasses(intervalStyle, interval.extraStyle)
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

		public SVG.StyleClass boxStyle = new SVG.StyleClass("plot-boxes-box");

		public double xmin = 0.0;
		public double xmax = 0.0;
		public double ymin = 0.0;
		public double ymax = 0.0;

		public BottomAxis xaxis = null;
		public LeftAxis yaxis = null;

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

		public BottomAxis makeXAxis() {
			SVGPlot.BottomAxis axis = new SVGPlot.BottomAxis();
			axis.min = 0.0;
			axis.max = Math.ceil(xmax/axisTicksOn)*axisTicksOn;
			axis.addTicksOn(axisTicksOn);
			return axis;
		}

		public LeftAxis makeYAxis() {
			SVGPlot.LeftAxis axis = new SVGPlot.LeftAxis();
			axis.min = 0.0;
			axis.max = Math.ceil(ymax/axisTicksOn)*axisTicksOn;
			axis.addTicksOn(axisTicksOn);
			return axis;
		}

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
				if (box.getDX() == 0 && box.getDY() == 0) {
					// draw a point
					d = svg.makePoint(
						box.xmin, box.ymin,
						0.01
					);
				} else {
					// draw a real rect
					d = svg.makeRect(
						box.xmin, box.xmax,
						box.ymin, box.ymax
					);
				}
				d.setStyleClasses(boxStyle, box.extraStyle)
					.setId(box.id)
					.draw();
			}

			// make the axes if needed
			if (xaxis == null) {
				xaxis = makeXAxis();
			}
			xaxis.draw(svg);
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
