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
				.setStyleClass(lineStyle)
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
				.setStyleClass(lineStyle)
				.draw();

			// add the text
			svg.makeText(String.format(tickFormat, y))
				.setPos(x - tickLength - tickTextMargin, y)
				.setDY(tickTextDy, SVG.LengthUnit.px)
				.setStyleClass(tickTextStyle)
				.draw();
		}
	}

	public static class Intervals {

		public double barWidth = 8.0;
		public double barSpacing = 2.0;
		public double axisTicksOn = 10.0;
		public double minRectHeight = 0.2;

		public SVG.StyleClass barStyle = new SVG.StyleClass("plot-intervals-bar");

		public double xmin = 0.0;
		public double xmax = 0.0;
		public double ymin = 0.0;
		public double ymax = 0.0;

		public LeftVerticalAxis axis = null;

		private List<Double> mins = new ArrayList<>();
		private List<Double> maxs = new ArrayList<>();
		private List<String> ids = new ArrayList<>();

		public Intervals() {

			// config styles
			barStyle.setStrokeWidth(0.2);
			barStyle.setStrokeColor(0x666666);
			barStyle.setFillColor(0x999999);
		}

		public void addInterval(double min, double max) {
			addInterval(min, max, null);
		}

		public void addInterval(double min, double max, String id) {

			// save the interval
			mins.add(min);
			maxs.add(max);
			ids.add(id);

			// update bounds
			xmax = mins.size()*(barSpacing + barWidth);
			ymax = Math.max(ymax, max);
		}

		public LeftVerticalAxis makeAxis() {
			SVGPlot.LeftVerticalAxis axis = new SVGPlot.LeftVerticalAxis();
			axis.min = 0.0;
			axis.max = Math.ceil(ymax/axisTicksOn)*axisTicksOn;
			axis.addTicksOn(axisTicksOn);
			return axis;
		}

		public void draw(SVG svg) {

			svg.putStyleClasses(barStyle);

			int n = mins.size();
			for (int i=0; i<n; i++) {

				double x = (i+1)*barSpacing + i*barWidth;
				double y1 = mins.get(i);
				double y2 = maxs.get(i);

				// if the interval is big enough, draw a rect
				SVG.Drawable d;
				if (y2 - y1 >= minRectHeight) {
					d = svg.makeRect(x, x + barWidth, y1, y2);
				} else {
					d = svg.makeLine(
						x, y1,
						x + barWidth, y1
					);
				}
				d.setStyleClass(barStyle)
					.setId(ids.get(i))
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
}
