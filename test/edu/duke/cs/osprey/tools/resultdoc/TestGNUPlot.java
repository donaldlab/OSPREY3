package edu.duke.cs.osprey.tools.resultdoc;


import edu.duke.cs.osprey.tools.FileTools;
import org.junit.Test;

import java.io.File;
import java.util.Arrays;


public class TestGNUPlot {

	private static void writePlot(Plot plot, String filename) {

		// no idea how to automate test checking here...
		// I just look at the SVG file manually
		FileTools.writeFile(plot.renderSvg(), new File(filename + ".svg"));
		FileTools.writeFileBytes(plot.renderPng(), new File(filename + ".png"));
	}

	@Test
	public void pointsX() {

		Plot plot = new Plot();
		plot.xlabel = "Value";
		plot.ylabels = Arrays.asList("A", "B", "C", "D");

		Plot.PointsX data1 = plot.new PointsX();
		data1.name = "Points 1";
		data1.data = Arrays.asList(0.0, 1.0, 5.0, 6.2);
		data1.size = 2.0;

		Plot.PointsX data2 = plot.new PointsX();
		data2.name = "Points 2";
		data2.data = Arrays.asList(3.3, 4.6, 9.3, 1.5);
		data2.size = 2.0;

		writePlot(plot, "pointsx");
	}

	@Test
	public void linesX() {

		Plot plot = new Plot();
		plot.xlabel = "Value";
		plot.ylabels = Arrays.asList("A", "B", "C", "D");

		Plot.LinesX data1 = plot.new LinesX();
		data1.name = "Lines 1";
		data1.data = Arrays.asList(0.0, 1.0, 5.0, 6.2);

		Plot.LinesX data2 = plot.new LinesX();
		data2.name = "Lines 2";
		data2.data = Arrays.asList(3.3, 4.6, 9.3, 1.5);

		writePlot(plot, "linesx");
	}

	@Test
	public void intervalsX() {

		Plot plot = new Plot();
		plot.xlabel = "Value";
		plot.ylabels = Arrays.asList("A", "B", "C");

		Plot.IntervalsX data = plot.new IntervalsX();
		data.name = "Intervals";
		data.data = Arrays.asList(
			new Plot.Interval(1.0, 2.5),
			new Plot.Interval(5.3, 7.2),
			new Plot.Interval(2.0, 4.0)
		);

		writePlot(plot, "intervalsx");
	}
}
